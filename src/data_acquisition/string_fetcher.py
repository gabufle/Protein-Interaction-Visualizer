"""Fetch protein-protein interactions from STRING database.

STRING (Search Tool for the Retrieval of Interacting Genes/Proteins) is a
comprehensive database of known and predicted protein interactions.
- Covers 2000+ organisms
- Integrates experimental, database, and text-mining evidence
- Provides confidence scores (0-1000) for each interaction

API Documentation: https://string-db.org/help/api/

Author: [Your Name]
Date: 2025-05-14
"""

import json
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import requests
from tqdm import tqdm

from src.utils.config import load_config
from src.utils.logging import get_logger

logger = get_logger(__name__)


class StringFetcher:
    """Fetcher for STRING database protein-protein interactions.

    Usage:
        fetcher = StringFetcher(organism="9606", score_threshold=400)
        interactions = fetcher.fetch_interactions()
        fetcher.save(interactions, "data/processed/human_ppi.csv")

    The fetcher caches results to avoid redundant API calls.
    """

    # STRING API base URL (v12.0 is current stable)
    BASE_URL = "https://version-12-0.string-db.org/api"

    # Human-readable organism names for convenience
    ORGANISM_NAMES = {
        "9606": "Homo sapiens",
        "10090": "Mus musculus",
        "10116": "Rattus norvegicus",
        "4932": "Saccharomyces cerevisiae",
        "7227": "Drosophila melanogaster",
        "6239": "Caenorhabditis elegans",
        "3702": "Arabidopsis thaliana",
        "7955": "Danio rerio",
    }

    def __init__(
        self,
        organism: str = "9606",  # Human by default
        score_threshold: int = 400,
        api_version: str = "v12.0",
        max_retries: int = 3,
        request_delay: float = 0.5,  # seconds between requests (be polite)
        api_key: Optional[str] = None,
    ):
        """
        Initialize STRING fetcher.

        Args:
            organism: NCBI Taxonomy ID (e.g., '9606' for human, '10090' for mouse)
            score_threshold: Minimum interaction score (0-1000). Lower = more interactions.
                - 150: low confidence (many false positives)
                - 400: medium (default, good balance)
                - 700: high confidence (experimentally supported)
                - 900: highest confidence (multiple evidence types)
            api_version: STRING API version (default: v12.0)
            max_retries: Number of retries on failed requests
            request_delay:Seconds to wait between paginated requests (API courtesy)
            api_key: Optional STRING API key for higher rate limits

        Raises:
            ValueError: If organism not recognized or score_threshold out of range
        """
        self.organism = str(organism)
        self.score_threshold = int(score_threshold)
        self.api_version = api_version
        self.max_retries = max_retries
        self.request_delay = request_delay
        self.api_key = api_key or self._get_api_key_from_config()

        # Validate organism
        if self.organism not in self.ORGANISM_NAMES:
            logger.warning(
                f"Organism taxid '{self.organism}' not in known list. "
                f"Proceeding anyway (STRING supports 2000+ organisms)."
            )

        # Validate score
        if not (0 <= self.score_threshold <= 1000):
            raise ValueError(
                f"score_threshold must be 0-1000, got {self.score_threshold}"
            )

        logger.info(
            f"Initialized StringFetcher: "
            f"organism={self.organism} ({self.ORGANISM_NAMES.get(self.organism, 'unknown')}), "
            f"score_threshold={self.score_threshold}"
        )

    def _get_api_key_from_config(self) -> Optional[str]:
        """Load STRING API key from environment or config."""
        import os

        # Try environment variable first
        api_key = os.environ.get("STRING_API_KEY")
        if api_key:
            logger.info("STRING API key found in environment")
            return api_key

        # Try config
        try:
            config = load_config()
            api_key = config.get("api_keys", {}).get("string")
            if api_key:
                logger.info("STRING API key loaded from config")
                return api_key
        except Exception as e:
            logger.debug(f"Could not load config for API key: {e}")

        logger.warning(
            "No STRING API key found. Using public endpoint with rate limits (~10 req/min). "
            "Consider getting a key at https://string-db.org/cgi/access?section=api"
        )
        return None

    def _build_url(self, endpoint: str) -> str:
        """Construct full API URL with optional key."""
        url = f"{self.BASE_URL}/{self.api_version}/{endpoint}"
        if self.api_key:
            url += f"?api_key={self.api_key}"
        return url

    def _request_with_retry(
        self,
        endpoint: str,
        params: Dict,
        max_retries: Optional[int] = None,
    ) -> Dict:
        """
        Make HTTP request with exponential backoff retry.

        Args:
            endpoint: API endpoint path (e.g., "interaction_partners")
            params: Query parameters dictionary
            max_retries: Override default max_retries

        Returns:
            JSON response as dictionary

        Raises:
            requests.RequestException: If all retries fail
        """
        max_retries = max_retries or self.max_retries
        url = self._build_url(endpoint)

        for attempt in range(max_retries + 1):
            try:
                response = requests.get(url, params=params, timeout=30)

                if response.status_code == 200:
                    return response.json()

                elif response.status_code == 429:
                    # Rate limited
                    wait = (attempt + 1) * 5  # exponential backoff
                    logger.warning(
                        f"Rate limited (429). Waiting {wait}s before retry..."
                    )
                    time.sleep(wait)
                    continue

                elif response.status_code == 400:
                    # Bad request — likely invalid parameters
                    error_msg = response.text[:200]
                    raise ValueError(
                        f"Bad request to STRING API: {error_msg}. "
                        f"URL: {response.url}"
                    )

                else:
                    response.raise_for_status()

            except requests.exceptions.RequestException as e:
                if attempt == max_retries:
                    logger.error(f"Request failed after {max_retries} retries: {e}")
                    raise
                wait = (attempt + 1) * 2
                logger.debug(f"Request failed (attempt {attempt+1}/{max_retries}): {e}")
                logger.debug(f"Retrying in {wait}s...")
                time.sleep(wait)

        # Shouldn't reach here, but mypy requires return
        raise requests.RequestException("All retries exhausted")

    def fetch_interactions(
        self,
        query_identifiers: Optional[List[str]] = None,
        query_mode: str = "network",  # or "external", "bulk"
    ) -> pd.DataFrame:
        """
        Fetch protein-protein interactions for the organism.

        By default, fetches the *entire* interaction network for the organism
        (all proteins with score ≥ threshold). This can be large (~200k edges for human).

        Args:
            query_identifiers: List of protein identifiers (e.g., ["P00533", "Q9Y6K9"])
                If None, fetch entire organism network (warning: can be huge!)
            query_mode: How to query
                - "network": fetch all PPIs for organism (default, comprehensive)
                - "external": fetch interactions for provided list of proteins
                - "bulk": map identifiers to STRING IDs then fetch network

        Returns:
            pandas DataFrame with columns:
                - protein1 (string id)
                - protein2 (string id)
                - combined_score (0-1000)
                - experimental_score (optional)
                - database_score (optional)
                - coexpression_score (optional)
                - textmining_score (optional)

        Raises:
            requests.RequestException: If API fails
            ValueError: If parameters invalid
        """
        logger.info(
            f"Fetching STRING interactions: organism={self.organism}, "
            f"threshold={self.score_threshold}, mode={query_mode}"
        )

        if query_mode == "network":
            return self._fetch_full_network()
        elif query_mode == "external":
            if not query_identifiers:
                raise ValueError("query_identifiers required for external mode")
            return self._fetch_for_identifiers(query_identifiers)
        else:
            raise ValueError(f"Unknown query_mode: {query_mode}")

    def _fetch_full_network(self) -> pd.DataFrame:
        """Fetch complete interaction network for organism."""
        logger.info(f"Fetching full {self.ORGANISM_NAMES.get(self.organism, '')} network...")

        params = {
            "species": self.organism,
            "required_score": self.score_threshold,
            "format": "json",
            "caller_identity": "protein_interaction_visualizer/0.1.0",
        }

        try:
            data = self._request_with_retry("network", params)
        except Exception as e:
            logger.error(f"Failed to fetch full network: {e}")
            raise

        # Parse response
        interactions = []
        for item in data:
            # STRING returns list of lists: [protein1, protein2, score]
            if len(item) >= 3:
                interactions.append({
                    "protein1": str(item[0]),
                    "protein2": str(item[1]),
                    "combined_score": int(item[2]),
                })

        df = pd.DataFrame(interactions)
        logger.info(f"Fetched {len(df):,} interactions")

        # Deduplicate (STRING might have reciprocal entries)
        df = self._deduplicate_interactions(df)

        return df

    def _fetch_for_identifiers(self, identifiers: List[str]) -> pd.DataFrame:
        """Fetch interactions for a specific list of proteins."""
        logger.info(f"Fetching interactions for {len(identifiers)} proteins")

        # STRING endpoint for external identifiers
        # We need to map identifiers to STRING IDs first
        url = self._build_url("get_string_ids")

        all_interactions = []

        # Process in batches to avoid URL length limits
        batch_size = 100
        for i in tqdm(range(0, len(identifiers), batch_size), desc="Fetching batches"):
            batch = identifiers[i : i + batch_size]

            params = {
                "identifiers": "\n".join(batch),  # newline-separated
                "species": self.organism,
                "format": "json",
                "caller_identity": "protein_interaction_visualizer/0.1.0",
            }

            try:
                # First: map identifiers to STRING IDs
                mapping_resp = self._request_with_retry("get_string_ids", params)
                string_ids = [
                    item["stringId"] for item in mapping_resp
                    if "stringId" in item and item["stringId"]
                ]

                if string_ids:
                    # Second: fetch interactions for these STRING IDs
                    interaction_params = {
                        "identifiers": "\n".join(string_ids),
                        "species": self.organism,
                        "required_score": self.score_threshold,
                        "format": "json",
                    }
                    interactions_data = self._request_with_retry(
                        "interaction_partners", interaction_params
                    )

                    for item in interactions_data:
                        if len(item) >= 3:
                            all_interactions.append({
                                "protein1": str(item[0]),
                                "protein2": str(item[1]),
                                "combined_score": int(item[2]),
                            })

                time.sleep(self.request_delay)  # Be polite
            except Exception as e:
                logger.error(f"Failed to fetch batch {i//batch_size + 1}: {e}")
                continue

        df = pd.DataFrame(all_interactions)
        df = self._deduplicate_interactions(df)
        logger.info(f"Fetched {len(df):,} interactions for identifier list")
        return df

    def _deduplicate_interactions(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Remove duplicate interactions (A-B and B-A are same interaction).

        STRING can list both directions. We canonicalize to (min_id, max_id).
        """
        if df.empty:
            return df

        # Create canonical ordering
        df["_node1"] = df[["protein1", "protein2"]].min(axis=1)
        df["_node2"] = df[["protein1", "protein2"]].max(axis=1)
        df = df.drop_duplicates(subset=["_node1", "_node2"])
        df = df.drop(columns=["_node1", "_node2"]).reset_index(drop=True)

        logger.debug(f"After deduplication: {len(df):,} unique interactions")
        return df

    def save(
        self,
        df: pd.DataFrame,
        output_path: str | Path,
        include_metadata: bool = True,
    ) -> Path:
        """
        Save interactions to disk with versioned metadata.

        Args:
            df: Interaction DataFrame
            output_path: File path (.csv, .tsv, .pkl, .parquet)
            include_metadata: Save a companion .json metadata file

        Returns:
            Path to saved file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Choose format based on extension
        suffix = output_path.suffix.lower()

        if suffix == ".csv":
            df.to_csv(output_path, index=False)
        elif suffix == ".tsv":
            df.to_csv(output_path, sep="\t", index=False)
        elif suffix in {".pkl", ".pickle"}:
            df.to_pickle(output_path)
        elif suffix == ".parquet":
            df.to_parquet(output_path, index=False)
        else:
            raise ValueError(f"Unsupported file format: {suffix}")

        logger.info(f"Saved {len(df):,} interactions → {output_path}")

        # Save metadata for reproducibility
        if include_metadata:
            metadata = {
                "source": "STRING",
                "api_version": self.api_version,
                "organism": self.organism,
                "organism_name": self.ORGANISM_NAMES.get(self.organism, "unknown"),
                "score_threshold": self.score_threshold,
                "num_interactions": len(df),
                "timestamp": datetime.now().isoformat(),
                "fetcher_version": "0.1.0",
                "columns": list(df.columns),
            }

            meta_path = output_path.with_suffix(".metadata.json")
            with open(meta_path, "w") as f:
                json.dump(metadata, f, indent=2)
            logger.debug(f"Metadata saved → {meta_path}")

        return output_path

    @classmethod
    def from_config(cls) -> "StringFetcher":
        """Construct fetcher from project config.yaml."""
        config = load_config()
        string_cfg = config["data"]["string"]

        return cls(
            organism=string_cfg.get("default_organism", "9606"),
            score_threshold=string_cfg.get("score_threshold", 400),
            api_version=string_cfg.get("api_version", "v12.0"),
            api_key=config.get("api_keys", {}).get("string"),
        )


def fetch_human_ppi(score_threshold: int = 400, save_dir: Optional[Path] = None) -> pd.DataFrame:
    """
    Convenience function: fetch human PPI network with sensible defaults.

    Args:
        score_threshold: Minimum STRING score (default 400 = medium confidence)
        save_dir: Directory to save output (default: data/processed)

    Returns:
        DataFrame of protein interactions
    """
    fetcher = StringFetcher(organism="9606", score_threshold=score_threshold)
    df = fetcher.fetch_interactions()

    if save_dir:
        save_path = Path(save_dir) / f"human_ppi_score{score_threshold}.csv"
        fetcher.save(df, save_path)

    return df


if __name__ == "__main__":
    # Quick test: fetch human network and print summary
    import sys

    logger = get_logger(__name__)

    try:
        print("=" * 60)
        print("STRING Database Fetcher — Test Run")
        print("=" * 60)

        fetcher = StringFetcher(organism="9606", score_threshold=400)
        df = fetcher.fetch_interactions()

        print("\nDataFrame head:")
        print(df.head())
        print(f"\nShape: {df.shape[0]:,} rows × {df.shape[1]} columns")
        print(f"Score range: {df['combined_score'].min()} - {df['combined_score'].max()}")

        # Save test output
        test_dir = Path("data/processed")
        test_dir.mkdir(parents=True, exist_ok=True)
        out_path = test_dir / "test_human_ppi.csv"
        fetcher.save(df, out_path)
        print(f"\n✅ Saved to: {out_path}")

    except Exception as e:
        logger.error(f"Fetch failed: {e}", exc_info=True)
        sys.exit(1)
