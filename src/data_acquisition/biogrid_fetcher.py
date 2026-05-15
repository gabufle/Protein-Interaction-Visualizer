"""Fetch protein-protein interactions from BioGRID database.

BioGRID (Biological General Repository for Interaction Datasets) is a curated
database of physical and genetic interactions. All data are experimentally validated.

Key differences from STRING:
  - Only experimentally validated interactions (no predictions)
  - Broader organism coverage (all species)
  - Detailed experimental evidence types
  - No confidence scores (binary: interaction exists/doesn't)
  - Includes genetic interactions (not just physical PPI)

API: https://wiki.thebiogrid.org/doku.php/biogridrest

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


class BioGRIDFetcher:
    """Fetcher for BioGRID interaction database.

    Fetches experimentally validated protein-protein interactions.
    No API key required for basic access (but rate-limited).
    Registration gives higher limits.

    Example:
        fetcher = BioGRIDFetcher(organism=9606)  # Human
        df = fetcher.fetch_interactions()
        fetcher.save(df, "data/processed/biogrid_human.csv")
    """

    # BioGRID API base URL (v3.5+)
    BASE_URL = "https://webservice.thebiogrid.org"

    # Evidence types we're interested in (experimental physical interactions)
    EXPERIMENTAL_EVIDENCE = [
        "Experimental",
        "Co-Expression",
        "Co-Purification",
        "Co-Crystal Structure",
        "FRET",
        "PCA",
        "Two-Hybrid",
        "Co-Immunoprecipitation",
        "Pull Down",
        "Yeast Two-Hybrid",
        "Protein F complementation",
    ]

    def __init__(
        self,
        organism: Optional[int] = None,  # NCBI TaxID, None = all organisms
        evidence_types: Optional[List[str]] = None,
        search_names: bool = True,  # Also search by gene names
        include_secondary: bool = True,  # Include secondary IDs
        max_retries: int = 3,
        request_delay: float = 0.5,
        username: Optional[str] = None,
        password: Optional[str] = None,
    ):
        """
        Initialize BioGRID fetcher.

        Args:
            organism: NCBI Taxonomy ID (e.g., 9606 for human, 10090 for mouse).
                If None, fetches across all organisms (can be massive).
            evidence_types: List of BioGRID experimental system types to include.
                Default: ["Experimental"] plus common physical methods.
            search_names: If True, search by gene/protein names (slower but more coverage)
            include_secondary: Include secondary gene IDs (aliases)
            max_retries: HTTP request retry count
            request_delay: Seconds between paginated requests
            username: BioGRID account username (optional; increases rate limit)
            password: BioGRID password (optional)
        """
        self.organism = organism
        self.evidence_types = evidence_types or ["Experimental"]
        self.search_names = search_names
        self.include_secondary = include_secondary
        self.max_retries = max_retries
        self.request_delay = request_delay

        # Authentication
        self.username = username or self._get_username_from_config()
        self.password = password or self._get_password_from_config()

        if self.username and self.password:
            logger.info("BioGRID: Using authenticated access")
            self.auth = (self.username, self.password)
        else:
            logger.info("BioGRID: Using anonymous access (lower rate limit)")
            self.auth = None

    def _get_username_from_config(self) -> Optional[str]:
        """Load BioGRID username from environment or config."""
        import os

        username = os.environ.get("BIOGRID_USERNAME")
        if username:
            return username

        try:
            config = load_config()
            return config.get("api_keys", {}).get("biogrid_username")
        except Exception:
            return None

    def _get_password_from_config(self) -> Optional[str]:
        """Load BioGRID password from environment."""
        import os

        password = os.environ.get("BIOGRID_PASSWORD")
        if password:
            return password

        try:
            config = load_config()
            return config.get("api_keys", {}).get("biogrid_password")
        except Exception:
            return None

    def _request_with_retry(
        self,
        endpoint: str,
        params: Dict,
        max_retries: Optional[int] = None,
    ) -> Dict:
        """Make HTTP request with retry logic."""
        max_retries = max_retries or self.max_retries
        url = f"{self.BASE_URL}/{endpoint}"

        # Always request JSON format
        params["format"] = "json"
        params["caller_identity"] = "protein_interaction_visualizer/0.1.0"

        for attempt in range(max_retries + 1):
            try:
                response = requests.get(
                    url,
                    params=params,
                    auth=self.auth,
                    timeout=30,
                )

                if response.status_code == 200:
                    return response.json()

                elif response.status_code == 429:
                    wait = (attempt + 1) * 5
                    logger.warning(f"BioGRID rate limited (429). Waiting {wait}s...")
                    time.sleep(wait)
                    continue

                elif response.status_code in {400, 404}:
                    error_msg = response.text[:300]
                    raise ValueError(
                        f"BioGRID API error {response.status_code}: {error_msg}"
                    )

                else:
                    response.raise_for_status()

            except requests.exceptions.RequestException as e:
                if attempt == max_retries:
                    logger.error(f"Request failed after {max_retries} retries: {e}")
                    raise
                wait = (attempt + 1) * 2
                logger.debug(f"Request failed (attempt {attempt+1}): {e}")
                time.sleep(wait)

        raise requests.RequestException("All retries exhausted")

    def fetch_interactions(
        self,
        gene_list: Optional[List[str]] = None,
        max_interactions: Optional[int] = None,
    ) -> pd.DataFrame:
        """
        Fetch experimentally validated interactions from BioGRID.

        Args:
            gene_list: Optional list of gene symbols or protein IDs to fetch.
                If None, fetches ALL interactions for the organism (warning: huge!).
                Recommended: provide a gene list for organism of interest.
            max_interactions: Optional cap on total interactions to fetch
                (for testing; real analyses should get all available)

        Returns:
            DataFrame with columns:
                - bioGRID_interaction_id (unique interaction identifier)
                - EntrezGeneID_A, EntrezGeneID_B (NCBI gene IDs)
                - SystematicNameA, SystematicNameB (gene symbols)
                - ExperimentalSystem (e.g., "Two-Hybrid")
                - ExperimentalSystemType (physical/genetic)
                - Author, PubMedID
                - Throughput (Low/High)
                - Taxid_A, Taxid_B (organism IDs)
                - ... plus qualifiers
        """
        logger.info(
            f"Fetching BioGRID interactions: organism={self.organism}, "
            f"evidence={self.evidence_types}"
        )

        params = {
            "searchNames": str(self.search_names).lower(),
            "includeInteractors": "true",  # include both directions
            "selfLoop": "false",  # exclude self-interactions
            "interSpeciesExcluded": "true",  # only same organism
        }

        if self.organism is not None:
            params["taxId"] = str(self.organism)

        if gene_list:
            # Gene list mode: fetch interactions for specific genes
            # Format: gene1|gene2|gene3...
            params["geneList"] = "|".join(gene_list[:1000])  # API limit

        all_interactions = []
        page = 1
        page_size = 10000  # Max allowed
        has_more = True

        with tqdm(desc="Fetching BioGRID pages", unit="page") as pbar:
            while has_more:
                params["page"] = str(page)
                params["pageSize"] = str(page_size)

                try:
                    data = self._request_with_retry("interactions", params)
                except Exception as e:
                    logger.error(f"Failed to fetch page {page}: {e}")
                    break

                # Parse interactions (list of dicts)
                if isinstance(data, dict) and "interactions" in data:
                    interactions = data["interactions"]
                elif isinstance(data, list):
                    interactions = data
                else:
                    logger.warning(f"Unexpected response format: {type(data)}")
                    break

                if not interactions:
                    logger.info("No more interactions found")
                    has_more = False
                    continue

                all_interactions.extend(interactions)
                logger.debug(f"Page {page}: {len(interactions):,} interactions")

                # Progress
                pbar.update(1)
                pbar.set_postfix({"total": f"{len(all_interactions):,}"})

                # Check pagination
                if len(interactions) < page_size:
                    has_more = False
                else:
                    page += 1
                    time.sleep(self.request_delay)

                # Optional max cap
                if max_interactions and len(all_interactions) >= max_interactions:
                    logger.info(f"Reached max_interactions limit: {max_interactions}")
                    all_interactions = all_interactions[:max_interactions]
                    break

        # Convert to DataFrame
        if not all_interactions:
            logger.warning("No interactions fetched!")
            return pd.DataFrame()

        df = pd.DataFrame(all_interactions)
        logger.info(f"Fetched {len(df):,} raw interaction records")

        # Filter to desired evidence types if specified
        if self.evidence_types:
            before = len(df)
            df = df[
                df["ExperimentalSystem"].isin(self.evidence_types)
            ].copy()
            logger.info(
                f"Filtered by evidence types {self.evidence_types}: "
                f"{before:,} → {len(df):,}"
            )

        # Deduplicate (BioGRID may have reciprocal entries)
        df = self._deduplicate(df)

        # Standardize column names for compatibility with STRING format
        df = self._standardize_columns(df)

        return df

    def _deduplicate(self, df: pd.DataFrame) -> pd.DataFrame:
        """Remove duplicate interactions (both directions)."""
        if df.empty:
            return df

        # Use BioGRID interaction ID or create canonical ordering
        if "bioGRID_interaction_id" in df.columns:
            # Keep the first occurrence (BioGRID assigns separate IDs to reciprocal)
            df = df.drop_duplicates(subset=["bioGRID_interaction_id"])
        else:
            # Create canonical ordering by gene symbol or Entrez ID
            col_a = "SystematicNameA" if "SystematicNameA" in df.columns else "EntrezGeneID_A"
            col_b = "SystematicNameB" if "SystematicNameB" in df.columns else "EntrezGeneID_B"

            df["_node1"] = df[[col_a, col_b]].min(axis=1)
            df["_node2"] = df[[col_a, col_b]].max(axis=1)
            df = df.drop_duplicates(subset=["_node1", "_node2"])
            df = df.drop(columns=["_node1", "_node2"])

        logger.debug(f"After deduplication: {len(df):,} interactions")
        return df

    def _standardize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Standardize column names to be compatible with STRING format.

        We create a unified schema:
            - protein1, protein2: protein identifiers (STRING IDs or gene names)
            - combined_score: confidence (1-1000 scale; for BioGRID we use 1000 = experimental)
            - source: "BioGRID"
            - experimental_system: BioGRID's ExperimentalSystem value
            - pubmed_id: PubMed reference
        """
        # Determine identifier columns
        if "STRING_id_A" in df.columns and "STRING_id_B" in df.columns:
            # Already has STRING IDs (BioGRID maps to STRING sometimes)
            df = df.rename(columns={
                "STRING_id_A": "protein1",
                "STRING_id_B": "protein2",
            })
        elif "SystematicNameA" in df.columns and "SystematicNameB" in df.columns:
            df = df.rename(columns={
                "SystematicNameA": "protein1",
                "SystematicNameB": "protein2",
            })
        else:
            # Fall back to Entrez IDs
            df = df.rename(columns={
                "EntrezGeneID_A": "protein1",
                "EntrezGeneID_B": "protein2",
            })

        # Confidence score: BioGRID doesn't provide scores.
        # Use 1000 (max) for experimentally validated interactions.
        df["combined_score"] = 1000  # Perfect confidence for validated

        # Add source column
        df["source"] = "BioGRID"

        # Keep important metadata
        keep_cols = [
            "protein1", "protein2", "combined_score",
            "ExperimentalSystem", "ExperimentalSystemType",
            "PubMedID", "source", "bioGRID_interaction_id",
        ]
        keep_cols = [c for c in keep_cols if c in df.columns]

        return df[keep_cols].reset_index(drop=True)

    def save(self, df: pd.DataFrame, output_path: str | Path) -> Path:
        """Save DataFrame to file (CSV/TSV/Pickle)."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

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

        logger.info(f"Saved {len(df):,} BioGRID interactions → {output_path}")

        # Metadata
        metadata = {
            "source": "BioGRID",
            "organism": self.organism,
            "evidence_types": self.evidence_types,
            "num_interactions": len(df),
            "timestamp": datetime.now().isoformat(),
            "fetcher_version": "0.1.0",
            "columns": list(df.columns),
        }

        meta_path = output_path.with_suffix(".metadata.json")
        with open(meta_path, "w") as f:
            json.dump(metadata, f, indent=2)

        return output_path


def merge_string_biogrid(
    string_df: pd.DataFrame,
    biogrid_df: pd.DataFrame,
    prioritize: str = "biogrid",  # or "string"
) -> pd.DataFrame:
    """
    Merge STRING and BioGRID interaction datasets into a unified network.

    Strategy:
        - Deduplicate interactions that appear in both databases
        - For overlaps, use confidence from the preferred source
        - Add a 'source' column indicating which database(s) reported it

    Args:
        string_df: DataFrame from StringFetcher
        biogrid_df: DataFrame from BioGRIDFetcher
        prioritize: Which source to trust for overlapping edges?
            "biogrid": BioGRID-confirmed edges override STRING predictions
            "string": Keep STRING confidence scores even if also in BioGRID

    Returns:
        Merged DataFrame with columns:
            protein1, protein2, combined_score, source (STRING|BioGRID|both)
    """
    logger.info(
        f"Merging STRING ({len(string_df):,}) and BioGRID ({len(biogrid_df):,})"
    )

    # Canonicalize both datasets (order proteins consistently)
    def canonicalize(df):
        df = df.copy()
        df["_node1"] = df[["protein1", "protein2"]].min(axis=1)
        df["_node2"] = df[["protein1", "protein2"]].max(axis=1)
        return df

    string_df = canonicalize(string_df)
    biogrid_df = canonicalize(biogrid_df)

    # Identify overlaps (edges in both)
    string_key = set(zip(string_df["_node1"], string_df["_node2"]))
    biogrid_key = set(zip(biogrid_df["_node1"], biogrid_df["_node2"]))
    overlap = string_key & biogrid_key

    logger.info(f"Overlap: {len(overlap):,} edges appear in both databases")

    # Build merged set
    merged_rows = []

    # Add STRING-only edges
    string_only = string_df[~string_df.set_index(["_node1", "_node2"]).index.isin(overlap)]
    for _, row in string_only.iterrows():
        merged_rows.append({
            "protein1": row["protein1"],
            "protein2": row["protein2"],
            "combined_score": row["combined_score"],
            "source": "STRING",
        })

    # Add BioGRID-only edges
    biogrid_only = biogrid_df[~biogrid_df.set_index(["_node1", "_node2"]).index.isin(overlap)]
    for _, row in biogrid_only.iterrows():
        merged_rows.append({
            "protein1": row["protein1"],
            "protein2": row["protein2"],
            "combined_score": row["combined_score"],
            "source": "BioGRID",
        })

    # Handle overlaps
    if prioritize == "biogrid":
        # Use BioGRID's confidence (1000 = validated) for overlapping edges
        overlap_df = biogrid_df[biogrid_df.set_index(["_node1", "_node2"]).index.isin(overlap)]
        for _, row in overlap_df.iterrows():
            merged_rows.append({
                "protein1": row["protein1"],
                "protein2": row["protein2"],
                "combined_score": 1000,  # BioGRID validated
                "source": "both",
            })
    else:  # prioritize == "string"
        overlap_df = string_df[string_df.set_index(["_node1", "_node2"]).index.isin(overlap)]
        for _, row in overlap_df.iterrows():
            merged_rows.append({
                "protein1": row["protein1"],
                "protein2": row["protein2"],
                "combined_score": row["combined_score"],
                "source": "both",
            })

    merged = pd.DataFrame(merged_rows)
    logger.info(f"Merged network: {len(merged):,} unique interactions")
    logger.info(f"  Source breakdown:\n{merged['source'].value_counts().to_string()}")

    return merged


if __name__ == "__main__":
    # Quick test fetch (small organism for speed)
    import sys
    from src.utils.logging import get_logger

    logger = get_logger(__name__)

    print("=" * 60)
    print("BioGRID Fetcher — Test (Yeast/S. cerevisiae)")
    print("=" * 60)

    try:
        # Yeast has smaller network (~10k interactions)
        fetcher = BioGRIDFetcher(organism=4932, evidence_types=["Experimental"])
        df = fetcher.fetch_interactions(max_interactions=1000)  # Limit for test

        print(f"\nFetched {len(df):,} interactions")
        print("\nColumns:", list(df.columns))
        print("\nSample:")
        print(df.head(3).to_string())

        # Save
        out_path = Path("data/raw") / "test_biogrid_yeast.csv"
        fetcher.save(df, out_path)
        print(f"\n✅ Saved to {out_path}")

    except Exception as e:
        logger.error(f"Test failed: {e}", exc_info=True)
        sys.exit(1)
