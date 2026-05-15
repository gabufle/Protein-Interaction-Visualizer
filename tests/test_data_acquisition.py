"""Tests for src.data_acquisition modules.

Tests:
    - StringFetcher: initialization, parameter validation, API mock
    - BioGRIDFetcher: similar
    - merge_string_biogrid: merging logic

Run with: pytest tests/test_data_acquisition.py -v
"""

import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from src.data_acquisition import StringFetcher, BioGRIDFetcher, merge_string_biogrid
from src.utils.config import load_config


@pytest.fixture
def sample_string_response():
    """Mock STRING API response for a small network."""
    return [
        ["P00533", "P04637", 950],
        ["P00533", "Q9Y6K9", 720],
        ["P04637", "Q9Y6K9", 400],
        ["P04637", "O14920", 300],
    ]


class TestStringFetcher:
    """Tests for StringFetcher class."""

    def test_initialization_valid(self):
        """Fetcher initializes with valid parameters."""
        fetcher = StringFetcher(organism="9606", score_threshold=400)
        assert fetcher.organism == "9606"
        assert fetcher.score_threshold == 400

    def test_initialization_invalid_score(self):
        """Score threshold must be 0-1000."""
        with pytest.raises(ValueError):
            StringFetcher(score_threshold=1100)
        with pytest.raises(ValueError):
            StringFetcher(score_threshold=-10)

    def test_unknown_organism_warning(self, caplog):
        """Non-standard organism IDs log warning but don't fail."""
        fetcher = StringFetcher(organism="99999")
        assert fetcher.organism == "99999"
        # caplog can check warnings if logger configured

    @patch("requests.get")
    def test_fetch_full_network(self, mock_get, sample_string_response):
        """Fetch full network successfully with mocked API."""
        # Configure mock
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = sample_string_response
        mock_get.return_value = mock_resp

        fetcher = StringFetcher(organism="9606", score_threshold=400)
        df = fetcher.fetch_interactions()

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 4
        assert "protein1" in df.columns
        assert "protein2" in df.columns
        assert "combined_score" in df.columns

    def test_deduplication(self):
        """Duplicate interactions are removed (A-B and B-A considered same)."""
        df = pd.DataFrame({
            "protein1": ["A", "B", "A"],
            "protein2": ["B", "A", "C"],
            "combined_score": [900, 900, 800],
        })
        fetcher = StringFetcher()
        cleaned = fetcher._deduplicate_interactions(df)
        # (A,B) appears twice → keep 1; (A,C) remains
        assert len(cleaned) == 2

    def test_save_to_csv(self, tmp_path):
        """DataFrame saves to CSV with metadata."""
        df = pd.DataFrame({
            "protein1": ["P1", "P2"],
            "protein2": ["P2", "P3"],
            "combined_score": [900, 800],
        })
        fetcher = StringFetcher()
        out = tmp_path / "test.csv"
        fetcher.save(df, out, include_metadata=True)

        assert out.exists()
        saved = pd.read_csv(out)
        assert len(saved) == 2
        # Metadata file created
        assert out.with_suffix(".metadata.json").exists()

    def test_save_invalid_format(self, tmp_path):
        """Unsupported file format raises error."""
        df = pd.DataFrame({"a": [1]})
        fetcher = StringFetcher()
        with pytest.raises(ValueError):
            fetcher.save(df, tmp_path / "bad.xyz")


class TestBioGRIDFetcher:
    """Tests for BioGRIDFetcher class."""

    def test_initialization(self):
        """Fetcher initializes correctly."""
        fetcher = BioGRIDFetcher(organism=9606)
        assert fetcher.organism == 9606

    def test_standardize_columns(self):
        """Column standardization produces expected schema."""
        raw_df = pd.DataFrame({
            "bioGRID_interaction_id": [1, 2],
            "SystematicNameA": ["BRCA1", "TP53"],
            "SystematicNameB": ["TP53", "MDM2"],
            "ExperimentalSystem": ["Two-Hybrid", "Co-IP"],
            "ExperimentalSystemType": ["physical", "physical"],
            "combined_score": [1000, 1000],
        })
        fetcher = BioGRIDFetcher()
        std = fetcher._standardize_columns(raw_df)

        assert "protein1" in std.columns
        assert "protein2" in std.columns
        assert std["combined_score"].iloc[0] == 1000
        assert std["source"].iloc[0] == "BioGRID"

    def test_deduplication(self):
        """Duplicate reciprocal entries removed."""
        df = pd.DataFrame({
            "bioGRID_interaction_id": [1, 2],
            "SystematicNameA": ["A", "B"],
            "SystematicNameB": ["B", "A"],
        })
        fetcher = BioGRIDFetcher()
        cleaned = fetcher._deduplicate(df)
        assert len(cleaned) == 1


class TestMergeDatasets:
    """Tests for merging STRING and BioGRID data."""

    def test_merge_non_overlapping(self):
        """Merges two completely disjoint networks."""
        string_df = pd.DataFrame({
            "protein1": ["A", "B"],
            "protein2": ["B", "C"],
            "combined_score": [900, 800],
        })
        biogrid_df = pd.DataFrame({
            "protein1": ["X", "Y"],
            "protein2": ["Y", "Z"],
            "combined_score": [1000, 1000],
        })

        merged = merge_string_biogrid(string_df, biogrid_df, prioritize="biogrid")

        assert len(merged) == 4
        sources = set(merged["source"])
        assert sources == {"STRING", "BioGRID"}

    def test_merge_with_overlap(self):
        """Overlapping edge gets combined source='both'."""
        string_df = pd.DataFrame({
            "protein1": ["A", "B"],
            "protein2": ["B", "C"],
            "combined_score": [900, 800],
        })
        biogrid_df = pd.DataFrame({
            "protein1": ["A", "B"],  # (A,B) overlaps with STRING
            "protein2": ["B", "D"],
            "combined_score": [1000, 1000],
        })

        merged = merge_string_biogrid(string_df, biogrid_df, prioritize="biogrid")
        overlap_row = merged[merged["protein1"] == "A"]
        assert len(overlap_row) == 1
        assert overlap_row["source"].iloc[0] == "both"
        # BioGRID confidence (1000) should be used
        assert overlap_row["combined_score"].iloc[0] == 1000


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
