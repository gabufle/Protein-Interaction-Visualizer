"""Data acquisition module: fetch protein interactions from various databases.

Available fetchers:
    - StringFetcher: STRING database (comprehensive, with confidence scores)
    - BioGRIDFetcher: BioGRID database (experimentally validated only)
    - (Future: IntAct, DIP, MINT, HINT, etc.)

Usage:
    from src.data_acquisition import StringFetcher, merge_string_biogrid
"""

from .string_fetcher import StringFetcher, fetch_human_ppi
from .biogrid_fetcher import BioGRIDFetcher, merge_string_biogrid

__all__ = [
    "StringFetcher",
    "fetch_human_ppi",
    "BioGRIDFetcher",
    "merge_string_biogrid",
]
