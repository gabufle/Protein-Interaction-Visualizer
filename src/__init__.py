"""Protein Interaction Visualizer & Predictor

A comprehensive toolkit for fetching, analyzing, visualizing, and predicting
protein-protein interactions from real-world biological databases.

Modules:
    data_acquisition: Fetch PPI data from STRING, BioGRID, etc.
    network_analysis  : Graph metrics, community detection, hub identification
    visualization     : Interactive Plotly network graphs
    prediction        : ML models for interaction prediction
    utils             : Configuration, logging, helper functions

For detailed documentation, see README.md in the project root.
"""

__version__ = "0.1.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

from .utils import load_config, get_logger

__all__ = ["load_config", "get_logger"]
