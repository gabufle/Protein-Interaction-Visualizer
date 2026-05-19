# Protein Interaction Visualizer & Predictor
# =============================================
# A reproducible pipeline for fetching, analyzing, visualizing, and
# predicting protein-protein interactions from real biological databases.
#
# Project Structure:
#   src/                    # Source code (modular, testable)
#     data_acquisition/    # Fetch PPI data from STRING, BioGRID, etc.
#     network_analysis/    # Graph metrics, topology, community detection
#     visualization/       # Interactive Plotly graphs, static plots
#     prediction/          # ML models for interaction prediction
#     utils/              # Config handling, logging, helpers
#   tests/                # Unit tests (pytest)
#   notebooks/            # Jupyter notebooks (exploration, demos)
#   data/
#     raw/                # Raw downloads (never modify)
#     processed/          # Cleaned, merged, versioned datasets
#   results/              # Output figures and prediction results
#   models/               # Trained ML model checkpoints
#
# Usage:
#   1. Install: pip install -r requirements.txt
#   2. Fetch data: python -m src.data_acquisition.string_fetcher
#   3. Visualize: python -m src.visualization.network_viz
#   4. Predict:  python -m src.prediction.interaction_predictor
#
# Reproducibility:
#   - All dependencies pinned in requirements.txt
#   - Random seeds fixed (random_state=42)
#   - Data fetched from public APIs with versioned queries
#   - Results and figures saved with timestamps
