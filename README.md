```text

.
├── src/                   # Source code (modular, testable)
│   ├── data_acquisition/  # Fetch PPI data from STRING, BioGRID, etc.
│   ├── network_analysis/  # Graph metrics, topology, community detection
│   ├── visualization/     # Interactive Plotly graphs, static plots
│   ├── prediction/        # ML models for interaction prediction
│   └── utils/             # Config handling, logging, helpers
├── tests/                 # Unit tests (pytest)
├── notebooks/             # Jupyter notebooks (exploration, demos)
├── data/
│   ├── raw/               # Raw downloads — never modify
│   └── processed/         # Cleaned, merged, versioned datasets
├── results/               # Output figures and prediction results
└── models/                # Trained ML model checkpoints
