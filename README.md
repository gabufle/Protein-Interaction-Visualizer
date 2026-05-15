<div align="center">

# 🧬 Protein Interaction Visualizer & Predictor

**A reproducible pipeline for fetching, analyzing, visualizing, and predicting protein–protein interactions from real biological databases.**

[![Python](https://img.shields.io/badge/Python-3.9%2B-3776AB?style=flat-square&logo=python&logoColor=white)](https://python.org)
[![Plotly](https://img.shields.io/badge/Plotly-Interactive-3F4F75?style=flat-square&logo=plotly&logoColor=white)](https://plotly.com)
[![License](https://img.shields.io/badge/License-MIT-22c55e?style=flat-square)](LICENSE)
[![Tests](https://img.shields.io/badge/Tests-pytest-f59e0b?style=flat-square&logo=pytest&logoColor=white)](tests/)
[![Notebooks](https://img.shields.io/badge/Notebooks-Jupyter-F37626?style=flat-square&logo=jupyter&logoColor=white)](notebooks/)

</div>

---

## ✨ Overview

This project provides an end-to-end pipeline for working with protein–protein interaction (PPI) networks. From raw data acquisition to interactive visualization and machine-learning–based interaction prediction, everything is modular, testable, and reproducible.

**Data sources include:** STRING · BioGRID · and more

---

## 📁 Project Structure

```text
.
├── src/                        # Source code — modular & testable
│   ├── data_acquisition/       # 🌐  Fetch PPI data from STRING, BioGRID, etc.
│   ├── network_analysis/       # 📊  Graph metrics, topology, community detection
│   ├── visualization/          # 🎨  Interactive Plotly graphs & static plots
│   ├── prediction/             # 🤖  ML models for interaction prediction
│   └── utils/                  # 🔧  Config handling, logging, helpers
│
├── tests/                      # ✅  Unit tests (pytest)
├── notebooks/                  # 📓  Jupyter notebooks — exploration & demos
│
├── data/
│   ├── raw/                    # 🗃️  Raw downloads  (never modify)
│   └── processed/              # 🧹  Cleaned, merged, versioned datasets
│
├── results/                    # 📈  Output figures & prediction results
└── models/                     # 💾  Trained ML model checkpoints
```

---

## 🚀 Getting Started

```bash
# Clone the repository
git clone https://github.com/your-username/protein-interaction-visualizer.git
cd protein-interaction-visualizer

# Install dependencies
pip install -r requirements.txt

# Run the pipeline
python src/main.py
```

---

## 🧩 Modules at a Glance

| Module | Description |
|---|---|
| `data_acquisition` | Pulls raw PPI data from STRING, BioGRID, and other databases |
| `network_analysis` | Computes graph metrics, topology analysis, and community detection |
| `visualization` | Generates interactive Plotly graphs and publication-ready static plots |
| `prediction` | Trains and evaluates ML models to predict novel protein interactions |
| `utils` | Shared config management, logging utilities, and helper functions |

---

## 🧪 Running Tests

```bash
pytest tests/
```

---

## 📓 Notebooks

Explore the `notebooks/` directory for step-by-step demos and exploratory analyses. Great starting point if you're new to the project!

---

## 📌 Data Conventions

- **`data/raw/`** — Raw downloads. **Never modify these files.**
- **`data/processed/`** — Cleaned, merged, and versioned datasets ready for analysis.

---

<div align="center">

Made with 🧬 and curiosity

</div>