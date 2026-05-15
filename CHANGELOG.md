# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased] — In Development

### Added
- Initial project scaffold (packaging, config, logging)
- StringFetcher: STRING API v12.0 client with deduplication, retry
- BioGRIDFetcher: fetch physical interactions with BioGRID API
- GraphAnalyzer: comprehensive network metrics, community detection (Louvain), hubs
- NetworkVisualizer: interactive Plotly 2D/3D HTML graphs, degree distribution plots
- InteractionPredictor: ML pipeline (RF, feature engineering, evaluation)
- `run_pipeline.py`: end-to-end orchestration script
- `demo.py`: synthetic data proof-of-concept run
- Dockerfile for containerized environment
- Unit tests for data_acquisition, network_analysis, visualization modules
- Documentation: README, PROJECT_PLAN, QUICKSTART, CONTRIBUTING

### Changed
- N/A (initial release)

### Deprecated
- N/A

### Removed
- N/A

### Fixed
- Fixed missing `import os` in logging module
- Fixed Plotly palette reference for newer versions
- Fixed indentation in PNG export blocks
- Ensured feature vector padding for skipped features in ML

### Security
- No hardcoded API keys; uses env vars

---

## [0.1.0] — 2026-05-14 (Proof-of-Concept)

Initial public release. This is a work-in-progress proof-of-concept demonstrating the full pipeline.

**Features:**
- Data acquisition from STRING & BioGRID
- Network analysis (topology, communities, hubs)
- Interactive 3D visualization (Plotly HTML)
- Machine learning predictor (Random Forest)
- Reproducible environment (Docker, pinned dependencies)
- Demo script and Jupyter notebook for exploration

**Known Limitations:**
- STRING API key required for real data (synthetic demo works without)
- PNG export requires Kaleido (optional, falls back to HTML only)
- Some advanced ML features (biological features) are placeholders
- Community detection only runs on largest component
- No CLI yet (Python API only)

**Next Steps (v0.2.0):**
- CLI interface (`protein-vis` command)
- Additional database: IntAct, DIP
- GO/domain-based features
- Graph Neural Network predictor
- Streamlit dashboard
- CI/CD (GitHub Actions)

---

## [Planned] v0.2.0 — TBD

### Added
- Command-line interface (argparse)
- IntAct fetcher
- GO term similarity features (via QuickGO API)
- Node2Vec embeddings (skipgram on graph)
- Simple Streamlit web app

### Changed
- Default model: XGBoost (better performance)
- Switch from CSV to Parquet for faster I/O

### Fixed
- Edge cases on disconnected graphs
- Memory optimization for >50k nodes

---

## [Planned] v1.0.0 — TBD

### Added
- Docker Compose for full stack
- Neo4j export
- REST API (FastAPI)
- Multi-organism comparison tool
- Publication-ready figure export (SVG, PDF)

### Changed
- Move to pyproject.toml-only (remove setup.py)
- Adopt Hydra for configuration

---

**Note:** This changelog follows the principle that any change visible to users (API, results, files) gets an entry. Internal-only refactors that don't affect output are omitted.
