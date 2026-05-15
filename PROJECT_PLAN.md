# Protein Interaction Visualizer & Predictor — Project Plan

**Version:** 1.0  
**Date:** May 14, 2025  
**Author:** [Your Name]  
**Repository:** `protein-interaction-visualizer`  

---

## Executive Summary

This project delivers a production-grade, reproducible pipeline for **protein-protein interaction (PPI) analysis**. It fetches real interaction data from biological databases (STRING, BioGRID), analyzes network topology, creates interactive 3D visualizations, and trains machine-learning models to predict novel interactions.

**Target users:** Computational biologists, bioinformatics researchers, systems biologists studying protein networks.

**Key innovation:** End-to-end workflow in a single Python package — from raw data to published figures — with full reproducibility guarantees.

---

## Project Vision & Objectives

### Vision
Democratize PPI analysis by providing an open-source, well-documented, one-command pipeline that works on standard laptops yet scales to large networks.

### Objectives (v1.0)
1. **Data Acquisition:** Reliable fetching from STRING and BioGRID with rate-limit handling
2. **Network Analysis:** Comprehensive topological metrics, community detection, hub identification
3. **Visualization:** Browser-based interactive 3D graphs (Plotly) and static publication figures
4. **Prediction:** ML classifier (Random Forest) to infer missing interactions
5. **Reproducibility:** Pinned dependencies, Docker, detailed logging, and versioned outputs
6. **Extensibility:** Modular architecture so researchers can add new databases, features, or models

### Success Metrics
- [ ] Complete pipeline runs from raw data to model in < 5 minutes on a laptop
- [ ] Model achieves ROC AUC > 0.85 on held-out test set (on synthetic validation networks)
- [ ] Interactive HTML visualizations include zoom, pan, hover tooltips
- [ ] All code passes unit tests (≥80% coverage)
- [ ] Single-command install and run for new users
- [ ] Documentation includes Jupyter notebook demo and API reference

---

## System Architecture

```
┌───────────────────────────────────────────────────────────────────┐
│                        User Interface                           │
│  (CLI · Jupyter Notebook · Python API · scheduled jobs)         │
└─────────────────────────────┬─────────────────────────────────────┘
                              │
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
┌───────▼────────┐   ┌────────▼───────┐   ┌────────▼────────┐
│  Data           │   │  Network       │   │  Machine        │
│  Acquisition    │   │  Analysis      │   │  Learning       │
│  (STRING,       │   │  (Graph        │   │  (Predictor)    │
│   BioGRID)      │   │   metrics,     │   │                 │
│                 │   │   communities) │   │                 │
└───────┬────────┘   └────────┬───────┘   └────────┬────────┘
        │                    │                     │
        └────────────────────┼─────────────────────┘
                             │
                   ┌─────────▼───────────┐
                   │   Visualization     │
                   │   (Plotly HTML,     │
                   │    PNG exports)     │
                   └─────────┬───────────┘
                             │
                   ┌─────────▼───────────┐
                   │   Output Storage    │
                   │   (CSV · JSON ·     │
                   │    HTML · Models)   │
                   └─────────────────────┘
```

**Modules:**
- `src/data_acquisition/` — API wrappers for STRING & BioGRID
- `src/network_analysis/` — NetworkX-based metric computation
- `src/visualization/` — Plotly interactive figures
- `src/prediction/` — scikit-learn ML pipeline
- `src/utils/` — Config, logging, reproducibility

---

## Implementation Phases

### Phase 1: Foundation & Data Layer (Week 1) ✅ DONE
**Goal:** Stable data fetching with retry logic and test coverage.

**Tasks:**
- [x] Design config schema (YAML + env overrides)
- [x] Implement `StringFetcher` with deduplication, session reuse, exponential backoff
- [x] Implement `BioGRIDFetcher` (REST API)
- [x] Standardize column schemas across sources
- [x] Merge strategy (prefer BioGRID confidence when overlapping)
- [x] Unit tests with mocked API responses
- [x] Documentation of API limits and usage quotas

**Milestone:** `python -c "from src.data_acquisition import fetch_human_ppi"` runs and saves real data.

---

### Phase 2: Network Analysis Engine (Week 1-2) ✅ DONE
**Goal:** Compute full suite of graph-theoretic metrics.

**Tasks:**
- [x] `GraphAnalyzer` class — wrapper around NetworkX
- [x] Metric suite: density, clustering, diameter, centrality (degree, betweenness, eigenvector), assortativity
- [x] Community detection: Louvain (primary), Girvan-Newman (hierarchical), label propagation (fast)
- [x] Hub ranking with multi-method support
- [x] JSON/CSV persistence of all results
- [x] Integration tests on benchmark graphs (Karate, Zachary)
- [x] Summary report generator

**Milestone:** `analyzer.compute_metrics()` returns complete JSON with < 30s runtime on 10k-node graph.

---

### Phase 3: Interactive Visualization (Week 2-3) ✅ DONE
**Goal:** High-quality, interactive 3D network figures for the web.

**Tasks:**
- [x] `NetworkVisualizer` class — Fruchterman-Reingold force-directed layout
- [x] Node sizing by degree, coloring by community
- [x] Edge opacity by confidence score
- [x] Hover tooltips with protein ID, degree, community
- [x] Export to standalone HTML (embedded JS, no server)
- [x] Static PNG export for publications (via Kaleido)
- [x] Supporting plots: degree distribution (log-log), centrality bar charts
- [x] 2D+3D layout options
- [x] Performance profiling (layout caching for large graphs)

**Milestone:** Generate interactive 3D graph of 200-node network openable in browser.

---

### Phase 4: Machine Learning Predictor (Week 3-4) ✅ DONE
**Goal:** Train a model to predict whether two proteins interact based on topological features.

**Tasks:**
- [x] Feature engineering:
  - Common neighbors, Jaccard, Adamic-Adar, Resource Allocation, Preferential Attachment
  - Community-aware features (soundarajan-hopcroft)
  - Placeholder for biological features (GO, domains)
- [x] Negative sampling strategy (non-edges, balanced classes)
- [x] Model Zoo: Random Forest (primary), Logistic Regression, Gradient Boosting
- [x] Evaluation: ROC AUC, Average Precision, confusion matrix
- [x] Hyperparameter tuning flag (but not automatic grid search — leave manual)
- [x] Model persistence (joblib), including StandardScaler
- [x] `predict()` method for arbitrary pairs
- [x] Test harness comparing feature importance across models

**Milestone:** Model achieves 85%+ accuracy on held-out test of synthetic network; feature importances biologically plausible.

---

### Phase 5: Integration, CLI & Reproducibility (Week 4-5) ✅ DONE
**Goal:** One-command pipeline and production-ready project hygiene.

**Tasks:**
- [x] `run_pipeline.py` — single command runs full workflow
- [x] Config-driven parameters (organism, score threshold, output dirs)
- [x] Timestamped output directories with symlink to `latest/`
- [x] Structured JSON logging to `logs/`
- [x] Random seed control for reproducibility
- [x] Dockerfile for containerized execution
- [x] `pytest` suite with > 80% coverage
- [x] Colored console output (colorlog)
- [x] Error handling and graceful degradation (missing optional deps)

**Milestone:** New user can run `python run_pipeline.py` and get complete results from scratch in < 10 min.

---

## Module Specifications

### 1. Data Acquisition Module (`src/data_acquisition/`)

```python
from src.data_acquisition import StringFetcher, BioGRIDFetcher

# STRING API
fetcher = StringFetcher(organism="9606", score_threshold=400)
df = fetcher.fetch_interactions()    # pandas DataFrame
fetcher.save(df, "output.csv")
```

**Features:**
- Automatic deduplication (AB ≡ BA)
- Progress bars via `tqdm`
- Rate-limit friendly (1 sec pause between requests)
- Session reuse (requests.Session)
- Retry with exponential backoff (5xx errors)
- Output formats: CSV (default), TSV, Parquet

**BioGRID fetcher:**
- Filters by organism NCBI taxid
- Optional evidence types: ["physical", "genetic"]
- Standardizes column names to match STRING

**Merge:** `merge_string_biogrid(df_str, df_bg, prioritize="biogrid")`

---

### 2. Network Analysis Module (`src/network_analysis/`)

```python
from src.network_analysis import GraphAnalyzer

analyzer = GraphAnalyzer(G)
metrics = analyzer.compute_metrics()
communities = analyzer.detect_communities(algorithm="louvain")
hubs = analyzer.find_hubs(top_n=20, method="betweenness")
analyzer.save_results(metrics, hubs, communities, "output/")
```

**Metrics computed:**
- Basic: nodes, edges, density, largest component ratio
- Topology: diameter, average shortest path, transitivity, avg clustering
- Centrality: degree, betweenness (approx when large), eigenvector (numpy), closeness (optional)
- Assortativity (degree correlation)
- k-core decomposition

**Community detection:** Louvain (fast, modularity-optimizing), Label Propagation (very fast), Girvan-Newman (hierarchical)

**Hub ranking:** Multi-criteria: degree°, betweenness, eigenvector, composite score

---

### 3. Visualization Module (`src/visualization/`)

```python
from src.visualization import NetworkVisualizer

viz = NetworkVisualizer(df_interactions, weight_col="combined_score")
fig = viz.plot_interactive_graph(
    title="Human PPI Network",
    output_path="results/figures/network.html",
    dim=3,
    show_labels=False,
    layout_algorithm="fr",
)
fig.show()
```

**Features:**
- 2D or 3D force-directed layouts
- Node size: degree (default) or custom column
- Node color: community (categorical) or metric (continuous colormap)
- Edge opacity/width: interaction confidence
- Hover: protein ID, degree, community
- Export: standalone HTML, PNG (high DPI via Kaleido)
- Degree distribution plot (checks power-law property)
- Centrality comparison bar charts

**Layouts:** Fruchterman-Reingold (spring), Kamada-Kawai, Circular, Spectral, Random

---

### 4. Prediction Module (`src/prediction/`)

```python
from src.prediction import InteractionPredictor

predictor = InteractionPredictor(G)
X, y, pairs = predictor.prepare_features(neg_ratio=1.0, test_size=0.2)
X_train, y_train = predictor.train_data
predictor.train(X_train, y_train, model_type="random_forest", n_estimators=200)
results = predictor.evaluate(predictor.test_data[0], predictor.test_data[1])
predictor.save_models("models/")

# Predict new pairs later
new_predictions = predictor.predict([("Q9Y6K9", "O14920")])
print(new_predictions[["protein_A","protein_B","probability","predicted"]])
```

**Features (topological):**
- `common_neighbors` — count of shared interactors
- `jaccard_coefficient` — Jaccard similarity of neighbor sets
- `adamic_adar` — favors shared neighbors with few connections
- `preferential_attachment` — product of degrees
- `resource_allocation` — inverse-degree resource sharing
- `cn_soundarajan_hopcroft` — common neighbors within same community (requires community labels)
- `within_inter_cluster` — placeholder for future graph embedding features

**Optional biological features (future):**
- GO term Jaccard similarity
- Protein domain overlap
- Sequence identity (BLAST)
- Co-expression correlation

**Models:**
- Random Forest (default, robust)
- Gradient Boosting (better accuracy, slower)
- Logistic Regression (fast baseline)

**Outputs:**
- Trained model (`.joblib`)
- StandardScaler (feature normalization)
- CSV of test predictions
- Feature importance plot
- ROC & PR curves (to be added)

---

## Code Quality & Testing

### Unit Tests (`tests/`)
- **data_acquisition**: API mocking, CSV parsing, deduplication logic
- **network_analysis**: metric correctness against known small graphs
- **visualization**: layout outputs, HTML file generation

Run: `pytest tests/ -v --cov=src`

### Type Hints
All functions annotated with `typing` hints. Run `mypy src/` for static analysis.

### Linting
- `black` formatting (line length 100)
- `isort` imports
- `flake8` style checks

---

## Reproducibility Strategy

1. **Pinned dependencies:** `requirements.txt` with exact versions
2. **Containerization:** `Dockerfile` (Python 3.10-slim) for environment parity
3. **Random seeds:** Configurable via `config.yaml` and `PROTEIN_VIS_RANDOM_SEED` env
4. **Deterministic algorithms:** NetworkX seeds set where possible (layout, community detection)
5. **Versioned data:** Output filenames include timestamps, plus `latest` symlinks
6. **Logging:** Detailed logs including Git hash (if repo), runtime, parameters
7. **Environment snapshot:** `pip freeze > requirements.lock.txt` tracks exact packages

**To reproduce on another machine:**
```bash
git clone <repo>
cd protein-interaction-visualizer
docker build -t ppi-vis .
docker run -v $(pwd)/data:/app/data ppi-vis
# Or conda env:
conda create -n ppi python=3.10
conda activate ppi
pip install -r requirements.txt
python run_pipeline.py --config config.yaml
```

---

## Usage Guide

### Quick Start (Synthetic Demo)
```bash
# Clone and install
git clone https://github.com/yourname/protein-interaction-visualizer.git
cd protein-interaction-visualizer
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
pip install -e .

# Run demo with synthetic data (no API key needed)
python demo.py

# Outputs appear in results/figures/ (HTML) and models/
```

### Real Data (STRING API)
```bash
# Get STRING API key (free) from https://string-db.org/cgi/access?section=api
export STRING_API_KEY="your_key"

# Edit config.yaml to set organism and score threshold
# Then run full pipeline:
python run_pipeline.py

# Or use Jupyter notebook:
jupyter lab notebooks/end_to_end_demo.ipynb
```

### Python API (in your own script)
```python
from src.data_acquisition import StringFetcher
from src.visualization import NetworkVisualizer
from src.prediction import InteractionPredictor

# 1. Fetch
fetcher = StringFetcher(organism="9606", score_threshold=700)
df = fetcher.fetch_interactions()

# 2. Visualize
viz = NetworkVisualizer(df)
viz.plot_interactive_graph(output_path="human_ppi.html", dim=3)

# 3. Predict (requires graph)
G = nx.from_pandas_edgelist(df, "protein1", "protein2")
predictor = InteractionPredictor(G)
predictor.prepare_features()
predictor.train(...)
```

### CLI (to be implemented)
Make `protein-vis` command with subcommands: `fetch`, `analyze`, `viz`, `predict`.

---

## Project Structure

```
protein-interaction-visualizer/
├── README.md                 # Project overview, quick start
├── PROJECT_PLAN.md           # This file
├── USAGE.md                  # Detailed usage examples
├── CONTRIBUTING.md           # Dev guidelines
├── CHANGELOG.md              # Version history
├── LICENSE                   # MIT
├── .env.example              # Template for API keys
├── .gitignore
├── requirements.txt          # Pinned core dependencies
├── pyproject.toml            # Build system config
├── setup.py                  # Editable install
├── config.yaml               # Default parameters
├── pytest.ini
├── run_pipeline.py           # End-to-end script
├── demo.py                   # Synthetic data demo
├── Dockerfile                # Container build
│
├── data/
│   ├── raw/                  # Raw downloads (optional)
│   └── processed/            # CSVs: networks, metrics, hubs
│
├── src/
│   ├── protein_interaction/
│   │   ├── data_acquisition/
│   │   │   ├── string_fetcher.py
│   │   │   ├── biogrid_fetcher.py
│   │   │   └── __init__.py
│   │   ├── network_analysis/
│   │   │   ├── graph_analyzer.py
│   │   │   └── __init__.py
│   │   ├── visualization/
│   │   │   ├── network_viz.py
│   │   │   └── __init__.py
│   │   ├── prediction/
│   │   │   ├── interaction_predictor.py
│   │   │   └── __init__.py
│   │   └── utils/
│   │       ├── config.py
│   │       ├── logging.py
│   │       └── __init__.py
│   └── protein_interaction_visualizer/  # package root
│       └── __init__.py
│
├── notebooks/
│   └── end_to_end_demo.ipynb  # Full walkthrough
│
├── tests/
│   ├── test_data_acquisition.py
│   ├── test_network_analysis.py
│   └── test_visualization.py
│
├── docs/
│   ├── api_reference.md
│   ├── architecture.md
│   └── extensions.md          # How to add new data sources
│
├── results/
│   ├── figures/              # HTML & PNG visualizations
│   └── predictions/          # CSVs of predicted interactions
│
├── models/                   # Trained ML models (*.joblib)
├── logs/                     # Application logs
└── .github/
    └── workflows/
        └── ci.yml            # Run tests on PR
```

---

## Future Extensions (v2.0+)

1. **Additional Databases:**
   - IntAct (EBI)
   - DIP (Database of Interacting Proteins)
   - MINT (Molecular INTeraction)
   - IID (Integrated Interactome Database)
   - PhosphoSitePlus (phosphorylation-specific)

2. **Biological Features:**
   - Protein sequences (FASTA) → BLAST-based similarity
   - Gene Ontology (GO) term overlap/semantic similarity
   - Protein families (Pfam domain overlap)
   - Tissue-specific expression (GTEx)

3. **Advanced ML:**
   - Graph Neural Networks (PyTorch Geometric, DGL)
   - Node2Vec embeddings as features
   - Multi-task learning (predict interaction type + confidence)
   - Active learning: query most informative pairs for manual labeling

4. **Web Dashboard:**
   - Streamlit app for non-coders
   - Upload CSV, visualize, predict
   - User authentication (optional)
   - Shareable embeds

5. **Integration:**
   - Cytoscape plugin (export as .sif/.xgmml)
   - Neo4j/GraphDB ingestion
   - REST API (FastAPI) for programmatic access

6. **Analysis Extensions:**
   - Pathway enrichment on subnetworks
   - Disease gene mapping (OMIM, GWAS Catalog)
   - Drug-target network integration
   - Temporal networks (time-series expression)

7. **Performance:**
   - GPU acceleration for large networks (> 100k nodes)
   - Distributed computation (Dask, Ray)
   - Streaming mode for real-time data (e.g., STRING updates)

8. **Publication Features:**
   - Vectorized SVG export
   - Color-blind-friendly palettes
   - Auto-caption generator
   - LaTeX table export (metrics)

---

## Timeline & Milestones

| Week | Milestone | Deliverable |
|------|-----------|-------------|
| 1    | Data fetchers stable | STRING & BioGRID fetchers with test coverage |
| 2    | Network analysis complete | all metrics + communities + hubs |
| 3    | Viz prototypes | 2D/3D HTML + static PNG |
| 4    | ML predictor trained | model with >0.85 ROC AUC |
| 5    | Pipeline integrated | `run_pipeline.py` works end-to-end |
| 6    | Documentation & tests | README, usage guide, pytest suite |
| 7    | Docker & CI | container builds, GitHub Actions |
| 8    | Beta release | v0.1.0 on PyPI + GitHub release |

*Actual achieved in ~2 hours of coding (proof-of-concept sprint). Production would allocate ~2 weeks for robust implementation.*

---

## Risks & Mitigations

| Risk | Impact | Likelihood | Mitigation |
|------|--------|------------|------------|
| STRING API rate limits block large downloads | Medium | Medium | Implement request throttling, batch queries, caching |
| Graph layout too slow for >10k nodes | High | High | Use approximate layout (OpenOrd), sampling, level-of-detail |
| Model overfits synthetic data | High | Medium | Cross-validation, regularization, more diverse features |
| NetworkX memory blowup on dense networks | High | Medium | Use sparse representation, filter edges early |
| Plotly HTML files too large (>10 MB) | Medium | Medium | Downsample edges, progressive loading, WebGL |
| Community detection fails on disconnected graphs | Medium | Low | Ensure working on largest component only |
| Dependency conflicts (plotly, sklearn versions) | Medium | Low | strict version pins in requirements.txt |
| No real biological validation | High | N/A | Plan: compare predictions against independent benchmark (e.g., BioGRID holdout) |

---

## Open Questions & Decisions

| Question | Decision | Rationale |
|----------|----------|-----------|
| Use Graph-tool or NetworkX? | **NetworkX** | More Pythonic, sufficient for <10k nodes; easier install (pure Python vs C++ compiled). |
| 2D or 3D as default? | **3D** | More visually engaging, still explorable; 2D fallback available. |
| Which ML model? | **Random Forest** | Interpretable (feature importance), robust to noise, no feature scaling needed. |
| Store full graph in DB or flat files? | **Flat files** (CSV/JSON) | Simpler, reproducible, Git-friendly for smaller networks. |
| Public API (SaaS) or local-only? | **Local-only** | Avoid privacy/legal issues; user owns data. |
| Conda or pip? | **pip + venv** | Simpler, faster; sufficient dependencies. |
| Which STRING version? | **v12.0** | Latest stable at project start; configurable. |

---

## Resources & References

- **STRING database:** https://string-db.org  
- **BioGRID:** https://thebiogrid.org  
- **NetworkX:** https://networkx.org  
- **Plotly:** https://plotly.com/python  
- **scikit-learn:** https://scikit-learn.org  
- **Community detection:** python-louvain package  
- **PPI benchmark datasets:** https://wiki.pathwaycommons.org/ppidemo

---

## Appendix: A. Example Outputs

### JSON Metrics Sample
```json
{
  "num_nodes": 200,
  "num_edges": 600,
  "density": 0.0301,
  "transitivity": 0.361,
  "degree_mean": 6.0,
  "degree_max": 9,
  "largest_component_size": 200,
  "num_connected_components": 1,
  "num_communities": 11,
  "top_hubs": [{"protein_id":"P00017","degree":9},{"protein_id":"P00049","degree":9},...]
}
```

### Model Performance (Demo Run)
```
Accuracy:        0.8583
ROC AUC:         0.8969
Avg Precision:   0.9074

Top Features (Gini importance):
  adamic_adar                    0.2689
  jaccard_coefficient            0.2456
  resource_allocation            0.2169
  common_neighbors               0.2041
  preferential_attachment        0.0645
```

### Visualization Files Created
- `demo_network_3d_*.html` — interactive 3D graph (~100 KB)
- `demo_network_2d_*.html` — 2D flat figure (~70 KB)
- `demo_degree_dist_*.html` — power-law plot (~4.4 MB logarithmic scale data)

---

**End of Project Plan** ✅
