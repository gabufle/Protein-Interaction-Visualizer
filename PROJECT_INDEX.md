# Project Index & Reference Card

**Protein Interaction Visualizer & Predictor** — v0.1.0 (POC)

> A reproducible Python pipeline for fetching, analyzing, visualizing, and predicting protein-protein interactions.

---

## Quick Links

| What you need | Where to look |
|----------------|---------------|
| Getting started | `QUICKSTART.md` |
| Full project plan | `PROJECT_PLAN.md` |
| API reference | `docs/api_reference.md` *(if created)* |
| Demo notebook | `notebooks/end_to_end_demo.ipynb` |
| Run pipeline | `python run_pipeline.py` |
| Run demo (synthetic) | `python demo.py` |
| Edit config | `config.yaml` |
| Add dependencies | `requirements.txt` |
| Docker | `Dockerfile` |
| Tests | `tests/` |
| Make shortcuts | `Makefile` → `make demo`, `make test`, etc. |

---

## File-By-File Guide

### Top-Level (Project Root)
| File | Purpose |
|------|---------|
| `README.md` | Project overview, key features, badges |
| `PROJECT_PLAN.md` | Comprehensive architecture, phases, decisions |
| `QUICKSTART.md` | 5-minute guide to running the code |
| `CONTRIBUTING.md` | Guidelines for collaborators |
| `CHANGELOG.md` | Version history |
| `LICENSE` | MIT License |
| `requirements.txt` | Pinned Python dependencies |
| `pyproject.toml` | Build system config (PEP 621) |
| `setup.py` | Editable install wrapper |
| `config.yaml` | Tunable parameters (organism, thresholds, viz options) |
| `.env.example` | Template for STRING_API_KEY |
| `Makefile` | Shortcut commands (`make demo`, `make test`, `make docker-build`) |
| `.pre-commit-config.yaml` | Code quality hooks (black, isort, flake8, mypy) |
| `.gitignore` | Excludes: `__pycache__`, `*.pyc`, `data/raw/*`, logs, models |

### Source Code (`src/`)
Organized as a Python package `protein_interaction_visualizer` (install via `pip install -e .`).

| Module | Key Files | Classes / Functions | Description |
|--------|-----------|---------------------|-------------|
| **data_acquisition** | `string_fetcher.py`<br>`biogrid_fetcher.py` | `StringFetcher`<br>`BioGRIDFetcher`<br>`merge_string_biogrid()` | Fetches PPI from STRING API and BioGRID REST API |
| **network_analysis** | `graph_analyzer.py` | `GraphAnalyzer`<br>`load_network_from_csv()` | NetworkX-based metrics, communities, hubs |
| **visualization** | `network_viz.py` | `NetworkVisualizer`<br>`visualize_from_csv()` | Interactive Plotly 3D/2D HTML + static PNG |
| **prediction** | `interaction_predictor.py` | `InteractionPredictor`<br>`load_and_predict()` | ML pipeline: feature engineering, training, eval |
| **utils** | `config.py`<br>`logging.py` | `load_config()`<br>`get_logger()`<br>`silence_third_party_loggers()` | YAML config loader, colored logging, reproducibility helpers |

### Scripts
| Script | Purpose | Usage |
|--------|---------|-------|
| `demo.py` | Synthetic data demo (no API key) | `python demo.py` |
| `run_pipeline.py` | Full pipeline with real data | `python run_pipeline.py` |

### Notebooks
| Notebook | Purpose |
|----------|---------|
| `notebooks/end_to_end_demo.ipynb` | Step-by-step walkthrough with explanations, plots inline |

### Tests (`tests/`)
| Test file | Coverage |
|-----------|----------|
| `test_data_acquisition.py` | STRING/BioGRID fetchers, merging, deduplication |
| `test_network_analysis.py` | GraphAnalyzer metrics, community detection, hub ranking |
| `test_visualization.py` | NetworkVisualizer, layout, HTML generation |

Run: `pytest tests/ -v --cov=src`

### Outputs (Generated at Runtime)

| Directory | Files | Meaning |
|-----------|-------|---------|
| `data/processed/` | `*.csv` network edge lists<br>`*_metrics.json`<br>`*_hubs.csv`<br>`*_communities.json` | Processed data + analysis results |
| `results/figures/` | `*_2d.html`, `*_3d.html`<br>`*_degree_dist.html` | Interactive visualizations (open in browser) |
| `results/predictions/` | `*.csv` | Predictions of novel protein pairs |
| `models/` | `ppi_predictor_*.joblib`<br>`scaler_*.joblib`<br>`predictor_metadata_*.json` | Trained ML models and metadata |
| `logs/` | `demo_report_*.json`<br>`protein_vis.log` | JSON pipeline summary, detailed logs |

---

## Common Tasks (Cheatsheet)

### First time setup
```bash
# Clone
git clone <url>
cd protein-interaction-visualizer

# Create venv (optional but recommended)
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install
pip install -r requirements.txt
pip install -e .

# Test
python demo.py
```

### Fetch real human PPI network (STRING)
```bash
# 1. Get API key at https://string-db.org/cgi/access?section=api
# 2. Set it:
export STRING_API_KEY="your_key"
# OR edit config.yaml

# 3. Run pipeline
python run_pipeline.py
```

### Jupyter Lab
```bash
pip install jupyterlab
jupyter lab
# Open notebooks/end_to_end_demo.ipynb
```

### Run tests
```bash
pytest tests/ -v
# with coverage:
pytest tests/ -v --cov=src --cov-report=html
# then open htmlcov/index.html
```

### Pre-commit checks
```bash
pre-commit install
# Now automatic on git commit. To run manually:
pre-commit run --all-files
```

### Format code
```bash
make format
# or directly:
black src/ tests/
isort src/ tests/
```

### Docker
```bash
docker build -t ppi-vis .
docker run -p 8888:8888 -v $(pwd)/data:/app/data ppi-vis
# Open http://localhost:8888 for Jupyter Lab
```

### Makefile shortcuts
```bash
make help          # show all targets
make demo          # run synthetic demo
make test          # run pytest
make lint          # check style
make docker-build  # build container
make clean         # remove cache files
```

### Use as library in your script
```python
import sys
sys.path.insert(0, "src")

from src.data_acquisition import StringFetcher
from src.visualization import NetworkVisualizer
from src.prediction import InteractionPredictor

# fetch
fetcher = StringFetcher(organism="9606", score_threshold=700)
df = fetcher.fetch_interactions()

# viz
viz = NetworkVisualizer(df)
viz.plot_interactive_graph(output_path="my_ppi.html")

# predict
G = nx.from_pandas_edgelist(df, "protein1", "protein2")
predictor = InteractionPredictor(G)
predictor.train(...)
```

---

## Output Interpretation

### Network Metrics JSON
```json
{
  "num_nodes": 200,
  "num_edges": 600,
  "density": 0.030,
  "transitivity": 0.36,
  "degree_mean": 6.0,
  "degree_max": 9,
  "largest_component_size": 200,
  "num_connected_components": 1,
  "diameter": 8,
  "avg_shortest_path": 3.2,
  "degree_centrality": {...}
}
```

### Model Metrics JSON
```json
{
  "accuracy": 0.858,
  "roc_auc": 0.897,
  "average_precision": 0.907,
  "classification_report": {
    "0": {"precision": ..., "recall": ..., "f1-score": ...},
    "1": {...}
  },
  "confusion_matrix": [[TN, FP], [FN, TP]],
  "feature_importance": [
    {"feature": "adamic_adar", "importance": 0.269},
    ...
  ]
}
```

### HTML Interactive Figures
- Open in any modern browser
- Drag nodes to rearrange
- Scroll to zoom; right-click drag to rotate (3D)
- Hover to see protein ID, degree, community
- Click camera icon to save screenshot as PNG

---

## Reproducing Results

### To re-run exact demo (same synthetic data):
```bash
python demo.py  # uses fixed random seed (42) in script
# Outputs will have same timestamp pattern; check logs/demo_report_*.json for config
```

### To replicate on your machine:
```bash
git clone <repo>
cd protein-interaction-visualizer
conda create -n ppi python=3.11 -y
conda activate ppi
pip install -r requirements.txt
pip install -e .
python demo.py
# All timestamps will differ but graphs should be isomorphic (same seed)
```

### To reproduce with real data (STRING human):
1. Get STRING key
2. `config.yaml` → `organism: "9606"`, `score_threshold: 700`
3. `python run_pipeline.py`
4. Results in `data/processed/`, `results/figures/`, `models/`

---

## FAQ

**Q: Can I analyze a network from my own CSV file?**  
A: Yes! Use `viz = NetworkVisualizer(pd.read_csv("my_edges.csv"))`.

**Q: How to change number of communities detected?**  
A: In `config.yaml` set `community_algorithm: "label_propagation"` for fewer, coarser communities, or `"girvan_newman"` for more granular.

**Q: My graph is huge (100k nodes). Will it work?**  
A: Layout will be slow. Set `layout.iterations: 100`, `layout.algorithm: "kk"`, `compute_betweenness: false`. Consider subsetting: `df = df.nlargest(5000, "combined_score")`.

**Q: I want to predict with my own proteins not in network.**  
A: `InteractionPredictor.predict()` requires proteins to be in graph `G` to compute topological features. If novel, can't predict (by definition features undefined). Either add them to graph with edges to some neighbors or use biological features if available.

**Q: Can I export to Cytoscape?**  
A: Yes, manually: `nx.write_edgelist(G, "network.sif", delimiter="\t")` format is SIF; Cytoscape imports that. Future: `export_to_cytoscape()` helper.

**Q: How do I cite this?**  
A: If you use this software in research, please cite:  
   **Your Name** (2025). *Protein Interaction Visualizer* (v0.1.0). https://github.com/yourusername/protein-interaction-visualizer  
   Also cite the underlying databases: STRING (Szklarczyk et al., NAR 2021), BioGRID (Oughtred et al., NAR 2021).

**Q: I got an error about missing 'community' module.**  
A: Install `python-louvain`: `pip install python-louvain`.

**Q: Plotly HTML file shows blank white page.**  
A: Ensure you're opening via HTTP server for some browsers; try `python -m http.server` and navigate. Or try a different browser.

**Q: Model training is too slow for my big graph.**  
A: Reduce `n_estimators` in config, or use `model: "logistic_regression"` for linear baseline. Also, sample edges: `df = df.sample(50000)`.

---

## Glossary

- **PPI** — Protein-Protein Interaction
- **STRING** — Search Tool for the Retrieval of Interacting Genes/Proteins (database)
- **BioGRID** — Biological General Repository for Interaction Datasets
- **NetworkX** — Python graph library
- **Louvain** — Community detection algorithm maximizing modularity
- **Force-directed layout** — Graph drawing where edges are springs
- **ROC AUC** — Area under receiver operating characteristic curve
- **Average Precision** — Area under precision-recall curve
- **Joblib** — Python persistence format (`.joblib` files)
- **Kaleido** — Static image export for Plotly

---

## Support

- **Issues:** https://github.com/yourusername/protein-interaction-visualizer/issues
- **Email:** your.email@example.com
- **Discord/Slack:** (if you have community)

---

*Generated 2026-05-14 | v0.1.0 proof of concept*
