# Quick Start Guide

Get up and running with the Protein Interaction Visualizer in 5 minutes.

---

## Prerequisites

- **Python 3.10+** (tested on 3.11)
- **Git** (to clone repository)
- (Optional) **Docker** if you prefer containerized execution

---

## Installation

### Option A: Local Installation (Recommended)

```bash
# 1. Clone the repo
git clone https://github.com/yourusername/protein-interaction-visualizer.git
cd protein-interaction-visualizer

# 2. Create and activate virtual environment
python -m venv venv
# On Windows: venv\Scripts\activate
# On Mac/Linux/WSL: source venv/bin/activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Install package in editable mode
pip install -e .

# 5. Verify installation
python -c "from src.data_acquisition import StringFetcher; print('OK')"
```

### Option B: Docker (No local Python needed)

```bash
# Build image
docker build -t ppi-visualizer .

# Run (mounts ./data to persist outputs)
docker run -p 8888:8888 -v $(pwd)/data:/app/data ppi-visualizer

# Then open browser to http://localhost:8888 for Jupyter Lab
```

---

## Running the Demo (Synthetic Data)

First, try the synthetic demo — no API keys required:

```bash
python demo.py
```

**What it does:**
1. Generates a random small-world protein network (200 nodes)
2. Computes topology metrics and detects communities
3. Generates interactive HTML visualizations (2D + 3D)
4. Trains a Random Forest predictor
5. Saves a model and sample predictions

**Outputs:**
```
results/figures/
  demo_network_2d_TIMESTAMP.html   ← interactive 2D view
  demo_network_3d_TIMESTAMP.html   ← interactive 3D view
  demo_degree_dist_TIMESTAMP.html  ← degree distribution

models/
  ppi_predictor_random_forest_TIMESTAMP.joblib  ← trained model
  scaler_TIMESTAMP.joblib
  predictor_metadata_TIMESTAMP.json

data/processed/
  demo_ppi_TIMESTAMP.csv              ← generated edges
  demo_analysis_TIMESTAMP_*.json/csv   ← metrics, hubs, communities
```

Open any `.html` file in a web browser to explore the network (drag nodes, zoom, hover for details).

---

## Fetching Real Data from STRING

For real-world analysis, fetch actual human protein interactions:

### 1. Get a STRING API Key (free)

1. Register at https://string-db.org
2. Visit https://string-db.org/cgi/access?section=api
3. Create an "app" → copy the `access_key`

### 2. Configure

Create a `.env` file in the project root:

```bash
cp .env.example .env
# Edit .env and set:
# STRING_API_KEY=your_key_here
```

Or set environment variable:

```bash
export STRING_API_KEY="your_key_here"
```

### 3. Edit `config.yaml` (optional)

```yaml
_data:
  string:
    default_organism: "9606"    # Homo sapiens (NCBI Taxonomy ID)
    score_threshold: 700        # 0-1000; higher = more confident (recommended: 700+)
    api_version: "v12.0"
    api_endpoint: "https://version-12-0.string-db.org"
```

Other organisms:
- `"4932"` — *Saccharomyces cerevisiae* (yeast)
- `"10090"` — *Mus musculus* (mouse)
- `"9606"` — *Homo sapiens* (human)
- `"3702"` — *Arabidopsis thaliana* (plant)
- Full list: https://string-db.org/cgi/input.pl?session_id=&input_page_active_form=species

### 4. Run the Pipeline

```bash
python run_pipeline.py
```

**What it does:**
- Fetches all interactions for the configured organism from STRING
- Optionally merges with BioGRID (set `include_biogrid: true` in config)
- Analyzes the network
- Generates visualizations
- Trains a predictor
- Saves everything with timestamps

**First run may take 2-5 minutes** depending on network size (~200k edges for human).

---

## Using the Jupyter Notebook

For an interactive walkthrough:

```bash
# Install Jupyter Lab (if not already)
pip install jupyterlab

# Start Jupyter
jupyter lab
```

Open `notebooks/end_to_end_demo.ipynb` and run all cells.

The notebook includes:
- Live code with explanations
- Inline visualization previews (if Kaleido installed)
- Metrics tables
- Model evaluation plots
- Exportable results

---

## Quick API Examples

### Fetch Interactions
```python
from src.data_acquisition import StringFetcher

fetcher = StringFetcher(organism="9606", score_threshold=700)
df = fetcher.fetch_interactions()
print(df.head())
#    protein1  protein2  combined_score
# 0  ENSP00000269391  ENSP00000354500              999
```

### Load Existing CSV and Visualize
```python
from src.visualization import NetworkVisualizer
import pandas as pd

df = pd.read_csv("data/processed/human_ppi_700.csv")
viz = NetworkVisualizer(df, weight_col="combined_score")
viz.plot_interactive_graph(
    title="Human PPI (score≥700)",
    output_path="results/figures/my_network.html",
    dim=3
)
```

### Analyze Network Metrics
```python
from src.network_analysis import GraphAnalyzer
import networkx as nx

G = nx.from_pandas_edgelist(df, "protein1", "protein2")
analyzer = GraphAnalyzer(G)
metrics = analyzer.compute_metrics()
communities = analyzer.detect_communities()
hubs = analyzer.find_hubs(top_n=20)

print(f"Density: {metrics['density']:.6f}")
print(f"Top hub: {hubs.iloc[0]['protein_id']} (degree={hubs.iloc[0]['degree']})")
```

### Train a Predictor
```python
from src.prediction import InteractionPredictor

predictor = InteractionPredictor(G)
X, y, pairs = predictor.prepare_features(neg_ratio=1.0, test_size=0.2)
predictor.train(predictor.train_data[0], predictor.train_data[1])
results = predictor.evaluate(*predictor.test_data)
print(f"ROC AUC = {results['roc_auc']:.3f}")

# Predict new pairs
new_pairs = [("P00533", "P04637"), ("Q9Y6K9", "O14920")]
preds = predictor.predict(new_pairs)
print(preds[["protein_A","probability","predicted"]])
```

### Load a Saved Model
```python
from src.prediction import InteractionPredictor

predictor = InteractionPredictor.load_models(
    model_path="models/ppi_predictor_latest.joblib",
    scaler_path="models/scaler_latest.joblib",
    graph=G  # optional; needed for extracting new features
)
new_preds = predictor.predict([("Q9Y6K9", "P00533")])
```

---

## Configuration Reference

Full list of options in `config.yaml`:

```yaml
# Data acquisition
data:
  string:
    default_organism: "9606"
    score_threshold: 700      # 0-1000
    api_version: "v12.0"
    batch_size: 100           # proteins per API call
  biogrid:
    evidence_types: ["physical"]
    include: false            # merge with STRING?

# Network analysis
analysis:
  compute_betweenness: true   # expensive; set false for >10k nodes
  betweenness_samples: 200    # approximate betweenness (k parameter)
  community_algorithm: "louvain"

# Visualization
visualization:
  layout:
    algorithm: "fr"           # fr, kk, circle, spectral
    iterations: 1000
  nodes:
    min_size: 10
    max_size: 60
    show_labels: false        # set true for < 100 nodes
  edges:
    width_scale: 2.0
    color: "rgba(150,150,150,0.3)"
  figures:
    width: 1200
    height: 800
    bgcolor: "#0d1117"        # dark theme

# Machine learning
prediction:
  training:
    neg_ratio: 1.0            # 1.0 = balanced (equal pos/neg)
    test_size: 0.15
    validation_size: 0.15
  model:
    algorithm: "random_forest"
    n_estimators: 200
    max_depth: 20
    min_samples_split: 10
  features:
    use_go_similarity: false   # require protein annotations
    use_domain_features: false

# Debugging
debug:
  random_seed: 42
  log_level: "INFO"           # DEBUG, INFO, WARNING, ERROR
  tqdm_disable: false
```

---

## Troubleshooting

### `ModuleNotFoundError: No module named 'src'`
Make sure you ran `pip install -e .` from the project root. Or use `python -m src...` with PYTHONPATH set.

### Cambios en el editum `pip install` says "no module named pip"
The `pip` module is missing. Run `python -m ensurepip --upgrade` first, then retry.

### Plotly HTML blank or doesn't load
Make sure you're opening the file via HTTP server, not as `file://` on some browsers due to CORS. Try:
```bash
python -m http.server 8000
# Then open http://localhost:8000/results/figures/demo_network_3d_*.html
```
But stand-alone Plotly HTML usually works as file://. If not, it's a WebGL/driver issue.

### Network layout takes forever for >5000 nodes
Adjust config: `visualization.layout.iterations: 300` and use `layout_algorithm: "kk"` or `"circle"`.

### Score threshold too high → too few edges
Lower the threshold in config or CLI: `--score 400`. STRING scores: 0-1000; <150 low, 400-700 medium, >700 high.

### `NetworkXError: v is not in the graph.`
You passed a protein ID that wasn't in the fetched network. Use only IDs from your data (check `df["protein1"].unique()`).

### Model predicts all zeros or all ones
Check class balance: if data is highly imbalanced, set `class_weight: "balanced"` in config.

### Docker build fails on `apt-get`
Some base images require `apt-get update` first or you're offline. Use `--network=host` or ensure internet access.

### Tests fail because of mocked API
Make sure you have `responses` library installed: `pip install responses`.

---

## Support & Contributing

- **Issues:** https://github.com/yourusername/protein-interaction-visualizer/issues
- **Discussions:** https://github.com/yourusername/protein-interaction-visualizer/discussions
- **Email:** your.email@example.com

Contributions welcome! See `CONTRIBUTING.md` for guidelines.

---

## Next Steps

1. **Real Data Run:** Follow "Fetching Real Data" section above to analyze human PPI.
2. **Customize:** Modify config or code for your organism/tissue of interest.
3. **Extend:** Add a new data source (see docs/extensions.md).
4. **Publish:** Use generated figures in your paper (cite this repo & STRING).
5. **Share:** Tweet your findings with #ProteinViz!

---

**Happy network exploring!** 🧬🔬
