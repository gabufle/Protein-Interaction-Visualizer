#!/usr/bin/env python3
"""
Demo: Proof-of-Concept Pipeline with Synthetic Data

This demonstrates all pipeline steps without requiring API keys:
  1. Create synthetic PPI network (small-world graph)
  2. Analyze network metrics
  3. Generate interactive 3D visualization
  4. Train ML predictor on synthetic features
  5. Save all outputs

Run:
    python demo.py

Outputs:
    data/processed/demo_network.csv
    results/figures/demo_3d.html
    models/demo_predictor.joblib
"""

import sys
from pathlib import Path
import json

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

import networkx as nx
import pandas as pd
import numpy as np
from datetime import datetime

# Import project modules
from src.network_analysis import GraphAnalyzer
from src.visualization import NetworkVisualizer
from src.prediction import InteractionPredictor
import random
from src.utils.logging import get_logger, silence_third_party_loggers

# Silence external loggers for cleaner output
silence_third_party_loggers()
logger = get_logger(__name__)

# Project paths
PROJECT_ROOT = Path(__file__).parent
DATA_PROC = PROJECT_ROOT / "data" / "processed"
RESULTS_FIG = PROJECT_ROOT / "results" / "figures"
RESULTS_PRED = PROJECT_ROOT / "results" / "predictions"
MODELS_DIR = PROJECT_ROOT / "models"
LOGS_DIR = PROJECT_ROOT / "logs"

for p in [DATA_PROC, RESULTS_FIG, RESULTS_PRED, MODELS_DIR, LOGS_DIR]:
    p.mkdir(parents=True, exist_ok=True)

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")


# ============================================
# STEP 1: Generate Synthetic PPI Dataset
# ============================================
print("\n" + "="*60)
print("STEP 1: Generating Synthetic Human PPI Dataset")
print("="*60)

np.random.seed(42)

# Create a small-world network (more realistic than random)
n_nodes = 200
avg_degree = 6

G_synthetic = nx.watts_strogatz_graph(
    n=n_nodes,
    k=avg_degree,
    p=0.15,
    seed=42,
)

# Convert to edge list DataFrame
edges = []
for u, v in G_synthetic.edges():
    # Simulate STRING confidence scores (0-1000) — biased toward higher for known edges
    score = np.random.randint(600, 999)
    edges.append({
        "protein1": f"P{u:05d}",
        "protein2": f"P{v:05d}",
        "combined_score": score,
    })

df_interactions = pd.DataFrame(edges)
print(f"✓ Generated {len(df_interactions):,} edges")
print(f"✓ Unique proteins: {len(set(df_interactions['protein1']) | set(df_interactions['protein2']))}")

# Save
csv_path = DATA_PROC / f"demo_ppi_{timestamp}.csv"
df_interactions.to_csv(csv_path, index=False)
print(f"✓ Saved to {csv_path}")


# ============================================
# STEP 2: Network Analysis
# ============================================
print("\n" + "="*60)
print("STEP 2: Network Topology Analysis")
print("="*60)

# Build NetworkX graph
G = nx.from_pandas_edgelist(
    df_interactions,
    source="protein1",
    target="protein2",
    edge_attr="combined_score",
    create_using=nx.Graph(),
)

analyzer = GraphAnalyzer(G)
metrics = analyzer.compute_metrics()
communities = analyzer.detect_communities(algorithm="louvain")
hubs = analyzer.find_hubs(top_n=10)

print(f"\nNetwork Properties:")
print(f"  Nodes:         {metrics['num_nodes']:,}")
print(f"  Edges:         {metrics['num_edges']:,}")
print(f"  Density:       {metrics['density']:.6f}")
print(f"  Transitivity:  {metrics['transitivity']:.4f}")
print(f"  Avg Degree:    {metrics['degree_mean']:.2f}")
print(f"  Max Degree:    {metrics['degree_max']}")
print(f"  Communities:   {len(communities)}")

print(f"\nTop 5 Hub Proteins (by degree):")
print(hubs[["protein_id", "degree", "centrality_score"]].head(5).to_string(index=False))

# Save analysis
analyzer.save_results(
    metrics, hubs, communities,
    output_dir=DATA_PROC,
    prefix=f"demo_analysis_{timestamp}",
)
print(f"\n✓ Analysis results saved to {DATA_PROC}")


# ============================================
# STEP 3: Interactive Visualization
# ============================================
print("\n" + "="*60)
print("STEP 3: Interactive 3D Network Visualization")
print("="*60)

viz = NetworkVisualizer(df_interactions)

# 2D view
fig2d = viz.plot_interactive_graph(
    title="Synthetic Human PPI Network (2D) — Demo",
    output_path=RESULTS_FIG / f"demo_network_2d_{timestamp}.html",
    dim=2,
    show_labels=False,
    layout_iterations=300,
    node_color_col="community",
)
print(f"✓ 2D visualization → {RESULTS_FIG / f'demo_network_2d_{timestamp}.html'}")

# 3D view
fig3d = viz.plot_interactive_graph(
    title="Synthetic Human PPI Network (3D) — Demo",
    output_path=RESULTS_FIG / f"demo_network_3d_{timestamp}.html",
    dim=3,
    show_labels=False,
    layout_iterations=300,
    node_color_col="community",
)
print(f"✓ 3D visualization → {RESULTS_FIG / f'demo_network_3d_{timestamp}.html'}")

# Degree distribution
deg_fig = viz.plot_degree_distribution(
    output_path=RESULTS_FIG / f"demo_degree_dist_{timestamp}.html",
)
print(f"✓ Degree distribution plot → {RESULTS_FIG / f'demo_degree_dist_{timestamp}.html'}")


# ============================================
# STEP 4: Machine Learning Prediction
# ============================================
print("\n" + "="*60)
print("STEP 4: Interaction Prediction (Machine Learning)")
print("="*60)

predictor = InteractionPredictor(G)

# Prepare features
X, y, all_pairs = predictor.prepare_features(
    neg_ratio=1.0,  # balanced
    test_size=0.2,
    val_size=0.1,
)

X_train, y_train = predictor.train_data
X_test, y_test = predictor.test_data

print(f"Feature matrix: {X.shape[0]:,} samples × {X.shape[1]} features")
print(f"Training set: {len(X_train):,}")
print(f"Test set: {len(X_test):,}")
print(f"Positive ratio: {y.mean():.2%}")

# Train model
print("\nTraining Random Forest...")
predictor.train(X_train, y_train, model_type="random_forest", n_estimators=100)

# Evaluate
results = predictor.evaluate(X_test, y_test)
print(f"\nModel Performance:")
print(f"  Accuracy:        {results['accuracy']:.4f}")
print(f"  ROC AUC:         {results['roc_auc']:.4f}")
print(f"  Avg Precision:   {results['average_precision']:.4f}")

# Feature importance
if predictor.feature_importance is not None:
    print(f"\nTop Features:")
    for _, row in predictor.feature_importance.head(5).iterrows():
        print(f"  {row['feature']:<30} {row['importance']:.4f}")

# Save model
predictor.save_models(MODELS_DIR)
print(f"\n✓ Model saved to {MODELS_DIR}")

# Example prediction on new pairs
print("\nPredicting on novel pairs:")
# Use real proteins from the graph
import random
random.seed(123)
all_nodes = list(G.nodes())
test_nodes = random.sample(all_nodes, min(6, len(all_nodes)))
test_new = [
    (test_nodes[0], test_nodes[1]),
    (test_nodes[2], test_nodes[3]),
    (test_nodes[4], test_nodes[5]),
]
predictions = predictor.predict(test_new)
print(predictions.to_string(index=False))

pred_csv = RESULTS_PRED / f"demo_predictions_{timestamp}.csv"
predictor.save_predictions(predictions, pred_csv)
print(f"\n✓ Predictions saved to {pred_csv}")


# ============================================
# STEP 5: Summary Report
# ============================================
print("\n" + "="*60)
print("STEP 5: Pipeline Summary Report")
print("="*60)

summary = {
    "pipeline": "Protein Interaction Visualizer — Demo Run",
    "timestamp": datetime.now().isoformat(),
    "dataset": {
        "type": "synthetic",
        "generator": "Watts-Strogatz small-world (n=200, k=6, p=0.15)",
        "num_interactions": len(df_interactions),
        "confidence_range": [int(df_interactions['combined_score'].min()),
                            int(df_interactions['combined_score'].max())],
    },
    "network": {
        "nodes": metrics['num_nodes'],
        "edges": metrics['num_edges'],
        "density": metrics['density'],
        "clustering": metrics['transitivity'],
        "avg_degree": metrics['degree_mean'],
        "communities": len(communities),
    },
    "model": {
        "type": predictor.model_type,
        "features": predictor.feature_names,
        "accuracy": results['accuracy'],
        "roc_auc": results['roc_auc'],
        "avg_precision": results['average_precision'],
    },
    "outputs": {
        "network_csv": str(csv_path.relative_to(PROJECT_ROOT)),
        "analysis_dir": str(DATA_PROC.relative_to(PROJECT_ROOT)),
        "figures": [str(p.relative_to(PROJECT_ROOT)) for p in RESULTS_FIG.glob(f"*{timestamp}*")],
        "model_dir": str(MODELS_DIR.relative_to(PROJECT_ROOT)),
        "predictions": str(pred_csv.relative_to(PROJECT_ROOT)),
    },
}

report_path = LOGS_DIR / f"demo_report_{timestamp}.json"
with open(report_path, "w") as f:
    json.dump(summary, f, indent=2)

print(f"\n✓ Full JSON report → {report_path}")

print("\n" + "="*60)
print("✅ DEMO PIPELINE COMPLETE")
print("="*60)
print("\nAll outputs saved to project directories:")
print(f"  📁 data/processed/      → network + analysis CSVs")
print(f"  📁 results/figures/     → HTML interactive visualizations")
print(f"  📁 results/predictions/ → ML predictions")
print(f"  📁 models/              → trained ML models")
print("\nNext: Open the HTML figures in a browser to explore!")
print(f"      Open: {RESULTS_FIG / f'demo_network_3d_{timestamp}.html'}")
print("\n" + "="*60)
