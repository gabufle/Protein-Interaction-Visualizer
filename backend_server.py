#!/usr/bin/env python3
"""
FastAPI backend server for Protein Interaction Visualizer & Predictor.

Serves the React frontend and provides REST API endpoints for:
- Network data (nodes, edges, layout positions)
- Network analysis metrics
- Community detection results
- Hub protein rankings
- ML model performance & predictions
- Degree distribution data
- Feature importance data

Run:
    python backend_server.py
    # Then open http://localhost:8000
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Any

import networkx as nx
import numpy as np
import pandas as pd
from fastapi import FastAPI, HTTPException, Query, UploadFile, File
import io
from itertools import combinations
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
import shutil
import subprocess
import tempfile
import os
import re
import requests

# Add src to path
PROJECT_ROOT = Path(__file__).parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from src.network_analysis.graph_analyzer import GraphAnalyzer, load_network_from_csv
from src.visualization.network_viz import NetworkVisualizer
from src.prediction.interaction_predictor import InteractionPredictor
from src.utils.config import load_config
from src.utils.logging import get_logger, silence_third_party_loggers

silence_third_party_loggers()
logger = get_logger(__name__)

app = FastAPI(
    title="PPI Visualizer API",
    description="Protein-Protein Interaction Network Analysis & Prediction",
    version="1.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ─── Global State ───────────────────────────────────────────────────────────
config = load_config()
G: Optional[nx.Graph] = None
df_edges: Optional[pd.DataFrame] = None
analyzer: Optional[GraphAnalyzer] = None
visualizer: Optional[NetworkVisualizer] = None
predictor: Optional[InteractionPredictor] = None
analysis_results: Dict[str, Any] = {}
model_results: Dict[str, Any] = {}


# ─── Data Loading ───────────────────────────────────────────────────────────

def generate_synthetic_data():
    """Generate synthetic PPI network for demo."""
    global G, df_edges, analyzer, visualizer, predictor, analysis_results, model_results

    np.random.seed(42)
    n_nodes = 200
    avg_degree = 6

    G_synthetic = nx.watts_strogatz_graph(n=n_nodes, k=avg_degree, p=0.15, seed=42)

    edges = []
    for u, v in G_synthetic.edges():
        score = int(np.random.randint(600, 999))
        edges.append({
            "protein1": f"P{u:05d}",
            "protein2": f"P{v:05d}",
            "combined_score": score,
        })

    df_edges = pd.DataFrame(edges)

    G = nx.from_pandas_edgelist(
        df_edges,
        source="protein1",
        target="protein2",
        edge_attr="combined_score",
        create_using=nx.Graph(),
    )

    logger.info(f"Loaded synthetic network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G, df_edges


def run_analysis():
    """Run full network analysis."""
    global analyzer, visualizer, analysis_results

    analyzer = GraphAnalyzer(G)
    metrics = analyzer.compute_metrics()
    communities = analyzer.detect_communities(algorithm="louvain")
    hubs = analyzer.find_hubs(top_n=20)

    # Build node-to-community mapping
    node_to_comm = {}
    for comm_id, members in communities.items():
        for node in members:
            node_to_comm[node] = comm_id

    # Degree distribution
    degrees = [d for _, d in G.degree()]
    degree_hist, bin_edges = np.histogram(degrees, bins=20)
    degree_dist = {
        "bins": bin_edges.tolist(),
        "counts": degree_hist.tolist(),
        "degrees": degrees,
    }

    # Community sizes
    comm_sizes = {str(k): len(v) for k, v in communities.items()}

    # Centrality data for all nodes
    deg_cent = nx.degree_centrality(G)
    betw_cent = nx.betweenness_centrality(G, k=min(200, len(G)), seed=42)

    analysis_results = {
        "metrics": {
            "num_nodes": metrics["num_nodes"],
            "num_edges": metrics["num_edges"],
            "density": round(metrics["density"], 6),
            "transitivity": round(metrics["transitivity"], 4),
            "avg_clustering_coefficient": round(metrics["avg_clustering_coefficient"], 4),
            "degree_mean": round(metrics["degree_mean"], 2),
            "degree_median": round(metrics["degree_median"], 2),
            "degree_max": metrics["degree_max"],
            "degree_std": round(metrics["degree_std"], 2),
            "num_connected_components": metrics.get("num_connected_components", 1),
            "largest_component_size": metrics["largest_component_size"],
            "largest_component_ratio": round(metrics["largest_component_ratio"], 4),
            "diameter": metrics.get("diameter"),
            "avg_shortest_path_length": round(metrics["avg_shortest_path_length"], 4) if metrics.get("avg_shortest_path_length") else None,
            "degree_assortativity": round(metrics["degree_assortativity"], 4) if metrics.get("degree_assortativity") is not None else None,
            "degree_centrality_mean": round(metrics.get("degree_centrality_mean", 0), 4),
            "degree_centrality_max": round(metrics.get("degree_centrality_max", 0), 4),
            "closeness_centrality_mean": round(metrics.get("closeness_centrality_mean", 0), 4),
            "num_communities": len(communities),
        },
        "communities": comm_sizes,
        "hubs": hubs.to_dict(orient="records"),
        "degree_distribution": degree_dist,
        "node_attributes": {
            node: {
                "degree": G.degree(node),
                "community": node_to_comm.get(node, -1),
                "degree_centrality": round(deg_cent.get(node, 0), 6),
                "betweenness_centrality": round(betw_cent.get(node, 0), 6),
            }
            for node in G.nodes()
        },
    }

    # Create visualizer
    visualizer = NetworkVisualizer(df_edges, community_labels=node_to_comm)

    return analysis_results


def run_prediction():
    """Run ML prediction pipeline."""
    global predictor, model_results

    predictor = InteractionPredictor(G)
    X, y, all_pairs = predictor.prepare_features(neg_ratio=1.0, test_size=0.2, val_size=0.1)
    X_train, y_train = predictor.train_data
    X_test, y_test = predictor.test_data

    predictor.train(X_train, y_train, model_type="random_forest", n_estimators=100)
    results = predictor.evaluate(X_test, y_test)

    # Feature importance
    feature_importance = []
    if predictor.feature_importance is not None:
        feature_importance = predictor.feature_importance.to_dict(orient="records")

    # ROC curve points
    from sklearn.metrics import roc_curve, precision_recall_curve
    X_test_scaled = predictor.scaler.transform(X_test)
    y_proba = predictor.model.predict_proba(X_test_scaled)[:, 1]

    fpr, tpr, _ = roc_curve(y_test, y_proba)
    precision, recall, _ = precision_recall_curve(y_test, y_proba)

    # Confusion matrix
    y_pred = predictor.model.predict(X_test_scaled)
    from sklearn.metrics import confusion_matrix
    cm = confusion_matrix(y_test, y_pred)

    # Classification report
    from sklearn.metrics import classification_report
    report = classification_report(y_test, y_pred, output_dict=True)

    model_results = {
        "accuracy": round(results["accuracy"], 4),
        "roc_auc": round(results["roc_auc"], 4),
        "average_precision": round(results["average_precision"], 4),
        "feature_importance": feature_importance,
        "roc_curve": {
            "fpr": fpr.tolist(),
            "tpr": tpr.tolist(),
        },
        "pr_curve": {
            "precision": precision.tolist(),
            "recall": recall.tolist(),
        },
        "confusion_matrix": cm.tolist(),
        "classification_report": report,
        "training_samples": len(X_train),
        "test_samples": len(X_test),
        "n_features": X.shape[1],
        "feature_names": predictor.feature_names,
        "model_type": predictor.model_type,
    }

    return model_results


def compute_network_layout(dim: int = 3) -> Dict[str, List[float]]:
    """Compute node positions for the network visualization."""
    pos = visualizer._compute_layout(algorithm="fr", dim=dim, iterations=300, seed=42)
    return {node: list(coords) for node, coords in pos.items()}


# ── Static test data ──────────────────────────────────────────────────────────
TEST_FASTA_DIR = PROJECT_ROOT / "test_fastas"
if TEST_FASTA_DIR.exists():
    from fastapi.staticfiles import StaticFiles
    app.mount("/test_fastas", StaticFiles(directory=str(TEST_FASTA_DIR), html=False), name="test_fastas")

# ── FASTA ingestion ──────────────────────────────────────────────────────────

def _parse_fasta(fasta_text: str) -> dict[str, str]:
    """
    Return {header_without_gt: sequence} supporting multi-line records.
    """
    seqs: dict[str, str] = {}
    current_header: str | None = None
    parts: list[str] = []
    for raw_line in fasta_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                seqs[current_header] = "".join(parts).upper()
            current_header = line[1:].split()[0]  # first token, no spaces
            parts = []
        else:
            parts.append(line)
    if current_header is not None:
        seqs[current_header] = "".join(parts).upper()
    return seqs


def _build_network_from_fasta(
    seqs: dict[str, str], label: str = "FA"
) -> tuple[nx.Graph, pd.DataFrame]:
    """
    Deterministic builder: the graph structure is a function of each sequence's
    length so different proteins produce different node/edge counts and density
    values.
    """
    all_nodes: list[str]   = []
    all_seqs:  dict[str, str] = {}
    edges:     list[dict]  = []

    rng = np.random.default_rng(42)  # fixed seed for reproducibility

    for idx, (header, seq) in enumerate(seqs.items()):
        L = len(seq)
        n_nodes = np.clip(L // 10, 15, 50)
        prefix  = f"{label}_{idx}"
        seqs_local: dict[str, str] = {}

        for i in range(int(n_nodes)):
            new_len = max(10, L + rng.integers(-10, 11))
            new_seq = "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), size=new_len))
            nid = f"{prefix}_{i}"
            seqs_local[nid] = new_seq

        max_diff = max(3, L // 100)

        for a, b in combinations(seqs_local.keys(), 2):
            diff = abs(len(seqs_local[a]) - len(seqs_local[b]))
            if diff <= max_diff:
                score = int(100 * (1 - diff / max_diff))
                edges.append({"protein1": a, "protein2": b, "combined_score": max(score, 1)})

        all_nodes.extend(seqs_local.keys())
        all_seqs.update(seqs_local)

    G = nx.Graph()
    G.add_nodes_from(all_nodes)
    G.add_edges_from((e["protein1"], e["protein2"]) for e in edges)
    df = pd.DataFrame(edges)
    logger.info(f"Built network from {len(seqs)} FASTA entries: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G, df


def load_fasta_into_backend(
    fasta_text: str, label: str = "FA"
) -> dict[str, Any]:
    """Replace the global graph and re-run the full pipeline."""
    global G, df_edges, analyzer, visualizer, predictor, analysis_results, model_results

    seqs = _parse_fasta(fasta_text)
    if not seqs:
        raise ValueError("No sequences found in FASTA")

    G_new, df_new = _build_network_from_fasta(seqs, label)
    G        = G_new
    df_edges = df_new

    run_analysis()
    run_prediction()
    return {
        "status":    "ok",
        "sequences": len(seqs),
        "nodes":     G.number_of_nodes(),
        "edges":     G.number_of_edges(),
    }


# ── Request models ──────────────────────────────────────────────────────────

class FASTATextRequest(BaseModel):
    fasta_text: str
    label:      str = "FA"


class PDBTextRequest(BaseModel):
    pdb_text: str
    filename: str = "uploaded.pdb"


# ── New ingestion endpoints ─────────────────────────────────────────────────

@app.post("/api/load-fasta", response_model=dict)
async def api_load_fasta(req: FASTATextRequest):
    try:
        result = load_fasta_into_backend(req.fasta_text, req.label)
        return result
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))


@app.post("/api/load-fasta-file", response_model=dict)
async def api_load_fasta_file(
    label:      str       = Query("FA", description="Node-id prefix"),
    fasta_file: UploadFile = File(...),
):
    try:
        content = (await fasta_file.read()).decode("utf-8", errors="replace")
        result  = load_fasta_into_backend(content, label)
        return result
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))


@app.post("/api/upload-pdb", response_model=dict)
async def api_upload_pdb(req: PDBTextRequest):
    """Stub – accepts PDB text and echoes back acceptance."""
    _pdb_text = req.pdb_text
    _pdb_name = req.filename
    return {"status": "ok", "filename": _pdb_name}


# ─── Initialize on startup ─────────────────────────────────────────────────

@app.on_event("startup")
async def startup_event():
    logger.info("Initializing PPI Visualizer backend...")
    generate_synthetic_data()
    run_analysis()
    run_prediction()
    logger.info("Backend ready!")


# ─── API Endpoints ──────────────────────────────────────────────────────────

@app.get("/api/health")
async def health_check():
    return {"status": "ok", "nodes": G.number_of_nodes() if G else 0, "edges": G.number_of_edges() if G else 0}


@app.get("/api/network/summary")
async def network_summary():
    """Get network summary statistics."""
    return analysis_results.get("metrics", {})


@app.get("/api/network/nodes")
async def network_nodes():
    """Get all nodes with attributes."""
    return analysis_results.get("node_attributes", {})


@app.get("/api/network/edges")
async def network_edges(limit: int = Query(1000, ge=1, le=10000)):
    """Get edge list."""
    if df_edges is None:
        return []
    return df_edges.head(limit).to_dict(orient="records")


@app.get("/api/network/layout")
async def network_layout(dim: int = Query(3, ge=2, le=3)):
    """Get node positions for visualization."""
    layout = compute_network_layout(dim=dim)
    return layout


@app.get("/api/network/communities")
async def network_communities():
    """Get community assignments."""
    return analysis_results.get("communities", {})


@app.get("/api/network/hubs")
async def network_hubs(top_n: int = Query(20, ge=1, le=100)):
    """Get top hub proteins."""
    hubs = analysis_results.get("hubs", [])
    return hubs[:top_n]


@app.get("/api/network/degree-distribution")
async def degree_distribution():
    """Get degree distribution data."""
    return analysis_results.get("degree_distribution", {})


@app.get("/api/analysis/metrics")
async def analysis_metrics():
    """Get full network analysis metrics."""
    return analysis_results.get("metrics", {})


@app.get("/api/prediction/results")
async def prediction_results():
    """Get ML model results."""
    return model_results


@app.get("/api/prediction/feature-importance")
async def feature_importance():
    """Get feature importance from trained model."""
    return model_results.get("feature_importance", [])


@app.get("/api/prediction/roc-curve")
async def roc_curve():
    """Get ROC curve data."""
    return model_results.get("roc_curve", {})


@app.get("/api/prediction/pr-curve")
async def pr_curve():
    """Get Precision-Recall curve data."""
    return model_results.get("pr_curve", {})


@app.get("/api/prediction/confusion-matrix")
async def confusion_matrix():
    """Get confusion matrix."""
    return {"matrix": model_results.get("confusion_matrix", [])}


class PredictRequest(BaseModel):
    protein_a: str
    protein_b: str


class DockingRequest(BaseModel):
    protein_a: str
    protein_b: str


# ─── Docking Globals ───────────────────────────────────────────────────────────

VINA_BINARY      = shutil.which("vina")
PREP_RECEPTOR    = shutil.which("prepare_receptor4.py")
PREP_LIGAND      = shutil.which("prepare_ligand4.py")

# Regex: exactly 4 uppercase alphanumeric characters (classic PDB identifiers)
PDB_ID_RE       = re.compile(r"^[A-Z0-9]{4}$")

RCSB_PDB_URL    = "https://files.rcsb.org/download/{pdb_id}.pdb"

VINA_BOX_CENTER = {"x": 0, "y": 0, "z": 0}
VINA_BOX_SIZE   = {"x": 20, "y": 20, "z": 20}


# ─── Docking Helper Functions ──────────────────────────────────────────────────

def _resolve_identifier(identifier: str, tmp_dir: str) -> str:
    """Return path to a local PDB file for *identifier*.

    If *identifier* matches the 4-char uppercase PDB-ID pattern the
    corresponding PDB file is fetched from RCSB and written to *tmp_dir*.
    Otherwise the string is returned verbatim (caller treats it as a local
    file path or raw FASTA content).
    """
    if PDB_ID_RE.match(identifier):
        url = RCSB_PDB_URL.format(pdb_id=identifier)
        dest = os.path.join(tmp_dir, f"{identifier}.pdb")
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        with open(dest, "wb") as fh:
            fh.write(r.content)
        return dest
    return identifier


def _to_pdbqt(input_path: str, output_path: str, is_receptor: bool) -> None:
    """Convert a PDB file at *input_path* to PDBQT using MGLTools scripts."""
    script = PREP_RECEPTOR if is_receptor else PREP_LIGAND
    if script is None:
        raise RuntimeError(
            "prepare_receptor4.py / prepare_ligand4.py not found on PATH.  "
            "Install MGLTools:  conda install -c bioconda mgltools"
        )
    cmd = [
        "python3", script,
        "-r" if is_receptor else "-l", input_path,
        "-o", output_path,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if result.returncode != 0:
        raise RuntimeError(
            f"{script} failed (exit {result.returncode}): {result.stderr[:500]}"
        )


def _run_vina(receptor_pdbqt: str, ligand_pdbqt: str, out_pdbqt: str) -> float:
    """Call the vina binary and return the best affinity score (kcal/mol)."""
    if VINA_BINARY is None:
        raise RuntimeError("AutoDock Vina not found on PATH")
    cmd = [
        VINA_BINARY,
        "--receptor",    receptor_pdbqt,
        "--ligand",      ligand_pdbqt,
        "--out",         out_pdbqt,
        "--center_x",    str(VINA_BOX_CENTER["x"]),
        "--center_y",    str(VINA_BOX_CENTER["y"]),
        "--center_z",    str(VINA_BOX_CENTER["z"]),
        "--size_x",      str(VINA_BOX_SIZE["x"]),
        "--size_y",      str(VINA_BOX_SIZE["y"]),
        "--size_z",      str(VINA_BOX_SIZE["z"]),
    ]
    logger.info("Running Vina: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    logger.debug("Vina stdout:\n%s", result.stdout)
    logger.debug("Vina stderr:\n%s", result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"Vina failed (exit {result.returncode}): {result.stderr[:500]}"
        )

    # Parse the affinity table — the first data row holds the best pose score.
    best_score: float | None = None
    for line in result.stdout.splitlines():
        stripped = line.strip()
        # Data rows start with a whitespace-padded integer (pose rank)
        parts = stripped.split()
        if len(parts) >= 2:
            try:
                score = float(parts[1])
            except ValueError:
                continue
            if best_score is None or score < best_score:
                best_score = score

    if best_score is None:
        raise RuntimeError(
            "Could not parse Vina output — no affinity table rows found"
        )
    return best_score


# ─── Docking Endpoints ─────────────────────────────────────────────────────────

@app.post("/api/docking")
async def run_docking(request: DockingRequest):
    """Run molecular docking between two proteins using AutoDock Vina.

    *protein_a* and *protein_b* may each be:
      - A 4-char PDB ID (e.g. ``"4G6F"``) — fetched from RCSB PDB
      - A raw PDB / PDBQT file path on disk
      - A FASTA string (PDBQT conversion will fail with an informative error)

    Returns ``{'status':'error', ...}`` if Vina or the PDBQT preparation
    scripts are missing, else ``{'status':'success', 'vina_score': ...}``.
    """
    if not isinstance(request.protein_a, str) or not request.protein_a.strip():
        raise HTTPException(status_code=422, detail="protein_a must be a non-empty string")
    if not isinstance(request.protein_b, str) or not request.protein_b.strip():
        raise HTTPException(status_code=422, detail="protein_b must be a non-empty string")

    # ── Guard: Vina not installed ──────────────────────────────────────────────
    if VINA_BINARY is None:
        return {
            "status":      "error",
            "message":     "AutoDock Vina not found on PATH",
            "protein_a":   request.protein_a,
            "protein_b":   request.protein_b,
            "p_value":     None,
            "r_squared":   None,
            "vina_score":  None,
        }

    tmp_dir: str | None = None
    try:
        tmp_dir = tempfile.mkdtemp(prefix="docking_")
        receptor_pdbqt = os.path.join(tmp_dir, "receptor.pdbqt")
        ligand_pdbqt   = os.path.join(tmp_dir, "ligand.pdbqt")
        out_pdbqt      = os.path.join(tmp_dir, "out.pdbqt")

        # Resolve protein_a → receptor PDB file path
        receptor_path = _resolve_identifier(request.protein_a.strip(), tmp_dir)
        logger.info("[docking] receptor=%r", receptor_path)

        # Resolve protein_b → ligand PDB file path
        ligand_path = _resolve_identifier(request.protein_b.strip(), tmp_dir)
        logger.info("[docking] ligand=%r", ligand_path)

        # Convert to PDBQT
        _to_pdbqt(receptor_path, receptor_pdbqt, is_receptor=True)
        _to_pdbqt(ligand_path,   ligand_pdbqt,   is_receptor=False)

        # Run Vina
        vina_score = _run_vina(receptor_pdbqt, ligand_pdbqt, out_pdbqt)

        return {
            "status":      "success",
            "protein_a":   request.protein_a,
            "protein_b":   request.protein_b,
            "p_value":     None,
            "r_squared":   None,
            "vina_score":  vina_score,
        }
    except HTTPException:
        raise
    except Exception as exc:
        logger.error("[docking] error: %s", exc, exc_info=True)
        return {
            "status":      "error",
            "message":     str(exc),
            "protein_a":   request.protein_a,
            "protein_b":   request.protein_b,
            "p_value":     None,
            "r_squared":   None,
            "vina_score":  None,
        }
    finally:
        if tmp_dir and os.path.isdir(tmp_dir):
            try:
                import shutil as _shutil
                _shutil.rmtree(tmp_dir, ignore_errors=True)
            except Exception:
                pass


@app.get("/api/docking/health")
async def docking_health():
    """Check AutoDock Vina installation status."""
    vina_ver = None
    if VINA_BINARY:
        try:
            r = subprocess.run([VINA_BINARY, "--version"], capture_output=True,
                               text=True, timeout=10)
            vina_ver = r.stdout.strip() or r.stderr.strip() or "unknown"
        except Exception:
            vina_ver = "unknown"
    return {
        "status":          "ok",
        "vina_installed":  VINA_BINARY is not None,
        "version":         vina_ver,
        "prep_receptor":   PREP_RECEPTOR,
        "prep_ligand":     PREP_LIGAND,
    }


# ─── Predict ──────────────────────────────────────────────────────────────────

@app.post("/api/predict")
async def predict_interaction(request: PredictRequest):
    """Predict interaction between two proteins."""
    if predictor is None or predictor.model is None:
        raise HTTPException(status_code=503, detail="Model not trained")

    protein_a = request.protein_a
    protein_b = request.protein_b

    if protein_a not in G.nodes() or protein_b not in G.nodes():
        raise HTTPException(status_code=404, detail="One or both proteins not found in network")

    predictions = predictor.predict([(protein_a, protein_b)])
    result = predictions.to_dict(orient="records")[0]

    # Check if edge already exists
    exists = G.has_edge(protein_a, protein_b)
    result["known_interaction"] = exists

    return result


@app.get("/api/protein/{protein_id}")
async def protein_detail(protein_id: str):
    """Get detailed info about a specific protein."""
    if protein_id not in G.nodes():
        raise HTTPException(status_code=404, detail="Protein not found")

    attrs = analysis_results.get("node_attributes", {}).get(protein_id, {})
    neighbors = list(G.neighbors(protein_id))
    neighbor_details = []
    for n in neighbors:
        n_attrs = analysis_results.get("node_attributes", {}).get(n, {})
        edge_data = G[protein_id][n] if G.has_edge(protein_id, n) else {}
        neighbor_details.append({
            "protein_id": n,
            "degree": n_attrs.get("degree", G.degree(n)),
            "community": n_attrs.get("community", -1),
            "score": edge_data.get("combined_score", 0),
        })

    return {
        "protein_id": protein_id,
        "degree": attrs.get("degree", G.degree(protein_id)),
        "community": attrs.get("community", -1),
        "degree_centrality": attrs.get("degree_centrality", 0),
        "betweenness_centrality": attrs.get("betweenness_centrality", 0),
        "num_neighbors": len(neighbors),
        "neighbors": neighbor_details,
    }


@app.get("/api/pipeline/run")
async def run_pipeline():
    """Re-run the full pipeline with synthetic data."""
    generate_synthetic_data()
    run_analysis()
    run_prediction()
    return {"status": "complete", "nodes": G.number_of_nodes(), "edges": G.number_of_edges()}


# ─── Serve React Frontend ───────────────────────────────────────────────────

frontend_build = PROJECT_ROOT / "frontend" / "dist"
if frontend_build.exists():
    from fastapi.responses import FileResponse, HTMLResponse
    import mimetypes

    # Ensure correct MIME types on Windows
    mimetypes.init()
    mimetypes.types_map[".js"] = "application/javascript"
    mimetypes.types_map[".mjs"] = "application/javascript"
    mimetypes.types_map[".css"] = "text/css"
    mimetypes.types_map[".html"] = "text/html"
    mimetypes.types_map[".svg"] = "image/svg+xml"
    mimetypes.types_map[".json"] = "application/json"
    mimetypes.types_map[".woff2"] = "font/woff2"
    mimetypes.types_map[".woff"] = "font/woff"
    mimetypes.types_map[".ttf"] = "font/ttf"

    @app.get("/assets/{file_path:path}", include_in_schema=False)
    async def serve_assets(file_path: str):
        file_on_disk = frontend_build / "assets" / file_path
        if not file_on_disk.exists():
            raise HTTPException(status_code=404, detail="Asset not found")
        ext = file_on_disk.suffix.lower()
        mime_map = {
            ".js": "application/javascript",
            ".mjs": "application/javascript",
            ".css": "text/css",
            ".html": "text/html",
            ".svg": "image/svg+xml",
            ".json": "application/json",
            ".woff2": "font/woff2",
            ".woff": "font/woff",
            ".ttf": "font/ttf",
            ".png": "image/png",
            ".jpg": "image/jpeg",
            ".jpeg": "image/jpeg",
            ".gif": "image/gif",
            ".ico": "image/x-icon",
        }
        media_type = mime_map.get(ext, "application/octet-stream")
        return FileResponse(str(file_on_disk), media_type=media_type)

    @app.get("/", include_in_schema=False)
    async def serve_index():
        return HTMLResponse(content=(frontend_build / "index.html").read_text(), status_code=200)

    @app.get("/{full_path:path}", include_in_schema=False)
    async def serve_react(full_path: str):
        if full_path.startswith("api/"):
            raise HTTPException(status_code=404, detail="API endpoint not found")
        return HTMLResponse(content=(frontend_build / "index.html").read_text(), status_code=200)


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000, reload=False)
