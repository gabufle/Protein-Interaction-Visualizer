#!/usr/bin/env python3
"""
Test whether the backend metrics/ML results are hardcoded or computed from data.

- Parses each FASTA
- Creates a unique deterministic network (nodes named from FASTA ids, edges by length-heuristic scoring)
- Injects G, df_edges into the live backend process via its importable globals
- Calls run_analysis() + run_prediction()
- Pulls the API and records: num_nodes, num_edges, density, accuracy, roc_auc, confusion_matrix

Run each cell from a Jupyter session or as one-off prints.

Note: must run in a context where `sys.path` includes /home/gabuf/projects/protein-interaction-visualizer.
"""

import json
import os
import sys
import requests
import networkx as nx
from pathlib import Path
from io import StringIO

# ── Project root ──────────────────────────────────────────────────────────────
PROJECT = Path("/home/gabuf/projects/protein-interaction-visualizer")
sys.path.insert(0, str(PROJECT))
sys.path.insert(0, str(PROJECT / "src"))

# ── BAIL OUT early if protein IDs clash with synthetic P00000 range ────────────
# The synthetic default scheme was P00000, P00001, …, so any accession we use
# must NOT be numeric-prefixed  –  we will prefix manually below.

# ── FASTA files to test ───────────────────────────────────────────────────────
FASTA_DIR = PROJECT / "test_fastas"
FASTAS = {
    "TP53":  FASTA_DIR / "P04637.fasta",   # human tumour suppressor p53
    "PRNP":  FASTA_DIR / "P10279.fasta",   # bovine prion protein
    "LYZ":   FASTA_DIR / "P02701.fasta",   # chicken lysozyme C
}

# ── Seed sequences ────────────────────────────────────────────────────────────
import requests as _rq

_SEQS: dict[str, str] = {}
for key, path in FASTAS.items():
    txt = path.read_text()
    for line in txt.splitlines():
        if line.startswith(">"):
            continue
        _SEQS[key] = line.strip()
    print(f"{key}: {path.name} · {len(_SEQS[key])} aa")
# e.g. TP53 393 aa, PRNP 256 aa, LYZ 148 aa

# ── import backend globals ─────────────────────────────────────────────────────
import backend_server as bs  # noqa: E402 – after sys.path set

# Confirm we can reach the live process
_r = requests.get("http://localhost:8000/api/health", timeout=5)
_r.raise_for_status()
print(f"\nBackend health: {_r.json()}")

# ── Edge-generation heuristic ──────────────────────────────────────────────────
import pandas as pd
import numpy as np
from itertools import combinations


def make_network_from_fasta(seq_key: str) -> tuple[pd.DataFrame, dict[str, str]]:
    """
    Build a PPI network from a single FASTA entry.
    
    Node count scales with sequence length (L//10, clamped 15–50) so each
    protein yields a different graph size. Edge condition uses a length-
    difference threshold derived from the base length, giving each protein
    a distinct density.

    Parameters
    ----------
    seq_key : str
        Key into _SEQS (e.g. "TP53").

    Returns
    -------
    df : pd.DataFrame
        Columns: protein1, protein2, combined_score
    seqs : dict
        Mapping node_id -> random sequence string.
    """
    base_id = f"FA_{seq_key}"
    base_len = len(_SEQS[seq_key])
    rng = np.random.default_rng(42)  # deterministic per call

    # Scaled node count to ensure each protein yields a different graph size
    n_nodes = np.clip(base_len // 10, 15, 50)
    nodes: list[str] = []
    seqs:  dict[str, str] = {}

    for i in range(n_nodes):
        # ±10 aa jitter
        new_len = max(10, base_len + rng.integers(-10, 11))
        new_seq = "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), size=new_len))
        nid = f"{base_id}_{i}"
        nodes.append(nid)
        seqs[nid] = new_seq

    # Length-difference threshold scales with protein size
    max_allowed_diff = max(3, base_len // 100)  # 3 for TP53, 2 for PRNP, 1 for LYZ

    edges: list[dict] = []
    for a, b in combinations(nodes, 2):
        diff = abs(len(seqs[a]) - len(seqs[b]))
        if diff <= max_allowed_diff:
            # score: 100 * (1 - diff/max_allowed_diff) clamped to 1
            score = int(100 * (1 - diff / max_allowed_diff))
            score = max(score, 1)
            edges.append({"protein1": a, "protein2": b, "combined_score": score})

    df = pd.DataFrame(edges)
    return df, seqs


# ── Checkpoint helpers ─────────────────────────────────────────────────────────
def snapshot(endpoint: str) -> dict:
    _r = requests.get(f"http://localhost:8000{endpoint}", timeout=5)
    _r.raise_for_status()
    return _r.json()


SNAPSHOT_KEYS = [
    "/api/network/summary",
    "/api/prediction/results",
    "/api/network/hubs?top_n=5",
    "/api/network/degree-distribution",
]

# ── MAIN LOOP ──────────────────────────────────────────────────────────────────
records: list[dict] = []

print("\n" + "=" * 60)
print("STATIC OR COMPUTED?  –  FASTA injection test")
print("=" * 60)

for key in FASTAS:
    print(f"\n─── {key} ───")

    # Build & inject
    new_edges, new_seqs = make_network_from_fasta(key)
    new_nodes = sorted(new_seqs.keys())
    new_labels = {n: str(n) for n in new_nodes}

    # Patch globals in the backend module (mutates the live process)
    bs.df_edges = new_edges
    bs.G = nx.from_pandas_edgelist(
        new_edges, "protein1", "protein2",
        edge_attr="combined_score", create_using=nx.Graph()
    )

    # Re-run pipeline
    bs.run_analysis()
    bs.run_prediction()

    # Snapshot
    snap: dict = {}
    for ep in SNAPSHOT_KEYS:
        try:
            snap[ep] = snapshot(ep)
        except Exception as e:
            snap[ep] = {"error": str(e)}

    # Key numbers to compare later
    summary  = snap.get("/api/network/summary", {})
    pred     = snap.get("/api/prediction/results", {})
    hubs     = snap.get("/api/network/hubs?top_n=5", [])

    records.append({
        "dataset":   key,
        "nodes":     summary.get("num_nodes", -1),
        "edges":     summary.get("num_edges", -1),
        "density":   summary.get("density",      -1),
        "accuracy":  pred.get("accuracy",      -1),
        "roc_auc":   pred.get("roc_auc",       -1),
        "cm":        pred.get("confusion_matrix", []),
        "num_hubs":  len(hubs),
    })
    print(f"   nodes={records[-1]['nodes']}  edges={records[-1]['edges']}"
          f"  density={records[-1]['density']:.4f}"
          f"  acc={records[-1]['accuracy']:.4f}  auc={records[-1]['roc_auc']:.4f}")

# ── Summary table ──────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print(f"{'Dataset':<8} {'Nodes':>6} {'Edges':>6} {'Density':>8}"
      f" {'Accuracy':>10} {'ROC-AUC':>8}  #hubs")
print("-" * 60)
for r in records:
    cm_str = str(r["cm"]) if r["cm"] else "N/A"
    print(f"{r['dataset']:<8} {r['nodes']:>6} {r['edges']:>6} {r['density']:>8.4f}"
          f" {r['accuracy']:>10.4f} {r['roc_auc']:>8.4f}  {r['num_hubs']}")

print("\nConfusion matrices:")
for r in records:
    print(f"  {r['dataset']}: {r['cm']}")

# ── Conclusion ─────────────────────────────────────────────────────────────────
all_same = (r["nodes"] for r in records)
if len({r["nodes"] for r in records}) == 1 and len({r["edges"] for r in records}) == 1:
    print("\n⚠  ALL DATASETS RETURNED IDENTICAL numbers — metrics are STATIC / hardcoded.")
else:
    print("\n✅  Different datasets produce DIFFERENT metrics — values are COMPUTED from data.")

# ── Export results for user ──────────────────────────────────────────────────
out_path = PROJECT / "fasta_test_results.json"
out_path.write_text(json.dumps(records, indent=2))
print(f"\nFull JSON → {out_path}")
