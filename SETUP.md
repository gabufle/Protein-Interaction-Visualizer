# Setup Guide

## Prerequisites

| Component | Version | Notes |
|-----------|---------|-------|
| Python | 3.10+ | Tested on 3.10–3.12 |
| pip | latest | `python -m ensurepip --upgrade` if missing |
| Git | any | To clone the repository |

---

## Installation

```bash
# 1. Clone the repo
cd protein-interaction-visualizer

# 2. Create and activate a virtual environment
python -m venv venv
source venv/bin/activate   # On Windows: venv\Scripts\activate

# 3. Install Python dependencies
pip install -r requirements.txt

# 4. Install the package in editable mode
pip install -e .
```

---

## Dependencies

### Python (CRITICAL - pinned in `requirements.txt`)

| Package | Purpose | Version pin |
|---------|---------|-------------|
| `numpy` | Numerical arrays | `>=1.24,<2.0` |
| `pandas` | Data frames | `>=2.0,<3.0` |
| `networkx` | Graph/network algorithms | `>=3.0,<4.0` |
| `scikit-learn` | ML predictor | `>=1.3,<2.0` |
| `plotly` / `matplotlib` / `seaborn` | Visualizations | various |
| `requests` | STRING/BioGRID API fetching | `>=2.30,<3.0` |
| `biopython` | Protein sequence parsing | `>=1.81,<2.0` |
| `fastapi` / `uvicorn` | REST API server | `>=0.100,<1.0` |
| `pyyaml` / `joblib` / `tqdm` | Config, serialization, progress | various |

### Python (development)

| Package | Purpose |
|---------|---------|
| `pytest` | Unit tests |
| `black` | Code formatting |
| `flake8` | Linting |
| `mypy` | Type checking |
| `jupyterlab` | Interactive notebooks |

---

## External Tool: AutoDock Vina (Phase 3 — integrated)

> **Current state:** The `/api/docking` endpoint performs **real** docking when
> `vina` is on `PATH`.  When Vina is absent (typical for fresh installs, since
> `sudo` is unavailable), both `/api/docking` and `/api/docking/health` return
> `{"status":"error","message":"AutoDock Vina not found on PATH"}` and
> `{"vina_installed":false}`, respectively — no stub values are returned.

### Install Vina

#### Linux / WSL / macOS (Homebrew)

```bash
# macOS
brew install autodock-vina

# Ubuntu / Debian / WSL  (may not be available in all repos)
sudo apt-get update
sudo apt-get install -y autodock-vina
```

#### From Source (works without sudo — installs to $HOME)

```bash
# 1. Download latest release:
#    https://github.com/ccsb-scripps/AutoDock-Vina/releases
cd ~/src                      # or any writable directory
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_src.tar.gz
tar xzf vina_1.2.5_src.tar.gz
cd vina_1.2.5_src

# 2. Build
make -j$(nproc)

# 3. Prepend build dir to PATH (add to ~/.bashrc / ~/.profile to make permanent)
export PATH="$HOME/src/vina_1.2.5_src:$PATH"
```

#### Verify

```bash
vina --version
# → AutoDock Vina 1.2.5  (or newer)
```

### MGLTools — PDBQT preparation scripts

Before Vina can run, PDB files must be converted to PDBQT format.  The
integration looks for `prepare_receptor4.py` and `prepare_ligand4.py` on `PATH`
(these ship with **MGLTools / AutoDockTools**).

```bash
# Conda (no sudo needed if you use --user or a user env)
conda install -c bioconda mgltools -y

# The scripts will then be found automatically:
#   prepare_receptor4.py  — converts receptor PDB → PDBQT
#   prepare_ligand4.py    — converts ligand   PDB → PDBQT
```

Without these scripts the docking endpoint reports an error before calling Vina.

### Vina CLI Reference

```bash
vina \
  --receptor protein_a.pdbqt \
  --ligand   protein_b.pdbqt \
  --out      docked.pdbqt \
  --center_x 0 --center_y 0 --center_z 0 \
  --size_x 20 --size_y 20 --size_z 20
```

---

## Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `STRING_API_KEY` | no (but recommended for real data) | STRING DB API key |
| `BIOGRID_USER_KEY` | no | BioGRID API key |

Create a `.env` file in the project root:

```bash
cp .env.example .env
# Edit .env and add your keys.
```

---

## Running the Backend

```bash
# Development (auto-reload)
uvicorn backend_server:app --host 0.0.0.0 --port 8000 --reload

# Or simply:
python backend_server.py

# Open http://localhost:8000
```

---

## Quick API Check

```bash
# General health
curl http://localhost:8000/api/health

# Docking stub health
curl http://localhost:8000/api/docking/health

# Docking stub request
curl -X POST http://localhost:8000/api/docking \
  -H "Content-Type: application/json" \
  -d '{"protein_a": "P04637", "protein_b": "P10279"}'
```
