#!/usr/bin/env bash
# download_pdb.sh
#
# Downloads a protein structure file (PDB format) from the RCSB PDB archive
# into the project's temp/ directory.
#
# Usage:   ./download_pdb.sh <PDB_ID>
# Example: ./download_pdb.sh 1crn
#
# Requires: curl, mkdir

set -euo pipefail

# ── Helpers ───────────────────────────────────────────────────────────────────

usage() {
  cat <<'EOF'
Usage:  ./download_pdb.sh <PDB_ID>

  <PDB_ID>  4-character PDB accession code, e.g. 1crn, 6isb …

Downloads the file into:  temp/<PDB_ID>.pdb
EOF
  exit 1
}

# ── Argument validation ───────────────────────────────────────────────────────

if [[ $# -lt 1 ]]; then
  echo "Error: missing PDB ID argument." >&2
  usage
fi

PDB_ID="${1^^}"   # uppercase the ID

if [[ ! "$PDB_ID" =~ ^[0-9A-Z]{4}$ ]]; then
  echo "Error: '$PDB_ID' does not look like a valid 4-character PDB ID." >&2
  exit 1
fi

# ── Download ──────────────────────────────────────────────────────────────────

RCSB_URL="https://files.rcsb.org/download/${PDB_ID}.pdb"
OUT_DIR="$(dirname "$0")/temp"
OUT_FILE="${OUT_DIR}/${PDB_ID}.pdb"

mkdir -p "$OUT_DIR"

echo "Fetching ${PDB_ID} from ${RCSB_URL} …"
curl -fL --retry 3 --retry-delay 2 -o "$OUT_FILE" "$RCSB_URL"

echo "Saved to: ${OUT_FILE}"
