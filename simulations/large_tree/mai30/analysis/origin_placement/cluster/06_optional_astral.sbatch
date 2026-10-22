#!/bin/bash
# Optional, validation-gated full-tree analysis. Never submitted automatically.
#SBATCH --job-name=mai30-astral
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
mkdir -p "${SCRIPT_DIR}/logs"
: "${MSRC_PYTHON:?set MSRC_PYTHON before submission}"
: "${ASTRAL:?set ASTRAL explicitly; no local default is allowed}"
"$MSRC_PYTHON" "$SCRIPT_DIR/mai30_bridges2_driver.py" --mode analyze-existing --stage astral --task-id 0 --output-dir "${MAI_CLUSTER_OUT:?set MAI_CLUSTER_OUT}"
