#!/bin/bash
# Existing-checkpoint persistence summary; no simulator import is required.
#SBATCH --job-name=mai30-persistence
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
mkdir -p "${SCRIPT_DIR}/logs"
: "${MSRC_PYTHON:?set MSRC_PYTHON before submission}"
"$MSRC_PYTHON" "$SCRIPT_DIR/mai30_bridges2_driver.py" --mode analyze-existing --stage persistence --task-id 0 --output-dir "${MAI_CLUSTER_OUT:?set MAI_CLUSTER_OUT}"
