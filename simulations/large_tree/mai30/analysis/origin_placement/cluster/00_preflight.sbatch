#!/bin/bash
# Bridges2-only dependency and allocation preflight. No simulator work.
#SBATCH --job-name=mai30-preflight
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
mkdir -p "${SCRIPT_DIR}/logs"
: "${MSRC_PYTHON:?set MSRC_PYTHON before submission}"
: "${MSRC_SIM_ROOT:?set MSRC_SIM_ROOT before submission}"
"$MSRC_PYTHON" "$SCRIPT_DIR/mai30_bridges2_driver.py" --mode preflight --stage forward
