#!/bin/bash
#SBATCH --job-name=mai30-forward
#SBATCH --array=0-7%2
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
mkdir -p "${SCRIPT_DIR}/logs"
: "${MSRC_PYTHON:?set MSRC_PYTHON before submission}"
"$MSRC_PYTHON" "$SCRIPT_DIR/mai30_bridges2_driver.py" --mode run --stage forward --task-id "$SLURM_ARRAY_TASK_ID" --output-dir "${MAI_CLUSTER_OUT:?set MAI_CLUSTER_OUT}"
