#!/bin/bash
# One backward genealogy from each benchmark forward history; submit explicitly.
#SBATCH --job-name=mai30-bench-bwd
#SBATCH --array=0-7%2
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
mkdir -p "${SCRIPT_DIR}/logs"
: "${MSRC_PYTHON:?set MSRC_PYTHON before submission}"
"$MSRC_PYTHON" "$SCRIPT_DIR/mai30_bridges2_driver.py" --mode run --stage benchmark_backward --task-id "$SLURM_ARRAY_TASK_ID" --output-dir "${MAI_CLUSTER_OUT:?set MAI_CLUSTER_OUT}"
