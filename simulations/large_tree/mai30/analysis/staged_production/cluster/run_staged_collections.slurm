#!/bin/bash
# One task is one selected replicate/treatment collection: 90 tasks total.
# The array is prepared only; it has not been submitted from this workspace.
#SBATCH --job-name=mai30-msrc-stage
#SBATCH --array=0-89%6
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=mai30-stage-%A_%a.out
#SBATCH --error=mai30-stage-%A_%a.err

set -euo pipefail
STUDY_ROOT="${STUDY_ROOT:?set STUDY_ROOT to the msrc-study checkout}"
MSRC_PYTHON="${MSRC_PYTHON:?set MSRC_PYTHON to the msrc-sim environment Python}"
MANIFEST="${MANIFEST:-$STUDY_ROOT/simulations/large_tree/mai30/analysis/staged_production/staged_manifest.tsv}"
cd "$STUDY_ROOT"

"$MSRC_PYTHON" simulations/large_tree/mai30/analysis/scripts/run_production_collection.py \
  --manifest "$MANIFEST" --task-id "$SLURM_ARRAY_TASK_ID"
