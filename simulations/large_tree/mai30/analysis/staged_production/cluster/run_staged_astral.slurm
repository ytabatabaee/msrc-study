#!/bin/bash
# Submit only after collection_complete.json exists for every staged task.
#SBATCH --job-name=mai30-astral-stage
#SBATCH --array=0-89%4
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=mai30-stage-astral-%A_%a.out
#SBATCH --error=mai30-stage-astral-%A_%a.err

set -euo pipefail
STUDY_ROOT="${STUDY_ROOT:?set STUDY_ROOT to the msrc-study checkout}"
MSRC_PYTHON="${MSRC_PYTHON:?set MSRC_PYTHON to the msrc-sim environment Python}"
MANIFEST="${MANIFEST:-$STUDY_ROOT/simulations/large_tree/mai30/analysis/staged_production/staged_manifest.tsv}"
ASTRAL="${ASTRAL:?set ASTRAL to the ASTRAL-IV executable}"
cd "$STUDY_ROOT"

"$MSRC_PYTHON" simulations/large_tree/mai30/analysis/scripts/run_production_astral.py \
  --manifest "$MANIFEST" --task-id "$SLURM_ARRAY_TASK_ID" --astral "$ASTRAL"
