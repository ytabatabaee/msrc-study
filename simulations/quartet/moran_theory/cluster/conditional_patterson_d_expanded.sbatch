#!/bin/bash
# Prepared only; do not submit automatically.
# Ten independent array tasks × ten histories × 5,000 genealogies = 100 histories.
# Create the Slurm log directory before submission because Slurm opens log files first.
# Example: mkdir -p logs && sbatch --array=0-9%2 conditional_patterson_d_expanded.sbatch
#SBATCH --account=bio240015p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=08:00:00
#SBATCH --job-name=conditional-patterson-d
#SBATCH --array=0-9%2
#SBATCH --output=logs/conditional_patterson_d_%A_%a.out
#SBATCH --error=logs/conditional_patterson_d_%A_%a.err

set -euo pipefail

MSRC_REPO="${MSRC_REPO:?Set MSRC_REPO to the staged repository root}"
MSRC_SIM="${MSRC_SIM:?Set MSRC_SIM to the staged msrc-sim checkout}"
PYTHON_BIN="${PYTHON_BIN:?Set PYTHON_BIN to the validated Python executable}"
OUT_ROOT="${OUT_ROOT:-$MSRC_REPO/simulations/quartet/moran_theory/results/conditional_patterson_d_bridges}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array}"
TASK_OUT="$OUT_ROOT/task_${TASK_ID}"
HISTORIES_PER_TASK="${HISTORIES_PER_TASK:-10}"
GENEALOGIES="${GENEALOGIES:-5000}"
M="${M:-20}"
RHO="${RHO:-0.2}"
BASE_SEED="${BASE_SEED:-20261028}"
SEED=$((BASE_SEED + TASK_ID))

test -x "$PYTHON_BIN"
test -f "$MSRC_REPO/simulations/quartet/moran_theory/run_conditional_patterson_d_pilot.py"
test -d "$MSRC_SIM/msrcsim"
mkdir -p "$TASK_OUT"

echo "host=$(hostname)"
echo "python=$PYTHON_BIN"
"$PYTHON_BIN" --version
echo "repo=$MSRC_REPO"
echo "msrc_sim=$MSRC_SIM"
echo "job=${SLURM_JOB_ID:-unsubmitted} array_task=$TASK_ID"
echo "start=$(date -u +%Y-%m-%dT%H:%M:%SZ)"

"$PYTHON_BIN" -m simulations.quartet.moran_theory.run_conditional_patterson_d_pilot \
  --msrc-sim "$MSRC_SIM" \
  --output "$TASK_OUT" \
  --histories "$HISTORIES_PER_TASK" \
  --genealogies "$GENEALOGIES" \
  --M "$M" \
  --rho "$RHO" \
  --seed "$SEED"
