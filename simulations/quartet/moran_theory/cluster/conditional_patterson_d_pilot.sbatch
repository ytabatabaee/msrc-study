#!/bin/bash
# Prepared only; do not submit without reviewing the staged msrc-sim path.
# This job is for independent-history replication, not the small local pilot.
#SBATCH --account=bio240015p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --job-name=conditional-patterson-d
#SBATCH --output=logs/conditional_patterson_d_%A.out
#SBATCH --error=logs/conditional_patterson_d_%A.err

set -euo pipefail

MSRC_REPO="${MSRC_REPO:?Set MSRC_REPO to the staged repository root}"
MSRC_SIM="${MSRC_SIM:?Set MSRC_SIM to the staged msrc-sim checkout}"
PYTHON_BIN="${PYTHON_BIN:?Set PYTHON_BIN to the validated Python executable}"
OUT_DIR="${OUT_DIR:-$MSRC_REPO/simulations/quartet/moran_theory/results/conditional_patterson_d_bridges}"

test -x "$PYTHON_BIN"
test -f "$MSRC_REPO/simulations/quartet/moran_theory/run_conditional_patterson_d_pilot.py"
test -d "$MSRC_SIM/msrcsim"
mkdir -p "$OUT_DIR" logs

echo "host=$(hostname)"
echo "python=$PYTHON_BIN"
"$PYTHON_BIN" --version
echo "repo=$MSRC_REPO"
echo "msrc_sim=$MSRC_SIM"
echo "job=${SLURM_JOB_ID:-unsubmitted}"
date -u

"$PYTHON_BIN" -m simulations.quartet.moran_theory.run_conditional_patterson_d_pilot \
  --msrc-sim "$MSRC_SIM" \
  --output "$OUT_DIR" \
  --histories "${HISTORIES:-100}" \
  --genealogies "${GENEALOGIES:-1000}" \
  --M "${M:-20}" \
  --rho "${RHO:-0.2}" \
  --seed "${SEED:-20261028}"
