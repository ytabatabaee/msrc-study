#!/bin/bash
# Prepared validation workflow; do not submit without reviewing the manifest
# and the staged msrc-sim path.
#SBATCH --account=bio240015p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --array=0-10%2
#SBATCH --job-name=pectinate-msrc
#SBATCH --output=logs/pectinate_%A_%a.out
#SBATCH --error=logs/pectinate_%A_%a.err

set -euo pipefail

MSRC_REPO="${MSRC_REPO:?Set MSRC_REPO to the staged repository root}"
MSRC_SIM="${MSRC_SIM:?Set MSRC_SIM to the staged msrc-sim checkout}"
PYTHON_BIN="${PYTHON_BIN:?Set PYTHON_BIN to the validated Python executable}"
MANIFEST="${MANIFEST:-$MSRC_REPO/simulations/quartet/moran_theory/cluster/unbalanced_validation_manifest.tsv}"
OUT_DIR="${OUT_DIR:-$MSRC_REPO/simulations/quartet/moran_theory/results/pectinate_validation}"

test -x "$PYTHON_BIN"
test -f "$MSRC_REPO/simulations/quartet/moran_theory/validate_unbalanced_species_tree.py"
test -d "$MSRC_SIM/msrcsim"
test -f "$MANIFEST"
mkdir -p "$OUT_DIR" logs

row="$(awk -F '\t' -v task="$SLURM_ARRAY_TASK_ID" 'NR == task + 2 {print}' "$MANIFEST")"
if [[ -z "$row" ]]; then
    echo "No manifest row for task $SLURM_ARRAY_TASK_ID" >&2
    exit 2
fi
IFS=$'\t' read -r task M fraction rho sister intermediate root selection replicates seed <<< "$row"

echo "host=$(hostname)"
echo "python=$PYTHON_BIN"
"$PYTHON_BIN" --version
echo "repo=$MSRC_REPO"
echo "msrc_sim=$MSRC_SIM"
echo "task=$task M=$M fraction=$fraction rho=$rho sister=$sister intermediate=$intermediate root=$root selection=$selection replicates=$replicates seed=$seed"

"$PYTHON_BIN" "$MSRC_REPO/simulations/quartet/moran_theory/validate_unbalanced_species_tree.py" \
  --msrc-sim "$MSRC_SIM" \
  --M "$M" --fraction "$fraction" --rho "$rho" \
  --sister-length "$sister" --intermediate-length "$intermediate" \
  --root-length "$root" --selection "$selection" \
  --replicates "$replicates" --seed "$seed" \
  --output "$OUT_DIR/task_${task}.tsv"
