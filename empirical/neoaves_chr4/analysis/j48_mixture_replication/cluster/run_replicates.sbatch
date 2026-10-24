#!/bin/bash
#SBATCH --account=bio240015p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=14G
#SBATCH --time=24:00:00
#SBATCH --array=0-24%2
#SBATCH --job-name=j48rep
#SBATCH --output=/ocean/projects/cis260126p/shared/neoaves_stage4e/outputs/J48_MIX_REPLICATES/logs/%A_%a.out
#SBATCH --error=/ocean/projects/cis260126p/shared/neoaves_stage4e/outputs/J48_MIX_REPLICATES/logs/%A_%a.err
set -euo pipefail

BASE=/ocean/projects/cis260126p/shared/neoaves_stage4e
REPO="$BASE/repo"
OUT="$BASE/outputs/J48_MIX_REPLICATES"
INPUTS="$OUT/inputs"
RESULTS="$OUT/results"
MANIFEST="$REPO/empirical/neoaves_chr4/results/stage4e/j48_mixture_replication_manifest.tsv"
PYTHON_BIN="${PYTHON_BIN:-python3}"
ASTRAL_BIN="$BASE/software/astral4"
ID=$(awk -F '\t' -v n="$((SLURM_ARRAY_TASK_ID + 2))" 'NR==n {print $1}' "$MANIFEST")
[[ -n "$ID" ]] || { echo "manifest row missing for array index $SLURM_ARRAY_TASK_ID" >&2; exit 2; }
mkdir -p "$INPUTS" "$RESULTS" "$OUT/logs"

if [[ -d "$RESULTS/$ID" ]]; then
  "$PYTHON_BIN" "$REPO/empirical/neoaves_chr4/analysis/j48_mixture_replication/validate_results.py" \
    --root "$OUT" --manifest "$MANIFEST" --only "$ID" >/dev/null || { echo "existing output is incomplete or invalid; refusing overwrite: $ID" >&2; exit 3; }
  echo "validated existing result: $ID"; exit 0
fi

"$PYTHON_BIN" "$REPO/empirical/neoaves_chr4/analysis/j48_mixture_replication/prepare_manifest.py" \
  --manifest "$MANIFEST" --realize --job-index "$SLURM_ARRAY_TASK_ID" --input-root "$INPUTS"
"$PYTHON_BIN" "$REPO/empirical/neoaves_chr4/scripts/05_stage4e_run_aster.py" \
  --execute --allow-cluster-execution --mode topology_su --aster-bin "$ASTRAL_BIN" \
  --input "$INPUTS/$ID.tre" --output-dir "$RESULTS/$ID" --treatment "$ID" --threads 8 \
  --genelength 1000 --sulengthtype 1
"$PYTHON_BIN" "$REPO/empirical/neoaves_chr4/analysis/j48_mixture_replication/validate_results.py" \
  --root "$OUT" --manifest "$MANIFEST" --only "$ID" >/dev/null
echo "validated completed result: $ID"
