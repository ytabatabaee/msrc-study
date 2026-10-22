#!/bin/bash
# Bridges-2 only: build one deterministic mixture, then run one ASTRAL-IV job.
# No submission or successor launch is performed by this script.
#SBATCH --job-name=neoaves-j48-mix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --array=0-12%2

set -euo pipefail
: "${MSRC_REPO:?Set MSRC_REPO to the repository checkout}"
: "${NEOAVES_STAGE4E_ASTRAL:?Set NEOAVES_STAGE4E_ASTRAL to the pinned astral4 executable}"
MANIFEST="$MSRC_REPO/empirical/neoaves_chr4/results/stage4e/j48_threshold_experiment_manifest.tsv"
PYTHON="${PYTHON:-python}"
mkdir -p logs
row="$($PYTHON - "$MANIFEST" "$SLURM_ARRAY_TASK_ID" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
print('\t'.join(rows[int(sys.argv[2])].values()))
PY
)"
IFS=$'\t' read -r mixture_id kind epsilon seed _rest <<< "$row"
if [[ "$kind" != "coarse_grid" ]]; then echo "array row is not a coarse-grid mixture" >&2; exit 2; fi
cd "$MSRC_REPO"
"$PYTHON" empirical/neoaves_chr4/scripts/05_stage4e_prepare_threshold_mixtures.py \
  --mixture-id "$mixture_id" --epsilon "$epsilon" --seed "$seed"
"$PYTHON" empirical/neoaves_chr4/scripts/05_stage4e_run_aster.py --execute \
  --treatment "$mixture_id" --aster-bin "$NEOAVES_STAGE4E_ASTRAL" \
  --input "empirical/neoaves_chr4/results/stage4e/threshold_mixtures/$mixture_id.tre" \
  --output-dir "empirical/neoaves_chr4/results/stage4e/cluster/$mixture_id" --mode both
