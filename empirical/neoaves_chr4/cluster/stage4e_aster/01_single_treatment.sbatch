#!/usr/bin/env bash
# One full ASTRAL-IV inference per job.  Do not submit from this repository
# without first transferring inputs and setting the Bridges-2 paths.
#SBATCH --job-name=neoaves-aster-one
#SBATCH -p RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=16:00:00
# #SBATCH -A YOUR_ALLOCATION

set -euo pipefail
: "${MSRC_REPO:=$SLURM_SUBMIT_DIR}"
: "${NEOAVES_STAGE4E_DATA:=$MSRC_REPO/empirical/neoaves_chr4/data/stage4e}"
: "${NEOAVES_STAGE4E_OUT:=$PROJECT/neoaves_stage4e/outputs}"
: "${ASTER_BIN:?Set ASTER_BIN to the v1.25.4.8 astral4 executable}"
: "${PYTHON:=python3}"
: "${TREATMENT:?Export exactly one TREATMENT for this job}"

MANIFEST="$MSRC_REPO/empirical/neoaves_chr4/results/stage4e_input_manifest.tsv"
INPUT=$(awk -F '\t' -v t="$TREATMENT" 'NR>1 && $1==t {print $3}' "$MANIFEST")
if [[ -z "$INPUT" ]]; then echo "Unknown treatment: $TREATMENT" >&2; exit 2; fi
INPUT="$MSRC_REPO/$INPUT"
OUT="$NEOAVES_STAGE4E_OUT/$TREATMENT"
mkdir -p "$OUT"

exec "$PYTHON" "$MSRC_REPO/empirical/neoaves_chr4/scripts/05_stage4e_run_aster.py" \
  --execute --treatment "$TREATMENT" --aster-bin "$ASTER_BIN" \
  --input "$INPUT" --output-dir "$OUT" --mode topology_su \
  --genelength "${GENELENGTH:-1000}" --sulengthtype "${SULENGTHTYPE:-1}"
