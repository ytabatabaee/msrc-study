#!/bin/bash
#SBATCH --job-name=anopheles-patterson-d
#SBATCH --account=bio240015p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/anopheles_patterson_d_%j.out
#SBATCH --error=logs/anopheles_patterson_d_%j.err
set -euo pipefail
cd "${MSRC_REPO:?set MSRC_REPO to the staged repository}"
mkdir -p logs
"${PYTHON_BIN:-python3}" empirical/anopheles_2la/scripts/run_patterson_d.py \
  --run-full --url "${AG3_URL:-gs://vo_agam_release_master_us_central1}" --chunk-bp 1000000
