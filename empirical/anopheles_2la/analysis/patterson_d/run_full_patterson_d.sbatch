#!/bin/bash
#SBATCH --job-name=anopheles-patterson-d
#SBATCH --account=bio240015p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=14G
#SBATCH --time=24:00:00
#SBATCH --output=logs/anopheles_patterson_d_%j.out
#SBATCH --error=logs/anopheles_patterson_d_%j.err
set -euo pipefail
MSRC_REPO="${MSRC_REPO:?set MSRC_REPO to the staged repository}"
PYTHON_BIN="${PYTHON_BIN:?set PYTHON_BIN to the Linux authenticated Python interpreter}"
cd "$MSRC_REPO"
test -f "$MSRC_REPO/empirical/anopheles_2la/scripts/run_patterson_d.py"
printf 'hostname=%s\n' "$(hostname)"
printf 'python_executable=%s\n' "$PYTHON_BIN"
"$PYTHON_BIN" --version
"$PYTHON_BIN" -c 'import importlib.metadata; print("malariagen_data_version=" + importlib.metadata.version("malariagen-data"))'
printf 'working_directory=%s\n' "$PWD"
printf 'commit_sha=%s\n' "$(git -C "$MSRC_REPO" rev-parse HEAD 2>/dev/null || printf unavailable)"
printf 'malariagen_url=%s\n' "${AG3_URL:-gs://vo_agam_release_master_us_central1}"
printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID:-unset}"
printf 'start_timestamp=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'SSL_CERT_FILE_set=%s REQUESTS_CA_BUNDLE_set=%s\n' "${SSL_CERT_FILE:+yes}" "${REQUESTS_CA_BUNDLE:+yes}"
"$PYTHON_BIN" empirical/anopheles_2la/scripts/run_patterson_d.py \
  --run-full --url "${AG3_URL:-gs://vo_agam_release_master_us_central1}" --chunk-bp 1000000
