#!/bin/bash
# Static-only checks. This script never imports or invokes msrc-sim.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
for file in 00_preflight.sbatch 01_forward_history.sbatch 02_panel_selection.sbatch 03_backward_genealogy.sbatch 04_quartet_analysis.sbatch 05_independent_validation.sbatch 06_optional_astral.sbatch bench_forward.sbatch bench_backward.sbatch; do
  test -f "$HERE/$file"
  bash -n "$HERE/$file"
  grep -q -- '--mem=' "$HERE/$file"
  grep -q -- '--time=' "$HERE/$file"
  grep -q -- '--output=' "$HERE/$file"
  grep -q -- '--error=' "$HERE/$file"
done
grep -q -- '#SBATCH --array=0-7%2' "$HERE/01_forward_history.sbatch"
grep -q -- '#SBATCH --array=0-47%2' "$HERE/03_backward_genealogy.sbatch"
grep -q -- 'EXPECTED_SIMULATOR_COMMIT' "$HERE/mai30_bridges2_driver.py"
grep -q -- 'SLURM_JOB_ID' "$HERE/mai30_bridges2_driver.py"
grep -q -- 'prepare-only' "$HERE/mai30_bridges2_driver.py"
if grep -q '^import msrcsim\|^from msrcsim' "$HERE/mai30_bridges2_driver.py"; then
  echo 'driver imports simulator at module load' >&2
  exit 1
fi
awk -F '\t' 'NR == 1 { next } { count[$1]++ } END { if (count["forward"] != 8 || count["backward"] != 48 || count["benchmark_forward"] != 8 || count["benchmark_backward"] != 8) exit 1 }' "$HERE/mai30_bridges2_manifest.tsv"
grep -q $'\t0\t100\t' "$HERE/mai30_bridges2_manifest.tsv"
grep -q $'\t0.0001\t100\t' "$HERE/mai30_bridges2_manifest.tsv"
grep -q $'\t0.001\t100\t' "$HERE/mai30_bridges2_manifest.tsv"
grep -q $'\t0.002\t100\t' "$HERE/mai30_bridges2_manifest.tsv"
grep -q $'\t0.01\t100\t' "$HERE/mai30_bridges2_manifest.tsv"
grep -q $'\t0.02\t100\t' "$HERE/mai30_bridges2_manifest.tsv"
echo 'static workflow checks passed'
