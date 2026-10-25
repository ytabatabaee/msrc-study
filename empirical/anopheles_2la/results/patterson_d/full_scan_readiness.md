# Full-scan readiness

Decision: **READY_FOR_FULL_SCAN**, subject to a remote first-chunk preflight.

The code has a validated scalar reference, a vectorized production estimator,
atomic checkpointing, explicit row compatibility checks, frozen 988-window
input checksums, and bounded contiguous query regions. The full command is
deliberately not executed in this audit.

The recommended initial allocation remains 8 CPUs and 14 GB on
`RM-shared` under `bio240015p` (1.75 GB/core, below the 2,000 MB/core limit).
The historical pilot used about 742.4 MiB peak RSS. The synthetic vectorized
benchmark used 96.10 MiB; the allocation leaves substantial memory margin,
but remote array chunking and library overhead must be checked in the first
1-Mb query. Query wall time is the dominant unknown and cannot be inferred
reliably from CPU-only benchmarking.

From the verified Bridges-2 checkout, after checking the authenticated Linux
environment, run the isolated preflight first:

```bash
export MSRC_REPO=/ocean/projects/cis260126p/shared/neoaves_stage4e/repo
export PYTHON_BIN=/path/to/authenticated/linux/python
export AG3_URL=gs://vo_agam_release_master_us_central1
cd "$MSRC_REPO/empirical/anopheles_2la/analysis/patterson_d"
export MSRC_REPO="$(cd ../../../.. && pwd)"
test -f "$MSRC_REPO/empirical/anopheles_2la/scripts/run_patterson_d.py"
mkdir -p logs
sbatch --export=ALL,MSRC_REPO="$MSRC_REPO",PYTHON_BIN="$PYTHON_BIN",AG3_URL="$AG3_URL",SSL_CERT_FILE="${SSL_CERT_FILE:-}",REQUESTS_CA_BUNDLE="${REQUESTS_CA_BUNDLE:-}" preflight_patterson_d.sbatch
```

Inspect `bridges_preflight/manifest.json` and its 16-row TSV. Then use:

```bash
test -f "$MSRC_REPO/empirical/anopheles_2la/scripts/run_patterson_d.py"
mkdir -p logs
sbatch --export=ALL,MSRC_REPO="$MSRC_REPO",PYTHON_BIN="$PYTHON_BIN",AG3_URL="$AG3_URL",SSL_CERT_FILE="${SSL_CERT_FILE:-}",REQUESTS_CA_BUNDLE="${REQUESTS_CA_BUNDLE:-}" run_full_patterson_d.sbatch
```

The full branch writes `full_scan_report.md` and never overwrites the frozen
pilot report or pilot results. No submission was made during this audit.
