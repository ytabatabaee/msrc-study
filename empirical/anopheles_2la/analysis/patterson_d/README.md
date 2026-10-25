# Anopheles 2La Patterson's D workflow

This workflow computes Patterson's D directly from Ag3.10 nucleotide SNP
genotypes. It never uses Stage 2R quartet frequencies, local NJ trees, or
topology labels as substitutes for ABBA/BABA site patterns.

The design and numerical core are in
`scripts/run_patterson_d.py` and `scripts/patterson_d_core.py`. The default
command freezes the design and writes explicit pending tables. The eight-window
pilot is run only with `--pilot`; the 988-window analysis is run only with
`--run-full`. The isolated `--preflight-one-chunk` mode queries one frozen
1-Mb chunk and writes only to `results/patterson_d/bridges_preflight/`.

The validated Stage 1R environment is
`/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria310/bin/python` with
`malariagen-data 14.0.0`. The pilot queries each frozen 50-kb window
independently; it does not combine the eight distant windows into a 32-Mb
request. The completed pilot is recorded under `results/patterson_d/`.

Synthetic tests are independent of network access:

```bash
python -m pytest -q empirical/anopheles_2la/tests/test_patterson_d.py
```

After installing the repository's pinned MalariaGEN environment and restoring
credentials, run:

```bash
SSL_CERT_FILE=/path/to/certifi/cacert.pem \
REQUESTS_CA_BUNDLE=/path/to/certifi/cacert.pem \
python empirical/anopheles_2la/scripts/run_patterson_d.py --pilot \
  --url gs://vo_agam_release_master_us_central1
```

The full scan is deliberately a separate command:

```bash
python empirical/anopheles_2la/scripts/run_patterson_d.py --run-full \
  --url gs://vo_agam_release_master_us_central1 --chunk-bp 1000000
```

On Bridges-2, stage the repository and authenticated Python environment, run
the pilot first, inspect `pilot_results.tsv` and `pilot_report.md`, and submit
the full command only after reviewing outgroup callability and informative
site counts. The full mode checkpoints by completed window/design/filter and
uses the frozen Stage 2R grid. It does not redownload reads or alter Stage 1R
or Stage 2R outputs.

The exact Bridges-2 preflight submission command is:

```bash
export MSRC_REPO=/ocean/projects/cis260126p/shared/neoaves_stage4e/repo
export PYTHON_BIN=/path/to/authenticated/python
export AG3_URL=gs://vo_agam_release_master_us_central1
cd "$MSRC_REPO/empirical/anopheles_2la/analysis/patterson_d"
mkdir -p logs
test -f "$MSRC_REPO/empirical/anopheles_2la/scripts/run_patterson_d.py"
sbatch --export=ALL,MSRC_REPO="$MSRC_REPO",PYTHON_BIN="$PYTHON_BIN",AG3_URL="$AG3_URL" preflight_patterson_d.sbatch
```

The preflight is expected to produce `bridges_preflight/window_results.tsv`,
`manifest.json`, and `README.md`, containing 16 rows, verified sample order,
position bounds, access time, computation time, and peak RSS. It cannot write
the pilot or full result tables.

After reviewing the preflight, submit the full scan with:

```bash
export MSRC_REPO=/ocean/projects/cis260126p/shared/neoaves_stage4e/repo
export PYTHON_BIN=/path/to/authenticated/linux/python
export AG3_URL=gs://vo_agam_release_master_us_central1
cd "$MSRC_REPO/empirical/anopheles_2la/analysis/patterson_d"
export MSRC_REPO="$(cd ../../../.. && pwd)"
test -f "$MSRC_REPO/empirical/anopheles_2la/scripts/run_patterson_d.py"
mkdir -p logs
sbatch --export=ALL,MSRC_REPO="$MSRC_REPO",PYTHON_BIN="$PYTHON_BIN",AG3_URL="$AG3_URL",SSL_CERT_FILE="${SSL_CERT_FILE:-}",REQUESTS_CA_BUNDLE="${REQUESTS_CA_BUNDLE:-}" run_full_patterson_d.sbatch
```

The `mkdir -p logs` command must run before `sbatch`, because Slurm opens log
paths before executing the batch body. Monitor with `squeue -j JOBID` and
`sacct -j JOBID --format=JobID,State,Elapsed,MaxRSS,ExitCode`; the scripts print
host, interpreter, package version, working directory, commit, URL, job ID,
timestamp, and CA-bundle presence.

Run the pilot interactively first with the same `PYTHON_BIN`, certifi CA
variables, and Ag3 URL. The command above is the next full-scan command only
after reviewing the pilot; it has not been submitted.

The pilot completed with a corrected historical peak of 742.4 MiB (the old
macOS raw column was KiB despite its `peak_rss_mb` name) and about 8 minutes of
summed regional query time on the local authenticated client. New runs record
platform-normalized `peak_rss_mib`. The vectorized estimator agrees with the
scalar reference and gives a 21.7x synthetic CPU speedup. The full wrapper
requests 8 CPUs and 14 GB, which is 1.75 GB per core and respects the 2,000
MB/core RM-shared limit. Detailed preflight, checkpoint, design, and readiness
reports are under `results/patterson_d/`; the full scan remains unsubmitted.
