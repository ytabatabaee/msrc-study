# Anopheles 2La Patterson's D workflow

This workflow computes Patterson's D directly from Ag3.10 nucleotide SNP
genotypes. It never uses Stage 2R quartet frequencies, local NJ trees, or
topology labels as substitutes for ABBA/BABA site patterns.

The design and numerical core are in
`scripts/run_patterson_d.py` and `scripts/patterson_d_core.py`. The default
command freezes the design and writes explicit pending tables. The eight-window
pilot is run only with `--pilot`; the 988-window analysis is run only with
`--run-full`.

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

The explicit submission command is:

```bash
export MSRC_REPO=/ocean/projects/cis260126p/shared/neoaves_stage4e/repo
export PYTHON_BIN=/path/to/authenticated/python
export AG3_URL=gs://vo_agam_release_master_us_central1
cd "$MSRC_REPO/empirical/anopheles_2la/analysis/patterson_d"
mkdir -p logs
sbatch --export=ALL,MSRC_REPO="$MSRC_REPO",PYTHON_BIN="$PYTHON_BIN",AG3_URL="$AG3_URL" run_full_patterson_d.sbatch
```

Run the pilot interactively first with the same `PYTHON_BIN`, certifi CA
variables, and Ag3 URL. The command above is the next full-scan command only
after reviewing the pilot; it has not been submitted.

The pilot completed with approximately 0.77 GiB peak resident memory and about
8 minutes of summed regional query time on the local authenticated client.
The full wrapper requests 8 CPUs and 14 GB, which is 1.75 GB per core and
respects the 2,000 MB/core RM-shared limit. The full scan remains unsubmitted.
