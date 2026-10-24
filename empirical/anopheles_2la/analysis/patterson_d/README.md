# Anopheles 2La Patterson's D workflow

This workflow computes Patterson's D directly from Ag3.10 nucleotide SNP
genotypes. It never uses Stage 2R quartet frequencies, local NJ trees, or
topology labels as substitutes for ABBA/BABA site patterns.

The design and numerical core are in
`scripts/run_patterson_d.py` and `scripts/patterson_d_core.py`. The default
command freezes the design and writes explicit pending tables. The eight-window
pilot is run only with `--pilot`; the 988-window analysis is run only with
`--run-full`.

The current local environment lacks `malariagen_data`, so no real pilot was
run. Existing Stage 1A provenance records successful prior metadata access, but
the authenticated genotype-access environment must be restored before the
pilot command can query SNP arrays.

Synthetic tests are independent of network access:

```bash
python -m pytest -q empirical/anopheles_2la/tests/test_patterson_d.py
```

After installing the repository's pinned MalariaGEN environment and restoring
credentials, run:

```bash
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
sbatch --export=ALL,MSRC_REPO="$MSRC_REPO",PYTHON_BIN="$PYTHON_BIN",AG3_URL="$AG3_URL" run_full_patterson_d.sbatch
```

The pilot can be run interactively with the same `PYTHON_BIN` before that
submission. The initial wrapper requests 8 CPUs and 32 GB; there is no local
pilot runtime measurement because the required package is absent in the
current environment. Adjust the resource request only after the eight-window
pilot has reported actual peak usage.
