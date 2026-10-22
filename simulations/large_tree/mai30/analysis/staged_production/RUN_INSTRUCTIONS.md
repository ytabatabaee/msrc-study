# Staged Mai30 run instructions

The staged manifest is
`analysis/staged_production/staged_manifest.tsv`.  It has 90 tasks and writes
to `datasets/processed/mai30_staged_runs`.  It was selected from all 300
existing original Mai gene-tree collections using the original discordance
summary only: sort each condition by mean normalized RF and select the two
endpoints of each of five equal strata.  The validated four-replicate summary
was not overwritten; only replicates 042 and 007 happen to be selected by this
rule.

Before submission on Bridges2, set real paths for the checkout and the
authoritative `msrc-sim` Python environment:

```bash
export STUDY_ROOT=/path/to/msrc-study
export MSRC_PYTHON=/path/to/msrc-sim/.venv/bin/python
cd "$STUDY_ROOT"
```

Run a one-history benchmark for one task from each condition:

```bash
python simulations/large_tree/mai30/analysis/staged_production/scripts/profile_runtime.py \
  --manifest simulations/large_tree/mai30/analysis/staged_production/staged_manifest.tsv \
  --task-id 0 --histories 1
```

Use the existing environment's configured Python in place of `python` on the
cluster.  Repeat with a task from the 1.5 and 5 conditions, record walltime,
RSS, and storage, then adjust only SLURM resource requests if needed.

Submit the validated collection array:

```bash
sbatch --export=ALL,STUDY_ROOT="$STUDY_ROOT",MSRC_PYTHON="$MSRC_PYTHON" \
  simulations/large_tree/mai30/analysis/staged_production/cluster/run_staged_collections.slurm
```

Only validated `collection_complete.json` markers allow ASTRAL tasks to run.
After checking collection completion counts, submit:

```bash
sbatch --export=ALL,STUDY_ROOT="$STUDY_ROOT",ASTRAL=/path/to/astral \
  simulations/large_tree/mai30/analysis/staged_production/cluster/run_staged_astral.slurm
```

The full 900-collection production array is deliberately outside this task.
