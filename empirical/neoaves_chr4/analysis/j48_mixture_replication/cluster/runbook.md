# J48 replication runbook

This workflow prepares 25 additional mixtures and is ready for Bridges-2. It
does not submit jobs. The original 13 mixture outputs remain under the local
repository's `results/stage4e/cluster/J48_MIX_E*` directories and are not
inputs to this array.

The remote layout is:

```text
/ocean/projects/cis260126p/shared/neoaves_stage4e/
  repo/                         staged code and frozen manifests
  software/astral4              pinned executable
  outputs/J48_MIX_REPLICATES/
    inputs/                      generated gene-tree mixtures and sidecars
    results/                     ASTRAL outputs
    logs/                        SLURM logs
```

Set `BRIDGES_HOST` and `PSC_USER` locally. Transfer with SSH and tar:

```bash
export BRIDGES_HOST=bridges2.psc.edu PSC_USER=YOUR_PSC_USERNAME
export REMOTE=/ocean/projects/cis260126p/shared/neoaves_stage4e
tar -czf /tmp/j48-replication-code.tgz empirical/neoaves_chr4/analysis/j48_mixture_replication empirical/neoaves_chr4/scripts/05_stage4e_run_aster.py empirical/neoaves_chr4/scripts/stage4d_core.py empirical/neoaves_chr4/scripts/stage4e_annotations.py empirical/neoaves_chr4/data/stage4d empirical/neoaves_chr4/results/stage4d_locus_manifest.tsv
ssh "$PSC_USER@$BRIDGES_HOST" "mkdir -p '$REMOTE/repo' '$REMOTE/outputs/J48_MIX_REPLICATES/logs'"
scp /tmp/j48-replication-code.tgz "$PSC_USER@$BRIDGES_HOST:$REMOTE/"
ssh "$PSC_USER@$BRIDGES_HOST" "tar -xzf '$REMOTE/j48-replication-code.tgz' -C '$REMOTE/repo'"
```

On the cluster, load the site Python environment that provides DendroPy and
run the manifest preflight. Set the remote paths and `PYTHON_BIN` in that
session:

```bash
export REMOTE=/ocean/projects/cis260126p/shared/neoaves_stage4e
cd "$REMOTE/repo"
export PYTHON_BIN=python3
export ASTRAL_BIN="$REMOTE/software/astral4"
test -x "$ASTRAL_BIN"
"$PYTHON_BIN" empirical/neoaves_chr4/analysis/j48_mixture_replication/prepare_manifest.py --write-manifest
"$PYTHON_BIN" -m py_compile empirical/neoaves_chr4/analysis/j48_mixture_replication/*.py empirical/neoaves_chr4/scripts/05_stage4e_run_aster.py
"$ASTRAL_BIN" -h > "$REMOTE/outputs/J48_MIX_REPLICATES/astral_help.txt" 2>&1
```

Inspect the requested memory against the current Bridges-2 RM-shared policy
and compare with prior successful jobs before submission. Submit one test task
by using a temporary copy of the array script with `#SBATCH --array=0-0`:

```bash
mkdir -p "$REMOTE/outputs/J48_MIX_REPLICATES/logs"
sed 's/^#SBATCH --array=0-24%2$/#SBATCH --array=0-0/' empirical/neoaves_chr4/analysis/j48_mixture_replication/cluster/run_replicates.sbatch > /tmp/j48_rep_test.sbatch
sbatch /tmp/j48_rep_test.sbatch
squeue -u "$USER"
sacct -j TEST_JOB_ID --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
```

After the test result passes the validator and resource review, submit the
remaining array tasks. The script skips a validated result and refuses to
overwrite an existing incomplete or mismatched result:

```bash
sbatch empirical/neoaves_chr4/analysis/j48_mixture_replication/cluster/run_replicates.sbatch
squeue -u "$USER"
sacct -j ARRAY_JOB_ID --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
```

Validate all returned results before analysis:

```bash
"$PYTHON_BIN" empirical/neoaves_chr4/analysis/j48_mixture_replication/validate_results.py \
  --root "$REMOTE/outputs/J48_MIX_REPLICATES" \
  --manifest empirical/neoaves_chr4/results/stage4e/j48_mixture_replication_manifest.tsv
```

Return the result tree, sidecars, logs, and manifest with tar over SSH:

```bash
ssh "$PSC_USER@$BRIDGES_HOST" "tar -czf - -C '$REMOTE/outputs' J48_MIX_REPLICATES" > /tmp/J48_MIX_REPLICATES.tgz
tar -xzf /tmp/J48_MIX_REPLICATES.tgz -C empirical/neoaves_chr4/results/stage4e
```

Run local validation and postprocessing only after the archive is returned.
Then review the generated tables and commit locally. Do not call `sbatch`
until the test task and resource checks have been reviewed.
