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

From the local repository root, use the established account and transfer only
the workflow code and frozen manifest. The large Stiller source tree is not
included; validate the checksum-identified copy already staged on Bridges-2
before proceeding:

```bash
export BRIDGES_HOST=bridges2.psc.edu PSC_USER=ytabatab
export REMOTE=/ocean/projects/cis260126p/shared/neoaves_stage4e
tar -cf - empirical/neoaves_chr4/analysis/j48_mixture_replication empirical/neoaves_chr4/results/stage4e/j48_mixture_replication_manifest.tsv empirical/neoaves_chr4/scripts/05_stage4e_run_aster.py empirical/neoaves_chr4/scripts/stage4d_core.py empirical/neoaves_chr4/scripts/stage4e_annotations.py \
  | ssh "$PSC_USER@$BRIDGES_HOST" "mkdir -p '$REMOTE/repo' '$REMOTE/outputs/J48_MIX_REPLICATES/logs' && tar -xf - -C '$REMOTE/repo'"
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
test -f empirical/neoaves_chr4/data/stage4d/named_recollapsed.tre
test -f empirical/neoaves_chr4/results/stage4d_locus_manifest.tsv
test -f empirical/neoaves_chr4/data/stage4d/jarvis48_taxa.tsv
"$PYTHON_BIN" empirical/neoaves_chr4/analysis/j48_mixture_replication/prepare_manifest.py --validate-sources
"$PYTHON_BIN" -m py_compile empirical/neoaves_chr4/analysis/j48_mixture_replication/*.py empirical/neoaves_chr4/scripts/05_stage4e_run_aster.py
"$PYTHON_BIN" -c 'import dendropy; print(dendropy.__version__)'
"$ASTRAL_BIN" -h > "$REMOTE/outputs/J48_MIX_REPLICATES/astral_help.txt" 2>&1
scontrol show partition RM-shared
sacctmgr show qos rm format=Name,MaxTRES,MaxWall
```

Inspect the requested memory against the current Bridges-2 RM-shared policy
and compare with prior successful jobs before submission. The committed script
requests the established 8 CPUs, 14 GB total memory, and two concurrent array
tasks. Submit only array task 0 using an SLURM command-line override; do not
edit or regenerate the script:

```bash
mkdir -p "$REMOTE/outputs/J48_MIX_REPLICATES/logs"
TEST_JOB=$(sbatch --parsable --array=0-0%1 --export=ALL,PYTHON_BIN="$PYTHON_BIN" empirical/neoaves_chr4/analysis/j48_mixture_replication/cluster/run_replicates.sbatch)
echo "$TEST_JOB"
squeue -u "$USER"
sacct -j "$TEST_JOB" --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
"$PYTHON_BIN" empirical/neoaves_chr4/analysis/j48_mixture_replication/validate_results.py --root "$REMOTE/outputs/J48_MIX_REPLICATES" --manifest empirical/neoaves_chr4/results/stage4e/j48_mixture_replication_manifest.tsv --only J48_MIX_REP_E0p020_R01
```

After the test result passes the validator and resource review, submit the
remaining array tasks. The script skips a validated result and refuses to
overwrite an existing incomplete or mismatched result:

```bash
ARRAY_JOB=$(sbatch --parsable --array=1-24%2 --export=ALL,PYTHON_BIN="$PYTHON_BIN" empirical/neoaves_chr4/analysis/j48_mixture_replication/cluster/run_replicates.sbatch)
echo "$ARRAY_JOB"
squeue -u "$USER"
sacct -j "$ARRAY_JOB" --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
```

Validate all returned results before analysis:

```bash
"$PYTHON_BIN" empirical/neoaves_chr4/analysis/j48_mixture_replication/validate_results.py \
  --root "$REMOTE/outputs/J48_MIX_REPLICATES" \
  --manifest empirical/neoaves_chr4/results/stage4e/j48_mixture_replication_manifest.tsv
```

Return the result tree, sidecars, logs, and manifest with tar over SSH:

```bash
ssh "$PSC_USER@$BRIDGES_HOST" "tar -cf - -C '$REMOTE/outputs' J48_MIX_REPLICATES" \
  | tar -xf - -C empirical/neoaves_chr4/results/stage4e
```

Run local validation and postprocessing only after the archive is returned.
Then review the generated tables and commit locally. Do not call `sbatch`
until the test task and resource checks have been reviewed.
