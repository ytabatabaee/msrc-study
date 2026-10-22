# Mai30 Bridges2 quartet-validation runbook

This workflow is prepared for Bridges2 and deliberately refuses to run outside
a Slurm compute allocation. No simulation, benchmark, or stochastic analysis
was run while preparing these files. The optimized simulator checkout is not
present in this workspace; Bridges2 preflight therefore requires the exact
`msrc-sim` commit `bfba58d627f6f467af649cdb0b37f969349e0aab`.

## Configure the allocation and environment

Use the existing Bridges2 project configuration if it is available in the
user's shell or site module files. This repository contains no account or
partition name, so do not replace the variables below with guessed values.

```bash
cd /path/to/msrc-study/simulations/large_tree/mai30/analysis/origin_placement/cluster
export BRIDGES2_ACCOUNT='the-existing-account-from-your-Bridges2-config'
export BRIDGES2_PARTITION='the-existing-partition-from-your-Bridges2-config'
export MSRC_SIM_ROOT='/path/on/Bridges2/to/msrc-sim'
export MSRC_PYTHON='/path/on/Bridges2/to/the/configured/python'
export MAI_CLUSTER_OUT='/path/on/Bridges2/to/mai30-origin-bridges2'
```

`MSRC_SIM_ROOT` must be a Git checkout at the required commit. Each stage also
requires an explicit command variable, such as `MAI_FORWARD_COMMAND` or
`MAI_BACKWARD_COMMAND`. Set those only after checking the optimized checkout's
CLI/API and adapting the stage adapter to write per-history checkpoints. The
wrapper passes `MAI_STAGE`, `MAI_TASK_ID`, `MAI_TASK_JSON`, and
`MAI_OUTPUT_DIR` to the command. It does not guess an unavailable simulator
interface.

The first command is the non-computational manifest check:

```bash
"$MSRC_PYTHON" mai30_bridges2_driver.py --mode prepare-only
```

Then submit the allocation preflight, passing the account and partition at
submission time:

```bash
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" 00_preflight.sbatch
```

## Staged submission

Submit each stage only after inspecting the prior stage's completion markers.
The benchmark arrays are intentionally capped at two concurrent tasks and do
not expand automatically.

```bash
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" bench_forward.sbatch
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" bench_backward.sbatch
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" 01_forward_history.sbatch
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" 02_panel_selection.sbatch
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" 03_backward_genealogy.sbatch
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" 04_quartet_analysis.sbatch
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" 05_independent_validation.sbatch
```

Do not submit `06_optional_astral.sbatch` until independent validation has
passed and the original ASTRAL command, version, and input contract have been
recorded. It is an optional extension and is never submitted by another job.

Inspect jobs and resource usage with:

```bash
squeue -u "$USER"
sacct -j JOBID --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,MaxVMSize,ExitCode
find "$MAI_CLUSTER_OUT" -name '*.done.json' -print | sort
```

Interrupted tasks can be resubmitted with the same array command. Existing
valid completion markers are skipped; a marker is written atomically only
after the configured command returns zero. There is no automatic retry with
larger resources and no automatic scale-up from 100 to 1,000 histories.

## Resource gates

The initial limits are one CPU per simulation task, 8 GiB and 24 hours for
forward or backward production tasks, and `%2` array concurrency. Benchmarks
use one CPU, 4 GiB, 2 hours, and the same `%2` cap. Analysis uses 2 CPUs and
16 GiB; optional ASTRAL uses 4 CPUs and 16 GiB for 8 hours. Increase resources
only after reviewing `sacct` peak RSS, elapsed time, output size, completed
history count, and genealogies per second from the benchmark logs. Estimate
100- and 1,000-history costs from those measurements, then submit a revised
manifest explicitly; do not edit the current production scale in place.

The scientific replication unit remains an independent forward history. Each
forward history is reused for all six prescribed rates
`0, 0.0001, 0.001, 0.002, 0.01, 0.02`. Scenarios A and B are separate
configuration rows. The seed policy is recorded in the manifest and passed to
the stage adapter through `MAI_TASK_JSON`.
