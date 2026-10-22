# Mai30 Bridges2 runbook

This workflow is prepared for Bridges2. No simulation, benchmark, genealogy generation, or expensive analysis was run locally. Every computational stage requires a live Slurm compute allocation and the pinned `msrc-sim` checkout. The adapters are repository-maintained; no `MAI_*_COMMAND` variables are needed.

## Checkout and configure

Use the existing PSC account, partition, modules, and Python environment from the project configuration. This repository intentionally does not encode an account name.

```bash
git clone <the-project-repository-url> msrc-study
git clone <the-msrc-sim-repository-url> msrc-sim
git -C msrc-sim checkout bfba58d627f6f467af649cdb0b37f969349e0aab
cd msrc-study/simulations/large_tree/mai30/analysis/origin_placement/cluster
export BRIDGES2_ACCOUNT='<account from the existing PSC configuration>'
export BRIDGES2_PARTITION='<partition from the existing PSC configuration>'
export MSRC_SIM_ROOT='<absolute path to the msrc-sim checkout>'
export MSRC_PYTHON='<configured Bridges2 Python executable>'
export MAI_CLUSTER_OUT='<persistent project filesystem output directory>'
```

Verify the pinned checkout and prepare the manifest without starting a simulator:

```bash
git -C "$MSRC_SIM_ROOT" rev-parse HEAD
"$MSRC_PYTHON" mai30_bridges2_driver.py --mode prepare-only
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" 00_preflight.sbatch
```

The preflight job is the first cluster command. It verifies `git`, `sbatch`, the configured interpreter, the exact simulator commit, adapter files, and a real running compute allocation.

## Explicit stage submission

Wait for each stage and inspect its output manifests before submitting the next stage. No script submits another job, retries with larger resources, or scales the history count.

```bash
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" bench_forward.sbatch
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" bench_backward.sbatch
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" 01_forward_history.sbatch
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" 02_panel_selection.sbatch
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" 02_forward_persistence.sbatch
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" 03_backward_genealogy.sbatch
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" 04_quartet_analysis.sbatch
sbatch --account="$BRIDGES2_ACCOUNT" --partition="$BRIDGES2_PARTITION" 05_independent_validation.sbatch
```

Submit `06_optional_astral.sbatch` only after validation and after recording the site ASTRAL executable/version. The optional stage combines validated 30-taxon gene-tree collections and writes a checked species tree.

Use these commands to inspect progress and resource usage:

```bash
squeue -u "$USER"
sacct -j JOBID --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,MaxVMSize,ExitCode
find "$MAI_CLUSTER_OUT" \( -name '*.done.json' -o -name 'output_manifest.json' \) -print | sort
```

Resubmit an interrupted array with the same command. Valid checkpoints are reused. A partial checkpoint is never counted complete; the adapter validates the compressed JSON, IDs, hashes, and required metadata before continuing. The authoritative `.done.json` marker is written only after the driver revalidates the output manifest and every recorded file checksum.

## Resource gates

Benchmark forward/backward tasks request 4 GB and 2 hours. Production forward and backward tasks request 8 GB and 24 hours. Quartet analysis requests 16 GB; all arrays are capped at two concurrent tasks. Persistence, panel, and validation are bounded single tasks. Review `sacct` elapsed time, MaxRSS, output size, completed history count, and genealogies per second before considering a separately reviewed 100-history or 1,000-history manifest.

The independent forward history is the replication unit. Each forward history is reused across all six rates: `0, 0.0001, 0.001, 0.002, 0.01, 0.02`.
