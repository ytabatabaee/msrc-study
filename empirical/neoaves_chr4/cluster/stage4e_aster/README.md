# Stage 4E ASTRAL-IV Bridges-2 execution

Do not run the full Neoaves inputs locally. The production path uses one
ASTRAL-IV v1.25.4.8 inference per treatment with `--length SULength`,
`--sulengthtype 1`, and configurable `--genelength` (default 1000). The
annotated output contains both internal `CULength` and `SULength` values; the
runner extracts both after inference. Terminal CU values remain undefined.

`01_full363_primary.sbatch`, `02_full363_sensitivities.sbatch`, and
`03_j48_primary.sbatch` are retained as legacy sequential workflows. Use
`01_single_treatment.sbatch` for new runs. It executes exactly one treatment
and is restart-safe.

The expected Bridges-2 layout is:

```text
$PROJECT/neoaves_stage4e/
  repo/
  outputs/
  software/astral4
```

After transferring inputs and checking out the repository, set:

```bash
export MSRC_REPO=$PROJECT/neoaves_stage4e/repo
export NEOAVES_STAGE4E_DATA=$MSRC_REPO/empirical/neoaves_chr4/data/stage4e
export NEOAVES_STAGE4E_OUT=$PROJECT/neoaves_stage4e/outputs
export ASTER_BIN=$PROJECT/neoaves_stage4e/software/astral4
export PYTHON=python3
```

Run preflight first:

```bash
sbatch -A bio240015p empirical/neoaves_chr4/cluster/stage4e_aster/00_preflight.sbatch
```

Then submit the five primary treatments independently:

```bash
for T in T0 T_PNAS T_STRUCT_PRIMARY T_BLOCK_500K NONCHR4; do
  sbatch -A bio240015p --qos=rm -p RM-shared --cpus-per-task=32 \
    --time=16:00:00 --export=ALL,TREATMENT="$T" \
    empirical/neoaves_chr4/cluster/stage4e_aster/01_single_treatment.sbatch
done
```

Generate commands without submitting:

```bash
python empirical/neoaves_chr4/scripts/05_stage4e_make_commands.py --group primary
```

Collect returned outputs with `05_stage4e_collect_aster.py`. Existing legacy
SU-annotated outputs can be migrated without rerunning ASTRAL:

```bash
python empirical/neoaves_chr4/scripts/05_stage4e_run_aster.py \
  --extract-existing --input-annotated /path/to/T0/tree.su.nwk \
  --output-dir empirical/neoaves_chr4/results/stage4e/cluster/T0 \
  --treatment T0
```

Migration preserves the source checksum, refuses to overwrite a different
annotated tree, and writes `tree.annotated.nwk`, `tree.cu.nwk`, `tree.su.nwk`,
and `branch_lengths.tsv`. Derived files are written atomically.

Fixed-topology scoring is useful for branch-length comparisons on a common
topology, but it does not test topology selection:

```bash
astral4 -i MATCHED_INPUT.tre -C -c T0_TREE.nwk -u 2 \
  --length SULength --genelength 1000 --sulengthtype 1 \
  -o fixed.annotated.nwk
```

Use `sacct -j JOBID --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode`
for resource accounting. A preliminary Bridges-2 observation was about 9.5
hours for one full 363-species inference with 32 cores; this is a planning
measurement, not a runtime guarantee. Do not submit the 500 random-matched
inferences until deterministic outputs have been validated. Fixed-topology
and independent-topology comparisons must remain separate.
