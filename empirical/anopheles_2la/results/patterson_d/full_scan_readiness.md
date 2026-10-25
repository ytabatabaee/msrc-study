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

From the repository root, after checking the authenticated Linux environment:

```bash
cd empirical/anopheles_2la/analysis/patterson_d
mkdir -p logs
sbatch --export=ALL,MSRC_REPO="$PWD/../../..",PYTHON_BIN=/path/to/ag3-python,AG3_URL=gs://vo_agam_release_master_us_central1 run_full_patterson_d.sbatch
```

Do not use this command until the pilot reports and the remote credentials,
CA bundle, sample metadata, and source checksums have been verified. It is the
only command that can start the full scan; no submission was made here.
