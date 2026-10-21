# Methods

Primary treatments retain the frozen Stage 4E input manifests. The runner writes separate `tree.cu.nwk` and `tree.su.nwk` outputs for `--mode both`, records the exact commands and binary checksum, and requires SLURM for full-scale execution. Internal branches are matched by normalized unrooted bipartition; terminal branches by exact taxon label. Relative changes are omitted when the baseline is zero or near zero.
