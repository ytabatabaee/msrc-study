# Methods

Primary treatments retain the frozen Stage 4E input manifests. The production
runner performs one `--length SULength` ASTRAL-IV inference and extracts the
internal CU and SU annotations from that same annotated tree. It writes the
annotated source plus derived `tree.cu.nwk`, `tree.su.nwk`, and
`branch_lengths.tsv`, records the exact command and binary checksum, and
requires SLURM for full-scale execution. Internal branches are matched by
canonical unrooted bipartition; terminal branches by exact taxon label.
Terminal CU values remain missing because ASTRAL does not estimate them.
Relative changes are omitted when the baseline is zero or near zero. The
historical two-inference mode remains available only for legacy reproduction.
