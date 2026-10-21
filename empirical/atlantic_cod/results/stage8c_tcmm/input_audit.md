# Stage 8c input audit

Stage 8c used the fixed rooted 12-population Atlantic-cod topology from
`data/atlantic_cod/processed/baseline_population_tree_12taxa.nwk`. The source
file retains a long unary external stem inherited from the outgroup-containing
tree. Following the already validated Stage 8b induced-tree construction, the
TCMM query copy removes only that external stem by setting the root edge to
zero. All 12 labels, the rooted topology, and every within-cod branch length
are retained; no root was selected using MRCA measurements.

The reference trees are the 426 published 250-kb MCC-like window trees in
`data/atlantic_cod/processed/cod_window_trees.tsv`, derived from Nature Source
Data Fig. 4. The outside treatment contains 250 rows whose frozen
`inside_inversion` field is false; the all treatment contains all 426 rows.
Each reference has the same 12 population labels and finite branch lengths.

Per-window posterior `.trees` samples were not available locally and were not
listed as per-window files in the frozen Stage 7 inventory. Therefore this is
direct TCMM on published window MCC-like reference trees, with equal weight per
retained window, not a posterior-tree bootstrap analysis.

Stage 7 pairwise MRCA summaries are used only for fit diagnostics. They are
median MRCA times, whereas TCMM's objective matches tree metrics/patristic
distances. These quantities are reported separately and are not treated as
interchangeable.

The all-window reference set includes extreme LG07 branch-length observations.
Those observations are retained rather than clipped or excluded; their effect
on TCMM non-ultrametricity and fit diagnostics is reported.
