# Stage 8b methods

Stage 8b starts from the original rooted and time-calibrated Atlantic cod baseline Newick from Matschiner et al. The three outgroups (`Gadcha_spc`, `Gadmac_spc`, and `Gadoga_spc`) were pruned. This left two redundant unary root nodes above the induced 12-tip Atlantic-cod root; both were collapsed while preserving the within-cod branch lengths. No MRCA data were used to choose or optimize the root. The resulting original tree is `stage8b_baseline_original_times.nwk`.

Original node ages were calculated as the maximum descendant path length below each node. The pruned source tree was ultrametric within a maximum root-to-tip deviation of 1.2e-11; the root age of the induced Atlantic-cod tree was 0.0653792 in the published SNAPP time units represented by the Newick. The source tree is an independent reference and was not fitted to Stage-7 data.

Pairwise Stage-7 MRCA times were aggregated by calling the existing Stage-8 `aggregate_matrices()` function, retaining its median rule and boundary-window exclusion. Outside-only and all-window matrices were projected onto the fixed original rooted topology using the existing Stage-8 tree isotonic least-squares routine. Only node ages were fitted; topology, root, and original branch lengths were not used as projection constraints.

The three tree versions are therefore: the original study reference, the same rooted topology with outside-window aggregate ages, and the same rooted topology with all-window aggregate ages. Projected ages are deterministic summaries of local SNAPP MRCA times, not new posterior SNAPP divergence-time estimates.
