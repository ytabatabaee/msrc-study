# Methods

For each frozen local ML Newick tree, terminal and internal branch lengths were read directly from the source branch annotations. Patristic pairwise distances were calculated from the same trees. The primary all-sample `B(w)` is the mean over all non-*M. spretus* pairs. The leave-one-out version removes `Mmd_FRA4_18B.variant` from both numerator and denominator. Median and 5% trimmed means are separate robust summaries. Arrangement-informed `A_distance` and `A_relative` were recomputed with the existing Stage 7 within-subspecies t-versus-standard design. ALL_TIPS quartet support was recalculated with the frozen quartet definitions after removing the same tip from its group.

The branch detector uses the frozen Stage 13c shared configuration: robust MAD scaling, minimum segment size three, and the existing zero-to-four interval search. Reference inversion coordinates are used only for post-fit evaluation.
