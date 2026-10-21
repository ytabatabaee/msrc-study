# Atlantic cod Stage 8c report

Stage 8c successfully executed the published TCMM implementation at commit
`e6064b46ed0359ab7358d5110521fce1b1a825f7`. It used 250 outside-supergene and
426 all-window published 250-kb MCC-like reference trees, with lambda 0 and
0.001. The original 12-population rooted baseline topology and labels were
preserved. Stage 8 and Stage 8b outputs were not modified.

The source data contain branch-length-bearing published window trees, so direct
TCMM was possible. They do not contain per-window posterior tree samples; the
analysis is therefore based on published MCC-like reference trees rather than
posterior resampling. The first audit identified a long unary root stem in the
induced source Newick. The rerun removes only that stem in the TCMM query copy,
matching Stage 8b's validated induced baseline; the initial stem-retaining raw
run is preserved under `raw/legacy_root_stem_query/`.

Raw TCMM trees are finite and preserve the query topology, but they are not
ultrametric. Tiny negative edges occur at numerical solver scale. The all-window
set also contains extreme LG07 branch lengths, yielding large non-ultrametric
deviations and much poorer pairwise MRCA-metric diagnostics than the outside
set. These outputs cannot be interpreted as calibrated divergence-time trees
without an additional explicitly constrained method.

Thus TCMM provides a useful metric-matching sensitivity audit, but it does not
currently provide a better-supported time-tree estimator for the cod analysis.
The Stage 8b conclusion about supergene-associated temporal distortion remains
the calibrated result. Agreement or disagreement between TCMM and Stage 8b
cannot be reduced to a preferred numerical lambda because the two outputs
target different objects once TCMM loses ultrametricity.
