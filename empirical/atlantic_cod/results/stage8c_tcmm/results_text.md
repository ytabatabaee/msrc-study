# Stage 8c results

Published TCMM executed successfully for all four requested runs: outside and
all-window references at lambda 0 and 0.001. All 12 labels and the fixed rooted
baseline topology were preserved in every output tree.

The raw outputs are not ultrametric. At lambda 0, all-window outputs have a
maximum root-to-tip deviation of about 2.60, driven by extreme LG07 reference
branch lengths. The outside outputs have a maximum deviation of about 0.555.
Many outputs contain tiny negative edges of order 1e-22 to 1e-21, which are
solver numerical residues. They are retained in the raw outputs and counted in
the audit; no biological truncation was applied.

Consequently, the TCMM outputs are metric projections rather than validated
calibrated divergence-time trees. Stage 8b remains the valid calibrated
time-tree comparison. TCMM lambda 0 outside and all-window mean pairwise
MRCA-metric RMSE values are approximately 0.0358 and 0.399, respectively; the
lambda 0.001 values are approximately 0.0334 and 0.268. The large all-window
fit degradation reflects the retained extreme reference observations and the
metric objective, not a parameter adjustment to force agreement.

The direction of all-minus-outside shifts is therefore not directly comparable
as a calibrated-time result. TCMM shows substantial metric sensitivity in the
all-window treatment, while its non-ultrametric output prevents interpreting
that sensitivity as a robust divergence-time shift. The conclusion from Stage
8b is not replaced by Stage 8c.
