# Atlantic cod Stage 12 report

Stage 12 adds a boundary-aware interval detector with an explicit no-change option and left/right censoring. It compares topology-only, time-only, and joint segmentation with the unchanged Stage 11 Phlag and sign-rule baselines.

The selected primary intervals are:
Boundary topology: LG01=internal:9000001-26000000; LG02=right_censored:18250001-24000000; LG07=left_censored:1-12000000; LG12=left_censored:1-13250000
Boundary time A: LG01=internal:9250001-26000000; LG02=right_censored:18500001-24000000; LG07=internal:13750001-23000000; LG12=internal:1500001-13250000
Boundary joint: LG01=internal:9000001-26000000; LG02=right_censored:18500001-24000000; LG07=internal:13750001-23000000; LG12=right_censored:13500001-27250000
Boundary topology mean full-reference IoU=0.6929; mean F1=0.7381.
Boundary time A mean full-reference IoU=0.9314; mean F1=0.9682.
Boundary joint mean full-reference IoU=0.7250; mean F1=0.7559.

LG02 is right-censored under the topology, time, and joint primary fits, with the selected transition at 18.25/18.50 Mb and the observed endpoint at 24 Mb. LG07 is recovered accurately by time A but not by topology; the topology forced-interval sensitivity selects a short early outlier rather than the known long inversion interval.

Phlag diagnostics are in `phlag_diagnostics.md`. All saved Phlag posterior probabilities and state-divergence values are `nan`, final focal-edge estimates collapse to `1e-06`, and the all-zero state path is therefore not a calibrated null-state result.

The Stage 12 time signal remains a contrast of SNAPP-derived population-tree estimates. It is not treated as direct realized gene-lineage coalescence time.
