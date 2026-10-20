# Atlantic cod Stage 12b report

Stage 12b preserves Stage 12 and evaluates four prespecified boundary objectives. Model A is the exact Stage 12 joint detector. Model B uses robust median/MAD scaling with equal total topology/time weight. Model C requires positive inside-minus-outside shifts in both D(w) and A(w), based on the predeclared biological hypothesis that recombination suppression increases arrangement-associated topology and opposite-minus-same SNAPP population-tree divergence contrasts. Model D combines independent standardized topology and time profile improvements with equal weight.

The A(w) feature is a contrast of SNAPP-derived population-tree divergence-time estimates, not realized gene-lineage coalescence times.

## Primary model comparison

- Model A current: mean full-reference IoU=0.7250; mean F1=0.7559; mean observable boundary error=3276.4 kb.
- Model B robust: mean full-reference IoU=0.9458; mean F1=0.9780; mean observable boundary error=205.3 kb.
- Model C direction: mean full-reference IoU=0.7022; mean F1=0.7323; mean observable boundary error=414.0 kb.
- Model D profiles: mean full-reference IoU=0.7250; mean F1=0.7559; mean observable boundary error=3276.4 kb.

## LG12 diagnosis

The original joint detector is evaluated in `joint_model_diagnostics.tsv` and `candidate_interval_profiles.tsv`. These tables show each candidate's separate topology and time improvements, combined objective, penalty, segment means, and direction. The failure is assessed as a consequence of the unconstrained combined objective when the two modality shifts favor different intervals; Model C is a prespecified directional sensitivity rather than a coordinate-tuned repair.

- Model A current: right_censored 13500001–27250000; score=107.6; IoU=0.0314.
- Model B robust: left_censored 1–13750000; score=146; IoU=0.9147.
- Model C direction: left_censored 1–13250000; score=90.49; IoU=0.8967.
- Model D profiles: right_censored 13500001–27250000; score=40.59; IoU=0.0314.

## Arrangement-blind features

The arrangement-blind feature selected for Figure A was `blind_mean_internal_branch`, chosen by the largest pooled maximum interval score without reading inversion coordinates. For LG07 it selected `internal 13750001–23000000` with IoU=0.9566 and F1=0.9730. Full feature summaries and fitted intervals are in `arrangement_blind_features.tsv` and `arrangement_blind_predictions.tsv`. These summaries describe local SNAPP population-tree scale and branch-length variation; they are not direct coalescence-time observations.

## Detection significance and stability

The maximum-score p-values use a deterministic five-window block permutation that preserves within-block ordering but is heuristic for spatial dependence and the irregular window grid. Holm adjustment is applied across the four linkage groups within each method. Boundary perturbation ranges are descriptive and are not formal confidence intervals.

## Interpretation

The four linkage groups are too few for a strong generalization claim. The arrangement-informed time signal should be carried forward as a supplementary feature, with Model A retained as the frozen reference and Model B/D considered exploratory until tested on additional species. The arrangement-blind analysis is the relevant test for discovery of previously unknown regions.
