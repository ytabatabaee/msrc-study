# Atlantic cod Stage 10 report

## Implemented analysis

Stage 10 freezes a time-calibrated version of the independent collinear population-tree topology using outside-supergene SNAPP MRCA times, fits a two-background MSRC pairwise composite likelihood, and uses the fitted coalescence-time parameters to predict held-out quartet contrast in contiguous genomic blocks.

## Estimates and uncertainty

The primary estimates and sensitivity estimates are in `stage10_msrc_parameter_estimates.tsv`. The likelihood profiles are in `stage10_msrc_parameter_profiles.tsv`; bootstrap parameter replicates use seed 20261016. Interpret m as an effective ancestry-background switching rate in inverse SNAPP time units.

## Identifiability and prediction

The reported profile shapes and parameter bootstrap intervals should be used to judge whether lambda and each m are separately estimable. Held-out prediction summaries are in `stage10_prediction_summary.tsv`, with window-level predictions in `stage10_cv_predictions.tsv`. The primary result is the predictive comparison, not a p-value.

## Limitations

11673/27786 pair-window observations satisfy T < d and were excluded from the primary fit. The clipped fit is an explicit sensitivity analysis. Quartet prediction is limited to the exact symmetric 2:2 `A_m(t)` case available in the repository theory. No unsupported full quartet distribution was imputed. A good fit would support the statement that coalescence-time-estimated MSRC parameters predict held-out arrangement-associated quartet bias; it would not establish MSRC uniquely.
