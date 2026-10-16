# Stage 10 results

The primary full-data composite fit used 16113 compatible non-boundary pair-window observations. The held-out MSRC prediction used 176 windows and achieved MAE 0.6862, RMSE 0.7636, and correlation 0.4082 between predicted and observed eligible-quartet contrasts. The ordinary MSC comparator achieved MAE 0.8964, RMSE 0.9969, and correlation -0.4106; mean window-level log predictive probabilities were -0.9703 for MSRC and -1.289 for MSC.

The fit has 43 reported parameter rows, including the primary and linkage-group-specific sensitivity models. Profile values and block-bootstrap intervals are in the TSV outputs. Simulation validation covered six known CTMC parameter/state settings; the detailed rows are in `stage10_simulation_validation.tsv`.

Quartet eligibility was restricted to symmetric 2:2 background patterns supported by the exact available formula. Exclusion counts are {"arrangement_equals_baseline": 13399}. The result is therefore a parameter-to-prediction test for that mathematically supported quartet subset, rather than a claim of a full 12-population MSRC quartet distribution.
