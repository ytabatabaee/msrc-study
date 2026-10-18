# Atlantic cod Stage 10b report

Stage 10b fits the exact conditional symmetric 2:2 quartet model to quartet outcomes and evaluates predictions in held-out contiguous genomic blocks. It uses the checksum-verified rooted Stage-8b original time tree only to obtain baseline quartet intervals. Stage-7 SNAPP MRCA ages do not enter the likelihood.

## Eligibility and model

The analysis retained 6406 of 70266 local quartet rows in the eligibility audit, with 3284 eligible inside-supergene rows used for fitting and CV. Exclusions are tabulated in `stage10b_quartet_eligibility.tsv`. The objective is a quartet composite/pseudo-likelihood; windows are the natural observational unit and window-level D(w) metrics are reported separately.

The frozen exact-kernel validation maximum absolute discrepancy was 4.44e-16. Parameters are effective reduced-model rates in inverse Ma; rho=m/lambda is dimensionless. They are not meiotic recombination rates or full biological MSRC estimates.

## Held-out prediction

Pooled held-out MSRC MAE=0.079678, RMSE=0.178893, mean log predictive probability=-0.454239; MSC MAE=0.875328; training-frequency null MAE=0.0803213. Per-LG results are in `stage10b_model_comparison.tsv`.

## Identifiability and adequacy

True nuisance-reoptimized profiles are in `stage10b_parameter_profiles.tsv`; profile diagnostics and model-adequacy audits are in `stage10b_identifiability.tsv`. The reduced model has no physical genomic-position term. Any within-LG variation in predictions comes only from the eligible quartet set and its baseline interval t; where those are constant, the predicted distribution is mathematically constant.

## Interpretation

The defensible claim is limited to whether this reduced two-background conditional model provides predictive structure for held-out arrangement-associated quartet outcomes. It is not an independent second-modality validation and does not establish the full biological MSRC process.
