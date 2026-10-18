# Stage 10b methods

The fixed scaffold is `stage8b_baseline_original_times.nwk`, verified against `stage8b_manifest.json`. It is rooted, fully bifurcating, time-calibrated, and used only to calculate the internal baseline-tree interval t for each four-taxon quartet.

Quartet eligibility required four unique taxa, a 2:2 arrangement state pattern, an arrangement split different from the baseline split, a resolved local topology in {arrangement, baseline, third}, and t>0. No Stage-7 MRCA values, Stage-10 parameter estimates, re-rooting, or projected Stage-8 time trees were used.

For supported configuration 1010, the validated finite-state calculation gives A=lambda/(lambda+2m)[1-exp(-2(lambda+2m)t)] and q=((1+2A)/3,(1-A)/3,(1-A)/3). The maximum difference between this implementation and the frozen exact simulation grid was 4.44e-16. Fitting uses quartet composite/pseudo-likelihood; it is not a full iid-gene-tree likelihood. CV uses five deterministic contiguous folds within each linkage group. The training-frequency null is fitted from training windows only, with Jeffreys 0.5-count smoothing for finite held-out log scores.

Profiles fix one parameter and re-optimize all nuisance parameters. Simulation validation uses the empirical eligible quartet intervals and window structure under the same conditional model. m is an effective reduced-model switching rate in inverse Ma units, not a meiotic recombination rate.
