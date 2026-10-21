# Simulation plan before implementation

All experiments should use independent forward histories as the Monte Carlo
unit. Loci nested within one frequency history require history-level bootstrap
or clustering. Parameters and seeds must be frozen before evaluating real data.

| experiment | generative setup | predicted signal | observables | failure mode |
|---|---|---|---|---|
| MSC null | ordinary MSC, no structured background | no MSRC topology or event-time shift | quartet support, true coalescence times, estimated trees | false spatial segmentation |
| topology shift | finite-M Moran/Wright-Fisher, nontrivial `rho`, background-associated pattern | conditional `P(Q|Z)` and `P(Q,Z)` shift | topology contrasts and joint pattern law | marginal `q` too weak |
| time-focused candidate | vary structured duration, `M`, frequency regime, and event-time summaries while recording true genealogy times | test whether waiting times move with small topology change | true coalescence times first; estimated branch lengths second | no branch-length-only regime exists under tested parameters |
| joint shift | parameters producing both background sorting and altered waiting times | topology and continuous blocks change together | `P(Q,Z)`, event times, tree lengths | incompatible direction across intervals |
| multiple intervals | concatenate independent windows with two or more parameter regimes and gaps | unknown number of separated segments | segmentation metrics and common-support agreement | single-interval outer hull |
| frequency regimes | low, balanced, and high initial frequencies; Moran path variation | conditional sensitivity depends on `p(t)` | joint topology/pattern law | weak signal at absorbing endpoints |
| mutation controls | multiply all local branch lengths by a window factor; then add lineage-specific factors | ratios invariant only to common factor | raw contrasts, ratios, normalized within-tree features | false robustness under lineage-specific rates |
| tree-estimation error | perturb true gene trees through topology error, unresolved quartets, and branch-length noise | observation signal attenuates and can become biased | true versus estimated features | detector learns estimator artifacts |

The current `msrc-sim` theory validation establishes topology agreement for the
hierarchical Moran model. It does not establish the time-focused row. That row
is a required experiment, not an assumption. If it fails across a declared
parameter grid, continuous features should be treated as generic phylogenomic
covariates rather than MSRC-predicted observables.

Each experiment should report topology RMSE, conditional-law error, event-time
distribution error, branch-length calibration error, interval precision/recall/
IoU/F1, false-positive rate under MSC, and recovery as a function of missing
windows and unequal spacing. Structural annotations are not used to fit or
select the detector.

The real-data comparison should then use frozen House mouse Stage 13d and
Neoaves Stage 11C baselines, on matched domains, with topology-only,
continuous-only, and joint multi-interval configurations evaluated using the
same support rules.
