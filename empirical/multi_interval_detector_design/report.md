# MSRC theory-to-observable audit and detector design

## Conclusions

1. The directly predicted signal is the finite-M joint quartet law
   `P(Q,Z)`, including conditional topology changes caused by background
   sorting. The exact operator is documented in
   `simulations/quartet/moran_theory/` and has passed its recorded unit and
   simulator validations.
2. MSRC may change waiting times without a large integrated topology change,
   but the current repository has not derived or validated that general
   branch-length law. This is a simulation question, not a theorem in the
   present code.
3. Within-tree ratios such as `A_relative=A_distance/B` are invariant to a
   common multiplicative branch-length factor. Raw distances and means are not;
   lineage-specific and branch-specific rate variation remains confounding.
4. The joint detector can underperform because its modalities have different
   estimands, missingness, directions, noise, and spatial locations. Penalties,
   feature scaling, and the single-interval assumption amplify this mismatch.
5. The recommended detector is a penalized multivariate semi-Markov
   change-point model solved first by exact dynamic programming. It allows an
   unknown number of intervals, segment-specific unconstrained shifts, missing
   windows, unequal genomic spacing, and topology-only/continuous-only/joint
   modes.
6. Before implementation, essential experiments are MSC nulls, finite-M MSRC
   topology and joint regimes, a verified time-focused search, multiple
   separated intervals, frequency regimes, common and lineage-specific rate
   controls, and gene-tree estimation error.
7. The smallest defensible implementation for the paper is one frozen,
   multi-interval dynamic-programming detector with robust support-aware
   standardization, plus the predeclared simulation suite. It should be
   evaluated against frozen Phlag baselines and existing matched detectors;
   it should not be tuned to the structural annotations.

## What was inspected

The audit covers the finite-M quartet derivations and tests, the corrected
Moran/Wright--Fisher structured-to-panmictic boundary validation, the old
fixed-frequency `H_m/A_m` theorem, and the empirical feature/detector paths for
Atlantic cod, Neoaves, house mouse, Anopheles, and fire ants. The exact feature
definitions are tabulated in [`FEATURE_AUDIT.tsv`](FEATURE_AUDIT.tsv).

No frozen empirical results or detector implementation was modified for this
design task. The deliverables are methodological documentation and an
implementation plan. The next research priority is the controlled simulator
experiment that determines whether continuous genealogical observables are
actually MSRC-predicted before they are used as a main detector modality.
