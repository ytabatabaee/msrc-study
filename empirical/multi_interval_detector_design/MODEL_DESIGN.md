# Proposed general multi-interval detector

## Recommended model

Use a penalized multivariate semi-Markov change-point model implemented first
as dynamic programming over ordered genomic observations. Let

```text
X_w = (Q_w, C_w),
```

with a feature mask and reliability weight for every observed component. The
background state has a fitted mean and robust covariance. Each anomalous
segment has its own mean vector; its direction is unconstrained, so separate
intervals may have different topology and continuous shifts. A segment is
allowed only if it contains a minimum number of informative windows and a
minimum physical span. Gaps remain gaps; a transition may be penalized by
genomic distance or by the amount of observed support rather than silently
interpolating missing windows.

For a candidate segmentation `S`, use

```text
cost(S) = sum_segments robust_loss(X_segment; mu_segment, Sigma_segment)
         + penalty(number_of_segments, n, d),
```

with parameters selected without structural annotations. A practical first
version uses blockwise median/MAD standardization, diagonal or shrinkage
covariance, a Huber or trimmed quadratic loss, minimum segment size, and a
predeclared BIC-like penalty. Topology-only, continuous-only, and joint modes
are feature-mask restrictions of the same code path.

This is preferable to a topology-specific one-interval scan because it
supports an unknown number of intervals and heterogeneous anomaly regimes
while retaining an auditable objective. Generic dynamic programming costs
`O(K_max n^2 d)` for a bounded interval count and `O(n L_max d)` when a
maximum segment length is used. A semi-Markov HMM adds transition and state
probabilities but requires emission-model choices that are currently less
supported by the MSRC theory. PELT can reduce average cost, but its pruning
assumptions should be checked under spatial dependence before relying on it.

## Parameter selection

Freeze the standardization, minimum segment size, maximum interval count, and
penalty family before reading structural annotations. Use ordinary MSC and
permuted-coordinate controls to calibrate false segmentation. Structural
regions enter only in a final accuracy table. Do not choose feature weights,
interval count, or a direction constraint by maximizing structural overlap.

Report the selected segmentation, support mask, per-feature contribution,
segment means, penalty, and sensitivity to the predeclared penalty grid.
Report pairwise agreement and common-support agreement separately from
structural accuracy.

## Why not implement an HMM first?

An HMM would require a transition model for genomic distance, emission laws
for mixed topology and continuous features, and a rule for multiple anomaly
states. Those choices are not supplied by the current MSRC theorem. The
semi-Markov/dynamic-programming model exposes the assumptions directly and
can later be given probabilistic emissions after the simulation plan resolves
the observation model.
