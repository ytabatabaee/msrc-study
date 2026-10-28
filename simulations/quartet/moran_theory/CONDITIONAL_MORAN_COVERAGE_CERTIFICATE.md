# Conditional Moran local-coverage certificate

This audit uses the original finite-`M` Moran MSRC recursion, including
stochastic count propagation and binomial daughter-count inheritance. It does
not use the two-regime model and does not change any production solver.

## Analytic Jacobian

For `Z=0101`, `M=10`, initial count `K=5`, `rho=0.1`, sister durations
`(0.2, 0.3)`, terminal durations `(T1, 0.1, 0.2, 0.3)`, and root duration
`TA`, the audited map is evaluated at `(TA, T1) = (1, 0.2)`.

The derivative implementation propagates the exact identity

`d exp(A t)v / dt = A exp(A t)v`

through terminal, sister/intermediate, and root stages. It differentiates the
early-resolution mass and the conditional normalization as well as the root
survivor contribution.

| rooted shape | q1 | q2 | q3 | P(Z=0101) | det(J) |
|---|---:|---:|---:|---:|---:|
| balanced | 0.239918614 | 0.523140999 | 0.236940387 | 0.040281545 | -1.297436435e-4 |
| pectinate | 0.231133491 | 0.539834653 | 0.229031857 | 0.041043163 | -2.619600004e-4 |

Centered finite differences at steps `2e-4`, `1e-4`, and `5e-5` agree with
the analytic derivatives to at most `2.6e-10` for the balanced case and
`2.5e-10` for the pectinate case. The determinant remains negative across
those steps. Complete values are in
`results/conditional_moran_coverage/jacobian_certificate.tsv`.

This is a reproducible numerical certificate, not a formal interval-arithmetic
certificate. Consequently the inverse-function conclusion is conditional on
certifying the displayed determinant bound: analytically, a nonzero determinant
would imply a nonempty open conditional image in the two-dimensional simplex.
The repository does not claim that global coverage follows from this local
calculation.

## Image exploration and target checks

`image_samples.tsv` contains a bounded exact grid for the three patterns
`0011`, `0101`, and `0110` in both rooted shapes. It is descriptive and does
not establish an image boundary. `target_witnesses.tsv` records bounded
optimization attempts in a restricted balanced `M=4` two-duration family; the
noncentral targets were not matched to numerical tolerance. The center is
readily attained at zero root duration in the earlier audit, but a successful
center point does not establish interior coverage.

Thus the supported conclusions are local: both rooted shapes have stable
rank-two numerical Jacobians at the stated point, and all three patterns have
positive probability at the audited finite points. The union of sampled
images, complete interior coverage, and any global introgression-overlap
theorem remain unresolved. No unconditional claim follows because these are
arrangement-conditional results.

The figure `conditional_moran_global_image.pdf`/`.png` is a compact visual
diagnostic of exact solver evaluations, not an interpolated or simulated image.

## Reproduction

```bash
MPLCONFIGDIR=/tmp/msrc-mpl python -m \
  simulations.quartet.moran_theory.certify_conditional_moran_jacobian
MPLCONFIGDIR=/tmp/msrc-mpl python -m \
  simulations.quartet.moran_theory.explore_conditional_moran_global_image
```

The second command is intentionally a small bounded exploration and should
not be expanded into a broad parameter sweep without a separate computational
plan.
