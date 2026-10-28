# Rigorous conditional Moran coverage audit

## Pectinate discrepancy

The two earlier pectinate values use different chromosome-copy counts. The
historical image routine calls `pectinate_q(..., M=4)` to keep its exploratory
grid inexpensive. At `(T_A,T_1)=(1,0.2)`, `rho=0.1`, `Z=0101`, it gives

```text
q = (0.1380731890, 0.7286695969, 0.1332572141)
det(Df) = -6.4949035709e-4.
```

The requested common `M=10` configuration gives

```text
q = (0.2311334906, 0.5398346527, 0.2290318568)
det(Df) = -2.6196000042e-4.
```

The balanced value uses `M=10` in both analyses. The discrepancy is therefore
a finite-population parameter mismatch, not a topology-label or derivative
error. Regression tests now recompute both values and compare each against the
corresponding exact solver.

## What is certified

The analytic derivative implementation is correct for the complete finite-M
recursion: terminal derivatives are propagated through the lower branches,
early topology masses, root survivor propagation, and conditional pattern
normalization. At the common `M=10` point, finite differences at three step
sizes agree with the analytic Jacobian to below `3e-10` for both rooted shapes.

The current repository does **not** contain a formal interval-arithmetic or
validated-linear-algebra enclosure of the matrix exponential. High-precision
or finite-difference agreement cannot supply that missing proof. Accordingly,
the inverse-function theorem is recorded as a conditional mathematical result:
if the displayed determinant is enclosed by a certified interval excluding
zero, then the conditional image contains a nonempty open subset of the
simplex. The numerical results strongly support that premise but do not prove
it under a formal validated-computation standard.

## Global geometry

The bounded geometry audit varies root duration, taxon-1 terminal duration,
sister/intermediate durations, and switching rate at `M=4` for all three
patterns and both rooted shapes. It identifies parameter sensitivity but does
not establish a global boundary, convexity, or full-simplex coverage. The
target table from the restricted two-duration fits contains no noncentral
match below `1e-6`; this is a failed bounded search, not an impossibility
result. The center is available from the earlier exact zero-root calculation.

The three paired patterns approach different vertices, but vertex closure does
not imply interior coverage, and the conditional Moran image is not known to
be closed under convex mixtures. The global coverage question therefore
remains unresolved.

The hybridization line in `hybridization_line.tsv` and the SVG diagnostic are
illustrative. Because the local determinant is not formally interval-certified,
the repository does not claim a rigorous nontrivial overlap or a positive
length indistinguishability interval with the hybridization model.

## Reproduction

```bash
MPLCONFIGDIR=/tmp/msrc-mpl python -m \
  simulations.quartet.moran_theory.validate_conditional_moran_jacobian
MPLCONFIGDIR=/tmp/msrc-mpl python -m \
  simulations.quartet.moran_theory.conditional_moran_global_geometry
```

These are bounded exact calculations, not cluster jobs or broad optimization.
