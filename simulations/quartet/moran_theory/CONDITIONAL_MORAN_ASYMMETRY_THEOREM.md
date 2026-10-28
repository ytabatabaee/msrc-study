# Conditional Moran asymmetry audit

## Question

For `Z=0101`, write the conditional quartet law as

\[
q=((1-d)r,d,(1-d)(1-r)),\qquad d=q_2,
\quad r=\frac{q_1}{q_1+q_3}.
\]

The exact finite-`M` solver was used throughout. The calculation retains the
Moran count process, binomial daughter inheritance, and joint pattern mass.
It does not use the two-regime model or external mixtures.

## Findings

Root duration is a useful control of `d` in the paired-pattern construction:
longer root branches increase the opportunity for the dominant arrangement-
concordant resolution. This is a parameter sensitivity result, not a global
monotonicity theorem for all branch configurations.

The ratio `r` is exactly `1/2` in the audited symmetric zero-lower-branch
case. In the pectinate tree, asymmetric lower ancestral durations can move
the exact ratio to both sides of `1/2`. The committed certificate contains
examples with `r<1/2` and `r>1/2`, each with positive `P(Z=0101)`.

These examples establish that the model does not have a universal inequality
`q1 >= q3` for this pattern. They do **not** establish that every value of
`r` is attainable at fixed `d`, or that the required interval is connected.

## What remains open

No factorization or inverse of the form `d=f(u,v)`, `r=g(u,v)` has been
proved. In particular, the calculations do not show that an asymmetry-control
parameter can vary `r` while holding `d` fixed. The following implications
therefore remain open:

* coverage of the entire `Q2` dominance sector by `Z=0101`;
* coverage of the other two sectors by `Z=0011` and `Z=0110`;
* attainability of all six strict quartet orderings;
* full conditional simplex coverage or density;
* a conditional introgression-indistinguishability theorem.

The exact missing lemma is a global level-set result: for every feasible
`d`, the set of attainable `r` values must contain the full interval required
by the sector. A nonzero two-parameter Jacobian, opposite-side witnesses, or
vertex closure does not prove this.

## Reproducibility

Run:

```bash
python -m simulations.quartet.moran_theory.analyze_conditional_moran_asymmetry
python -m pytest simulations/quartet/moran_theory/tests/test_conditional_moran_asymmetry.py -q
```

The output table is
`results/conditional_moran_asymmetry/asymmetry_certificates.tsv`.
