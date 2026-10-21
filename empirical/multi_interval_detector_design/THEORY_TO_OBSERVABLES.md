# MSRC theory-to-observable audit

## Scope

This audit separates what the repository has derived from what it has only
evaluated numerically or observed in empirical tracks. It does not change the
frozen empirical results or implement a detector.

The exact finite-​M calculation is in
`simulations/quartet/moran_theory/moran_quartet_operator.py` and
`full_species_tree.py`. It averages a forward Moran count process with a
backward structured genealogy. For the shared population,

```text
J[z,r] = (pi_old kron h_r) exp(A4 T) b_z,
q_r = sum_z J[z,r],
P(Q_r | Z=z) = J[z,r] / P(Z=z).
```

The full `((1,2),(3,4))` recursion keeps the common parental count while
combining the two sister clades, uses binomial daughter inheritance, and
returns `P(Q,Z)`, `P(Z)`, `P(Q)`, and `P(Q|Z)`. Its tests and simulator
comparisons validate the stated hierarchical model. They do not prove that the
operator is the microscopic ancestral dual of one Moran reproduction process.
The effective switching intensity `rho` is therefore a model parameter, not
automatically a physical crossover rate.

## Evidence levels

| Level | Established result | Evidence |
|---|---|---|
| Analytically derived | Finite-​M Moran generator, `G1`, `G2`, `G4`, transpose/dual operator, terminal sampling, inheritance kernel, residual ordinary ancestor, and full joint quartet law | `moran_quartet_operator.py`, `full_species_tree.py`, `DERIVATION.md`, `FULL_SPECIES_TREE_DERIVATION.md` |
| Analytically derived special case | Fixed `p=1/2` `A_m`/`H_m` formulas and the symmetric constant-rate species-tree reduction | `simulations/quartet/full_species_tree_theory/` |
| Numerically evaluated theory | Sparse matrix values over finite `M`, branch lengths, frequencies, and `rho`; sensitivity of marginal versus conditional quartet distributions | `moran_theory/results/` |
| Simulation-supported | Agreement with the corrected `msrc-sim` Moran/full-tree simulator on the recorded validation grids | `moran_theory/results/POSTFIX_VALIDATION.md`, validation TSVs |
| Empirical/descriptive | `D(w)`, patristic-distance contrasts, SNAPP population-tree summaries, and local-tree summaries | `empirical/*/scripts/` and frozen result directories |
| Unresolved | A general analytic distribution for empirical branch lengths or realized coalescence times under the full finite-​M species-tree model; an observation model from genealogies to estimated trees; a proof of a microscopic Moran dual | `OPEN_THEORY_QUESTIONS.md` |

## Signals directly predicted by the current theory

The direct observable from the quartet theory is the joint law `P(Q,Z)`.
Topology shifts arise when the four lineage backgrounds have different
coalescence opportunities during the finite structured interval. Their size
depends on `M`, the frequency path, `rho`, structured branch durations, and
the terminal pattern. Conditional topology distributions can therefore change
substantially while the marginal `P(Q)` changes weakly; this is exactly why the
joint and conditional outputs are more informative than `q` alone.

The current theorem is topology-focused. The code records no general closed
form for an empirical branch-length distribution. A structured genealogy can
change waiting times even when integrated topology probabilities are close, but
that statement is a hypothesis for the current full observable pipeline until
it is tested with event-time output and a tree-estimation observation model.
The fire-ant and mouse/cod results show topology/length separation in empirical
tree summaries; they do not establish that MSRC caused it.

## Mutation-rate invariance

If every branch length in a local tree is multiplied by one common positive
factor `a_w`, every pairwise distance and every linear branch-length contrast
is multiplied by `a_w`. A ratio of two quantities from the same tree,
`A_distance / B`, is invariant to that common factor. A within-tree contrast
whose numerator and denominator have the same degree of homogeneity is also
invariant. Raw differences, means, total length, height, and mean internal
branch length are not invariant.

This cancellation does not hold for lineage-specific, branch-specific, or
class-specific rate factors. It also does not correct variation in tree
topology, missingness, SNP density, or estimator regularization. The audit
therefore treats `A_relative` and analogous ratios as more robust candidates
under a common window-wide rate multiplier, while retaining raw `A_distance`
and `B(w)` as descriptive signals. No feature is called mutation-rate robust
without stating this restricted invariance.

The exact feature definitions and units are in
[`FEATURE_AUDIT.tsv`](FEATURE_AUDIT.tsv).
