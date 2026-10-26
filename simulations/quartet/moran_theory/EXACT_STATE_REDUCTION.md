# Exact unconditional quartet computation and state reduction

## Conclusion

For arbitrary finite (M), switching rate, selection, branch durations,
population sizes, and species-tree inheritance kernels, the exact
unconditional quartet calculation does **not** close in an (M)-independent
matrix whose state contains only the three topology probabilities. The hidden
Moran count remains part of the exact value function.

The exact operator for (n) unresolved lineages is

\[
A_n=Q_M\otimes I_{d_n}+
\operatorname{blockdiag}(G_n(0)^T,\ldots,G_n(M)^T),
\]

with (d_1=2,d_2=6,d_3=11,d_4=19). The state dimension is therefore
((M+1)d_n). Unconditional marginalization reduces the number of terminal
outputs from 16 pattern-specific columns to 3 topology columns, but it does
not remove the count coordinate.

## Why the count does not close

At count (k), the genealogy generator contains switching rates

\[
\rho k/M,\qquad \rho(1-k/M),
\]

and same-background coalescence rates

\[
1/\max(k,1),\qquad 1/\max(M-k,1).
\]

Consequently, two count states with the same value of any finite collection of
ordinary moments generally have different instantaneous effects on topology
observables. Applying the generator repeatedly produces inverse-count and
count-history functions, so a finite moment family is not invariant in the
general model. Selection makes the count generator itself depend rationally on
(k), and speciation introduces products of daughter value functions evaluated
through a shared parental count.

The count partition is therefore not strongly lumpable for the topology
observables: rows of the joint generator differ in their transitions and
absorption rates as (k) varies. A lumping that identifies all (k) would
require equal outgoing aggregate rates to every lump, which fails as soon as
the background-dependent coalescence rates are distinct.

This is an obstruction to a general (M)-independent exact reduction. It is
not a claim that no individual parameter setting has a symmetry reduction.

## Hierarchical species-tree dependence

In a balanced tree, left and right daughter contributions are multiplied only
after conditioning on their shared parental count. In a pectinate tree, the
sister pair and taxon 3 are combined through a common intermediate count,
then taxon 4 enters at the root count. This dependence prevents replacing the
count-indexed functions by products of marginal means. Distinct count
distributions with the same mean can therefore produce different
unconditional species-tree laws.

For example, with (M=10), (ho=0.17), root length (0.5), sister lengths
((0.3,0.1)), terminal lengths ((0.2,0.4,0.1,0.3)), and neutral dynamics,
the root initial distributions (\delta_5) and
((\delta_2+\delta_8)/2) have the same mean count but give, in the balanced
recursion,

\[
(0.36468612,0.31765694,0.31765694)
\]

and

\[
(0.36434068,0.31782966,0.31782966),
\]

respectively. The difference is small in this example but is far above
floating-point error and demonstrates that the mean frequency is not a
sufficient state for the hierarchical calculation.

The single-population (M=10) example in which both initial distributions
give ((1/3,1/3,1/3)) is a symmetry result: with no species-tree structure,
the four sampled lineages are exchangeable after marginalization. Its
conditional pattern laws still differ, and it does not establish a general
species-tree reduction.

## Exact special cases

Several narrower reductions are exact:

- **Fixed frequency or fixed count:** if (k) is externally fixed on a
  branch, the count coordinate disappears and the operator is (G_n(k)) of
  dimension (d_n). If different branches have different fixed counts, each
  branch uses its own small generator and inheritance is applied at the stated
  fixed frequencies.
- **Zero switching:** background configurations do not move, but the count
  coordinate still matters through the inverse-frequency coalescence rates
  unless the count is also fixed or a symmetry makes the requested observable
  insensitive to it.
- **Neutral symmetric settings:** complement or taxon symmetries can force
  equal discordant probabilities, and fully exchangeable single-population
  cases can force a uniform quartet vector. These reduce the answer by
  symmetry, not by a general closed generator.
- **Output marginalization:** the unconditional calculation can avoid storing
  16 conditional output columns. Existing balanced and unbalanced recursions
  already exploit this, while retaining the count-indexed intermediate
  functions needed for shared inheritance.

## Matrix-free exact implementation

`matrix_free_quartet.py` implements the same finite-state exponential through
`scipy.sparse.linalg.LinearOperator`. It applies (Q_M\otimes I_d) and each
small (G_n(k)^T) block directly, without assembling the full sparse block
matrix. `matrix_free_shared_population_quartet_probabilities` returns the
same unconditional shared-population result as
`integrated_quartet_probabilities`.

This is an implementation optimization, not a state reduction: the operator
still has dimension ((M+1)d_n), and the full hierarchical balanced and
unbalanced recursions still retain count-indexed arrays at inheritance events.
It is exact up to the same matrix-exponential floating-point tolerance as the
sparse implementation.

The benchmark uses (n=4), (ho=0.17), (T=0.5), and a fixed midpoint
initial count:

| (M) | dimension | assembled nnz | matrix-free stored nnz | sparse seconds | matrix-free seconds | sum error |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 209 | 1,581 | 1,239 | 0.0584 | 0.4808 | (le 2.3\times10^{-16}) |
| 20 | 399 | 3,151 | 2,429 | 0.1391 | 0.4410 | (le 2.3\times10^{-16}) |
| 40 | 779 | 6,291 | 4,809 | 0.1228 | 0.8217 | (le 2.3\times10^{-16}) |

The small benchmark favors sparse assembly because the matrices are tiny. The
matrix-free path is intended for avoiding block-matrix materialization at
larger dimensions or for integrating operator application into a larger
workflow; it should not be assumed faster for every (M).

## Practical answer

1. No general exact (M)-independent unconditional matrix was found or is
   supported by the generator algebra.
2. Exact reductions require fixed counts/frequencies or special symmetries;
   the formula is then a (d_n)-state genealogy propagation with the stated
   branch kernels.
3. The preventing dependence is the count-specific switching/coalescence
   generator and the shared-count products at speciation.
4. Unconditional calculation is smaller than the full conditional output
   because pattern columns can be summed, but it does not eliminate hidden
   count states.
5. The most scalable exact approach is sparse or matrix-free exponential
   propagation with count-indexed vectors, reusing (G_n(k)) blocks and
   exploiting branch-specific dimensions. Approximations such as diffusion or
   trajectory Monte Carlo should be labeled as approximations.

The new tests compare the matrix-free and sparse shared-population laws for
neutral and selected dynamics, all (n=1,2,3,4) generator dimensions, and
distinct same-mean count distributions in the balanced species-tree model.
