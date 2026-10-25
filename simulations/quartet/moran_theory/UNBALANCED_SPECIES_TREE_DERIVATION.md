# Exact finite-M Moran-integrated law for `(((1,2),3),4)`

This document derives the implementation in `unbalanced_species_tree.py` for
the hierarchical model used by `msrc-sim`. The result is exact for that
hierarchical model: forward finite-population Moran count histories are
averaged with backward structured genealogies conditional on those histories,
and speciation uses conditionally independent binomial daughter counts. It is
not a claim that the transpose operator is the microscopic ancestral dual of a
single Moran reproduction process.

## Assumptions and notation

Let `M=2Ne` be the chromosome-copy count on a branch. Counts run forward from
the older endpoint to the younger endpoint with row generator `Q_M`, including
the optional genic selection transform already used by the balanced theory.
For a count `k`, write `p=k/M` and retain the simulator's finite-floor rates

\[
c_1(k)=1/\max(k,1),\qquad c_0(k)=1/\max(M-k,1).
\]

A lineage in background 0 switches to 1 at `rho p`, and a lineage in
background 1 switches to 0 at `rho(1-p)`. All durations use the simulator's
generation-equivalent scale.

The branch chromosome counts are `M_R` for the root branch, `M_I` for the
branch ancestral to 1,2,3, `M_S` for the 1,2 sister branch, and
`M_1,...,M_4` for terminal branches. Their durations are respectively
`T_A`, `T_I`, `T_S`, and `T_1,...,T_4`. The older root count law is
`pi_old`.

## The three-lineage generator

Order the eight transient states as

\[
z=(z_1,z_2,z_3)\in\{0,1\}^3
\]

using the same lexicographic product order as the existing operator code. Add
absorbing states `Q1=12|34`, `Q2=13|24`, and `Q3=14|23`. In the transient
block, each lineage switches independently at the rates above. The only
absorptions are

\[
R_3[z,1]=1_{z_1=z_2}c_{z_1}(k),\quad
R_3[z,2]=1_{z_1=z_3}c_{z_1}(k),\quad
R_3[z,3]=1_{z_2=z_3}c_{z_2}(k).
\]

The diagonal is minus the sum of all switching and absorption rates. Thus

\[
G_3(k)=\begin{pmatrix}S_3(k)&R_3(k)\\0&0\end{pmatrix}
\]

has dimension 11 and every row sums to zero, including at `k=0` and `k=M`.
The boundary values use exactly the existing simulator convention; they are
not replaced by a continuous-frequency limit.

The representation is sufficient. If any pair among 1,2,3 coalesces before
lineage 4 enters, that pair has a more recent MRCA than either remaining
lineage. Subsequent attachment and coalescence of lineage 4 therefore produce
the corresponding unrooted split: 12 gives `12|34`, 13 gives `13|24`, and 23
gives `14|23`. No later event can change that split. Consequently the
coalesced lineage does not need a background state after absorption.

For an `M`-copy branch, the dual operator is

\[
A_3(M)=Q_M\otimes I_{11}+
\operatorname{blockdiag}\left(G_3(0)^T,\ldots,G_3(M)^T\right),
\]

with dimension `11(M+1)` by `11(M+1)`. It has the same transpose convention
as `A_2` and `A_4`: a column of values at the younger boundary is propagated
to the older boundary by `exp(A_3 T)`.

## Recursion below the root

First propagate the terminal lineages 1 and 2 through their terminal branches,
combine their conditional background values at a common sister count, and
propagate with `A_2(M_S)` through `T_S`. Applying the binomial inheritance
kernel `H_{I,S}` gives a six-state array

\[
L_{12}(k_I,s),\qquad s\in\{00,01,10,11,C0,C1\},
\]

conditional on the count at the younger end of the intermediate branch. The
same terminal dual lift for taxon 3 followed by `H_{I,3}` gives
`L_3(k_I,b)` for `b\in\{0,1\}`.

At that node construct an 11-state younger-boundary array `Y_3`:

\[
Y_3(k_I,Q_1)=L_{12}(k_I,C0)+L_{12}(k_I,C1),
\]

and, for the eight unresolved states,

\[
Y_3(k_I,(z_1,z_2,z_3))
 =L_{12}(k_I,(z_1,z_2))L_3(k_I,z_3).
\]

The first line is essential: it carries sister-pair coalescences through the
intermediate count process instead of adding them as an independent scalar.
Propagate `Y_3` with `A_3(M_I)` through `T_I` and then apply `H_{R,I}`. The
result is `W_3(k_R,s)`, conditional on the root-node count. Its absorbing
columns contain all topology probabilities fixed below the root.

Let

\[
\pi_R=\exp(Q_{M_R}^T T_A)\pi_{old}.
\]

The early-topology contributions are

\[
c_i=\sum_{k_R}\pi_R(k_R)W_3(k_R,Q_i).
\]

This includes both coalescence of 1 and 2 on the sister branch and the first
coalescence among the three lineages on the intermediate branch.

## Root entry and final law

Propagate taxon 4 through its terminal branch and apply `H_{R,4}`, giving
`L_4(k_R,b_4)`. Define the root-entry unresolved array `F_pect` on the 19
states of `G_4` by

\[
F_{pect}(k_R,(z_1,z_2,z_3,z_4))
 =W_3(k_R,(z_1,z_2,z_3))L_4(k_R,z_4)
\]

for the 16 transient states, and zero for all three absorbed states. Thus
`F_pect` has dimension `19(M_R+1)` and contains only mass for which no
topology was fixed below the root. Its probability mass is evaluated after
weighting by `pi_R`; the unweighted sum over count-indexed conditional values
is not itself a probability.

The existing four-lineage operator and residual vectors can be reused without
change:

\[
A_4(M_R)=Q_{M_R}\otimes I_{19}+
\operatorname{blockdiag}\left(G_4(0)^T,\ldots,G_4(M_R)^T\right).
\]

With `h_i` equal to one third on each transient state and one on absorbed state
`Q_i`, the root contribution is

\[
r_i=(\pi_{old}\otimes h_i)^T
      \exp(A_4(M_R)T_A)F_{pect}.
\]

Therefore the exact pectinate law is

\[
\boxed{q_i=c_i+r_i,\qquad i=1,2,3.}
\]

The implementation uses `scipy.sparse.linalg.expm_multiply`; it never forms a
dense exponential for the branch operators.

## Probability conservation

For every count, `G_3` and `G_4` have zero row sums. The terminal sampling
vectors sum to one conditional on their count, and each inheritance kernel is
row stochastic. At the intermediate node, the sister-coalesced mass plus the
eight unresolved products sums to one. The `G_3` propagation partitions this
mass into its three absorbing columns and its eight transient columns. Hence
the early mass plus the unresolved root-entry mass is one after averaging over
`pi_R`.

The four-lineage residual vectors satisfy `h_1+h_2+h_3=1` on every transient
and absorbed state relevant to `F_pect`. Therefore the root contributions sum
to the unresolved root-entry mass, and `sum_i q_i=1`. The tests verify this
identity numerically at boundary and interior counts.

## Ordinary MSC limit and balanced-tree relationship

When every chromosome is in one background, `rho=0`, and the count is fixed
at `K=M`, the structured coalescence rate for every pair is `1/M`. With the
1,2 sister branch duration written as `T_S=M t` and no later topology bias,

\[
q_1=1-\frac23e^{-t},\qquad q_2=q_3=\frac13e^{-t},
\]

which is the ordinary MSC quartet law. The intermediate and root stages only
resolve otherwise unresolved lineages symmetrically in this limit.

The balanced recursion is recovered by replacing the pectinate three-lineage
stage with two independent sister-clade pair stages and entering `A_4`
directly. `A_4` and its residual vectors are identical; only the construction
of the root-entry unresolved vector and the early topology contributions
changes.

## Scope and limitations

The theorem is exact for the finite-state hierarchical simulator model,
including its binomial speciation inheritance and finite-floor boundary rates.
It is not a general theorem for every biological Moran model, nor does it
replace direct simulator validation. The validation driver compares independent
frequency-history/genealogy simulations from `msrc-sim` against the formula;
larger parameter grids are prepared but are not run automatically.
