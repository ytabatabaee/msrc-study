# Exact finite-M Moran MSRC quartet distribution on `((1,2),(3,4))`

This document describes the exact recursion implemented in
`full_species_tree.py` for the hierarchical model in `msrc-sim`:

1. forward Moran count histories on every population branch;
2. binomial daughter-count initialization at speciation;
3. backward structured genealogy conditional on those histories.

It does not claim a microscopic ancestral dual for the same Moran reproduction
events. The effective switching intensity `rho` is a background-switching
parameter, not by itself a physical recombination rate.

## Branch state spaces and time direction

On a branch with `M=2Ne`, the forward count `K` runs from older to younger
time with row generator `Q_M`. A backward genealogy starts at the younger
boundary. For fixed `K=k`, the background state changes `0 -> 1` at `rho*k/M`
and `1 -> 0` at `rho*(1-k/M)`. Same-background coalescence uses

```text
c1(k)=1/max(k,1),   c0(k)=1/max(M-k,1).
```

The code uses the dual operator

```text
A_n = Q_M kron I + blockdiag(G_n(0)^T, ..., G_n(M)^T)
```

with row-vector interpretation. The transpose is required because the Moran
history is forward while the genealogy is backward.

The implemented conditional genealogy spaces are:

```text
n=1: 0, 1
n=2: 00, 01, 10, 11, C0, C1
n=4: 16 transient background states plus Q1,Q2,Q3
```

In `G2`, `00 -> C0` occurs at `c0(k)` and `11 -> C1` at `c1(k)`. `C0` and
`C1` continue switching as one ancestral lineage; coalescence does not freeze
the background.

## Daughter inheritance and the shared parental count

At a speciation event, for a parent with total `M_P` and count `k`, each
daughter with total `M_D` starts with

```text
H_{M_P -> M_D}(k,a)
  = choose(M_D,a) (k/M_P)^a (1-k/M_P)^(M_D-a).
```

The left and right daughter counts are conditionally independent given the
common parent count, but are correlated after that count is marginalized. The
recursion therefore keeps the root-node count `k_P` explicit.

For each sister clade, terminal one-lineage dual lifts are combined at a fixed
sister count, propagated with `A_2`, and then mapped to the root-node count by
the inheritance matrix. This gives two six-state vectors

```text
L[k_P, state], R[k_P, state].
```

Only after conditioning on `k_P` are they multiplied. The four-survivor root
entry mass is

```text
F[k_P,(z1,z2,z3,z4)]
  = L[k_P,(z1,z2)] R[k_P,(z3,z4)]
```

for the 16 combinations in which neither sister pair coalesced. If either
sister pair is already in `C0` or `C1`, the unrooted quartet topology is fixed
to `Q1=12|34`; this contributes

```text
q1_sister = sum_k pi_root[k]
             [1 - s_L(k) s_R(k)],
```

where `pi_root` is the root-node count distribution after the forward root
Moran branch and `s_L,s_R` are the four-state survivor masses conditional on
that root count.

## Full operator formula

Let `pi_old` be the root count distribution at the older end of the root
shared branch. Let `F` be the conditional four-lineage root-entry vector just
defined, with zero mass in the absorbed root topology states. Let `h_r` be
the existing residual vectors: one third on transient states and unit mass on
absorbed topology `Q_r`. Then the root contribution is

```text
q_root[r] = (pi_old kron h_r) exp(A_4 T_A) F.
```

The full quartet law is therefore

```text
q_r = q_root[r] + 1[r=1] q1_sister.
```

This is an exact finite-state operator expression. There are no arbitrary
terminal-pattern weights and no independence approximation between the two
daughters after marginalizing their parent count.

For the full joint law, repeat the construction with terminal sampling
vectors conditioned on each pattern `z`. If `F_z` is the corresponding root
entry vector, then

```text
J[z,r] = (pi_old kron h_r) exp(A_4 T_A) F_z
         + 1[r=1] J_sister[z].
```

Here `J_sister[z]` is the mass for which at least one sister pair coalesces
before entering the root, and therefore has topology `Q1`. The implementation
returns

```text
P(Z=z)       = sum_r J[z,r]
P(Q_r)      = sum_z J[z,r]
P(Q_r|Z=z)  = J[z,r] / P(Z=z).
```

This derives `P(Q,Z)` from the same terminal, inheritance, sister-branch, and
root operators; it does not reconstruct the conditional law from arbitrary
pattern weights.

## Fixed-frequency regression

With all counts held at `K=M/2`,

```text
m = rho/2,   lambda = 2/M = 1/Ne,
```

the fixed-count conditional interface reproduces the old constant-rate
species-tree result. Averaging that conditional result over the 16 independent
terminal patterns gives the unconditional fixed-`p=1/2` quartet law. For the
specific conditional pattern `1010` and zero sister branches, it gives

```text
( (1-A_m)/3, (1+2A_m)/3, (1-A_m)/3 ).
```

The `A_m` statement is conditional on the pattern. An unconditional
fixed-`p=1/2` average over all terminal patterns is symmetric and generally
uniform.

## Structured origin and the ordinary ancestral phase

The genealogy simulator now carries the explicit rearrangement origin age

```text
origin_age = origin_branch.older_age
             - origin_time_from_branch_start.
```

While backward time is younger than this boundary, the structured process uses
frequency-dependent switching and same-background coalescence. Once the
boundary is crossed, the rearrangement is not yet present in the older
population: switching stops, background labels cease to constrain events, and
all lineage pairs coalesce at the ordinary pairwise rate `1/(2 Ne)`. The
frequency history is therefore not used as a frozen oldest-frequency process
after the origin. For a root origin at the oldest end of the root branch this
age equals `tree.root.age + tree.root_extension`; that equality is a special
configuration, not the definition of the boundary.

The operator theory has the same two-phase interpretation. Its structured
operator covers the finite interval where the rearrangement exists, and the
residual unresolved mass is resolved by the ordinary panmictic ancestor. This
is exact for the hierarchical simulator model and is not a claim about a
single microscopic Moran ancestral dual.

## A finite-M boundary-condition distinction

The literal finite-M `msrc-sim` species tree performs binomial daughter
initialization at every speciation, even when the branch duration is zero.
Consequently, setting both sister branch lengths to zero does not produce the
same terminal sampling law as the separate shared-population kernel, which
samples all four taxa independently conditional on one common terminal count.
The daughter count is shared within each sister clade before the terminal
samples, producing finite-M correlation.

Thus there are two distinct limits:

1. the **literal simulator species-tree limit**, retained by this module,
   which includes zero-duration binomial inheritance events;
2. the **collapsed shared-population kernel**, obtained by removing the
   speciation inheritance steps and represented by `integrated_quartet_probabilities`.

The tests make this distinction explicit rather than forcing the two operators
to agree. The direct simulator validation agrees with the literal recursion on
positive sister branches after the structured-to-panmictic boundary fix. This
zero-duration distinction remains intentional: it is a difference between
removing the speciation events in the collapsed kernel and retaining the
literal simulator's instantaneous inheritance events.
