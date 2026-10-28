# Arrangement-conditional simplex coverage under the finite-M Moran model

This audit uses the original count-indexed Moran operators and binomial
daughter-count inheritance. It does not use frozen frequencies, the two-regime
model, or externally chosen inheritance kernels.

## Proved conditional vertex-density result

Consider the balanced tree, neutral Moran dynamics, `rho=0`, equal finite copy
counts `M` on all branches, a fixed interior initial count `K`, zero terminal
and sister durations, and a finite root branch of duration `T`. Condition on
one of the paired patterns `0011`, `0101`, or `0110`.

For `0011`, the only homogeneous pairs are `(1,2)` and `(3,4)`, so the
species-tree topology is selected when either pair coalesces. For `0101`, the
corresponding selected topology is `Q2`; for `0110`, it is `Q3`. The exact
finite-M generator has a finite transient unresolved block. Its survival mass
`S_z(T)` tends to zero as `T` tends to infinity, while the older-boundary
residual law is uniform over the three topologies. Consequently,

\[
P(Q_i\mid Z=z)=
(1-S_z(T))\mathbf 1\{i=j(z)\}+S_z(T)/3.
\]

This proves that the original finite-M Moran model approaches all three
simplex vertices conditionally on suitable observed arrangement patterns. The
arrangement probability remains positive for every finite `T` in the reported
interior-count cases, although it can become small.

The same construction transfers to the pectinate implementation when the
sister and intermediate branches have zero duration: all four lineages enter
the root stage, and the exact pectinate certificates approach the same three
vertices. This is a valid limiting subfamily of both rooted shapes, not a
claim that arbitrary positive lower-branch durations have the same formula.

The result is a density statement for the vertices, not a proof that the union
of conditional images contains every interior point. The exact `S_z(T)` is
count-process dependent; replacing it by `exp(-2T)` is only a frozen-frequency
special case and is not valid for the stochastic Moran process in general.

## What was not proved

The exact finite-M conditional map has not been shown to be surjective for one
pattern or for the union over patterns. A rank-two local Jacobian, observed for
`Z=0101` when root duration and one terminal duration vary, establishes only a
local two-dimensional image. It does not establish global coverage.

The target vectors `(0.2,0.5,0.3)`, `(0.05,0.9,0.05)`, and `(0.7,0.1,0.2)` are
therefore recorded as unresolved targets rather than claimed matches. No
introgression-indistinguishability corollary for arbitrary unconditional
vectors follows.

For the balanced symmetric subfamily, taxon-exchange symmetry forces equality
of the two non-selected probabilities for each paired pattern. Breaking that
symmetry with valid branch durations, switching, selection, or population-size
parameters can produce a two-dimensional local image, but no global inverse
has been derived. Pectinate trees may have fewer such symmetries, but that does
not by itself prove surjectivity.

## Audit of the earlier argument

The earlier vertex construction was conditional on `Z`, not unconditional. Its
fixed-frequency exponential formula described a frozen-frequency approximation
or a limiting special case. The present certificates use the exact
finite-M Moran-integrated joint law and divide by the exact positive pattern
probability. Thus the vertex result transfers, but the stronger claim of full
simplex coverage does not.
