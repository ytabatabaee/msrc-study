# Quartet-simplex coverage: what is proved

The image of a quartet-frequency map depends on the parameter class. This
audit distinguishes the original finite-`M` Moran model from the 38-state
two-regime model, and finite parameters from limits.

## A constructive two-regime family

Allow arbitrary row-stochastic inheritance kernels, set `M=1`,
`p_low=0`, `p_high=1`, `rho=0`, set all lower branches to zero duration, and
use deterministic taxon-specific kernels to create one of the patterns
`0011`, `0101`, or `0110`. For root duration `T`, the exact balanced recursion
gives

\[
q^{(j)}(T)=\left(\frac{e^{-2T}}3,\frac{e^{-2T}}3,\frac{e^{-2T}}3\right)
 +\left(1-e^{-2T}\right)e_j,
\]

where `j=1,2,3` for the three paired patterns. The explanation is direct: the
two matching background pairs each coalesce at rate one; if either coalesces,
the corresponding quartet topology is fixed, while the no-coalescence mass
contributes the ordinary unresolved quartet distribution.

Therefore each simplex vertex is approached exponentially as `T -> infinity`.
At every finite `T`, the other two probabilities remain positive. This proves
that the unrestricted two-regime model has all three vertices in its closure,
but it does not prove that its finite-parameter image is the entire simplex or
that its closure contains every interior point.

More generally, under finite branch durations and finite generator rates, any
positive unresolved four-lineage mass entering the root contributes positive
mass to all three absorbing topology states because the residual vectors have
positive entries on every unresolved state. Thus a finite-parameter result
with positive root-entry unresolved mass cannot be a simplex vertex. This is a
rigorous obstruction to exact vertex coverage for the usual finite-rate model;
the vertex construction reaches the boundary only as a duration tends to
infinity.

The exact certificates are in `results/simplex_coverage/edge_certificates.tsv`.

## What is and is not established

**Proved:** the displayed one-parameter edge-to-vertex family; vertex density
for the unrestricted two-regime model; the formula's agreement with the exact
38-state calculation in the certificate tests.

**Numerically supported:** the previous unrestricted inheritance witness with
an unconditional anomaly; the targeted identity, symmetric, and founding
searches with no anomaly; local rank-two variation at the previously tested
interior point.

**Not proved:** full-simplex coverage, positive-area overlap with the network
MSC image, or density of the entire simplex for unrestricted inheritance.

The family changes the inheritance kernels taxon by taxon and uses boundary
frequencies and an infinite-duration limit. It is therefore not a biological
inheritance theorem for identity or founding inheritance.

## Biologically constrained models

For identity inheritance, the targeted grid found no unconditional anomaly,
but no general inequality `q1 >= max(q2,q3)` has been proved for the pectinate
model. In the balanced model, taxon-exchange symmetry gives `q2=q3` under
matched left/right parameters; species-tree dominance then follows in the
tested finite cases, not from the symmetry alone.

For symmetric and founding-count inheritance, the same distinction applies:
the finite search is not a surjectivity proof. No coverage theorem is claimed.

## Original finite-M Moran model

The displayed construction does not transfer to the original Moran model.
Its count process retains `M+1` states, and the inheritance events are
binomial count kernels rather than arbitrary two-state taxon-specific
assignments. Exact trajectory integration remains count-dependent. No full
simplex or full-closure theorem was established for the Moran model.

Because the request concerns all of `Delta_2`, exact full coverage is already
not established for either model. The strongest rigorous result is the
unrestricted two-regime vertex-closure construction above; the biologically
constrained and original Moran cases remain unresolved beyond the reported
obstructions and numerical certificates.
