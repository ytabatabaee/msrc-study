# Two-regime MSRC model

## Model definition

This is a new stochastic frequency model, separate from the finite-(M) Moran
model. Each species-tree branch carries a regime (E_b(t)in\{L,H\}) with

\[
R_b=\begin{pmatrix}-a_b&a_b\\b_b&-b_b\end{pmatrix},
\qquad p(E=L)=p_L,\quad p(E=H)=p_H.
\]

Conditional on the regime, background switching and same-background
coalescence use the existing genealogical generator (G_n(p,M,\rho)). The
copy number (M) remains a rate-scaling parameter, but it is not a stochastic
count state. The joint operator has dimensions 4, 12, 22, and 38 for one
through four lineages:

\[
\widetilde A_n=R_b\otimes I_{d_n}+
\operatorname{blockdiag}(G_n(p_L)^T,G_n(p_H)^T).
\]

The forward regime generator is not transposed; this follows the existing
Moran dual convention, in which the forward frequency-process generator is
combined with transposed backward genealogy blocks.

## Speciation inheritance

At a species-tree split, (H(e_p,e_d)) is a row-stochastic parent-to-child
regime kernel. Conditional on the parental regime, daughter regime draws are
independent. Identity inheritance is the default. Nonidentity kernels are
supported explicitly and are applied before multiplying daughter value
functions, preserving the parental-regime correlation.

For example, if (u_d(e_d,b)) is a daughter value at its older boundary,
the value conditional on the parental regime is

\[
u_p(e_p,b)=\sum_{e_d}H(e_p,e_d)u_d(e_d,b).
\]

Products are formed only after this conditioning is retained.

## Exact recursions

The balanced implementation constructs regime-indexed left and right sister
survivor arrays, applies root-to-daughter inheritance, and forms

\[
F_{\mathrm{bal}}(e,z_Lz_R)=L(e,z_L)R(e,z_R).
\]

The early term is (q_{\mathrm{sis}}(e)=1-s_L(e)s_R(e)), averaged over the
regime distribution at the younger end of the root branch. The root term is
obtained by applying the 38-state operator. Thus

\[
q_i^B=\delta_{i1}q_{\mathrm{sis}}+
(\pi_E\otimes h_i)^T e^{\widetilde A_4T_A}F_{\mathrm{bal}}.
\]

For the pectinate tree, the corresponding intermediate operator has 22 states
and the root operator has 38 states. Sister coalescence is seeded into the
(Q_1) absorbing state; the other two absorption states are generated on the
intermediate branch. After root-to-intermediate and root-to-taxon-4
inheritance,

\[
q_i^U=c_i+(\pi_E\otimes h_i)^T e^{\widetilde A_4T_A}F_{\mathrm{pect}}.
\]

The code also returns (P(Z=z,Q_i)), pattern probabilities, and
(P(Q_i\mid Z=z)). Pattern-specific early mass includes the sampling mass of
lineages that have not yet entered the genealogy when the topology resolves.

## Deterministic validation

The tests verify dimensions, nonnegativity, conservation, joint-to-marginal
identities, nonidentity inheritance, and the (p_L=p_H) special case. When
(p_L=p_H), regime switching factors out and the result agrees with the
prescribed fixed-frequency recursion. With zero regime switching, the result
equals the mixture of the two fixed-regime calculations weighted by the older
regime distribution. The latter agreement is exact to floating-point error.

The compact validation gives, for (p_L=0.2,p_H=0.8,
\rho=0.17), initial regime weights ((0.6,0.4)), and active regime
switching, approximately

\[
P(Q)=(0.36342249,0.31836560,0.31821190).
\]

The unconditional discordant probabilities are close but unequal. Several
terminal patterns show conditional anomalies; the largest observed conditional
anomaly margin in this diagnostic is approximately 0.2096, with pattern mass
retained in the output. These are model-specific exact calculations, not
claims about the Moran model or introgression.

## Comparison and interpretation

The regime model can reproduce conditional discordant asymmetry and
conditional anomalous quartets because regime paths and terminal sampling are
correlated. It does not automatically reproduce the finite-(M) Moran law:
matching a mean frequency, variance, or switching rate does not match the
full frequency-path distribution or shared-ancestry inheritance.

The current implementation validates the shared regime model and its complete
balanced and pectinate recursions analytically. A full empirical comparison to
the Moran model requires a predeclared parameter matching rule; it should be
reported as a sensitivity comparison, not calibration or equivalence.

The model is mathematically coherent as a finite-regime structured coalescent.
It is suitable as a secondary analytical model or sensitivity model. The
original finite-(M) Moran process remains the primary model when finite-count
drift and its induced frequency histories are scientifically required.
