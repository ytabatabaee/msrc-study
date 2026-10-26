# Exact unconditional low-dimensional reduction audit

## Result

For the original hierarchical finite-(M) Moran MSRC model, no general exact
unconditional quartet formula with dimension bounded independently of (M) was
found. More strongly, a necessary intermediate observable already has a number
of independent time modes that grows with (M).

This conclusion concerns an exact finite-dimensional time-homogeneous linear
closure that works for arbitrary initial count distributions, branch lengths,
switching rates, and inheritance settings. It does not rule out special
symmetry reductions, fixed-frequency calculations, or approximations.

## Backward equation

For a column-valued conditional value (u_i(k,z,t)), the joint backward
operator used by the existing exact implementation is

\[
A_n=Q_M\otimes I_{d_n}+
\operatorname{blockdiag}(G_n(0)^T,\ldots,G_n(M)^T).
\]

The count coordinate cannot generally be integrated by summing over (k),
because (G_n(k)) contains

\[
\rho k/M,\quad \rho(1-k/M),\quad
1/\max(k,1),\quad 1/\max(M-k,1).
\]

Applying the operator to any proposed finite collection of ordinary moments
produces new inverse-count-weighted observables. Selection also changes the
count generator, and speciation multiplies daughter value functions
conditional on a shared parental count. Thus the products needed by the
hierarchical recursion are not determined by marginal moments.

## A growing-mode obstruction

Set switching to zero and consider two lineages in background 0 on one branch.
While they remain unresolved, the count process is governed by the killed
Moran generator

\[
B_M=Q_M-\operatorname{diag}\left(\frac{1}{\max(M-k,1)}\right)_{k=0}^M.
\]

The survival probability is

\[
s_0(t)=\pi^T\exp(B_Mt)\mathbf 1.
\]

The interior block on (k=1,\ldots,M-1) is an irreducible tridiagonal
matrix. For the neutral rates used by `msrc-sim`, it is diagonally
symmetrizable and has a simple spectrum generically. Consequently, the
interior survival contribution contains (M-1) independent exponential
modes for generic initial distributions. Any exact linear state realization
for this family of survival functions therefore needs dimension at least
(M-1).

This survival mass is not an artificial auxiliary quantity: it is the sister
pair survival term in the balanced recursion,

\[
q_{\mathrm{early},1}=1-s_Ls_R.
\]

By varying the sister branch duration and the initial count distribution, the
unconditional quartet output inherits this growing family of functions. The
same hidden count is retained through the pectinate intermediate branch and
through the shared-parent inheritance products.

The script `analyze_observable_reduction.py` computes the interior spectrum
and records the same-mean-count species-tree counterexample. It is a finite
diagnostic supporting the analytic obstruction; the obstruction itself follows
from the killed-generator representation and does not rely on Monte Carlo.

## What does reduce

The following cases are exact but narrower:

* A fixed count or prescribed frequency removes the count coordinate and uses
  the 2/6/11/19-state genealogy operators.
* Exchangeable single-population settings can force (q_2=q_3), and often the
  full quartet vector is (1/3,1/3,1/3). This is an output symmetry, not a
  general closure.
* Marginalizing terminal arrangement patterns reduces the number of outputs to
  three, but not the hidden count states needed at branch boundaries.
* The matrix-free implementation avoids assembling the block matrix while
  retaining exactly the ((M+1)d_n) mathematical dimension.

The existing prescribed-frequency result is therefore correctly interpreted as
\[
P(Q_i)=\mathbb E_{\mathcal P}\left[q_{\mathrm{early},i}(\mathcal P)+
F(\mathcal P)^TP_4(\mathcal P)h_i\right],
\]
not as one unconditional small matrix.

## Separate low-dimensional alternative

`regime_switching_quartet.py` defines a different model. A regime
\(E\in\{L,H\}\) switches with generator

\[
R=\begin{pmatrix}-a&a\\b&-b\end{pmatrix},
\]

and each regime has fixed frequency (p_L) or (p_H). Its shared-population
operator is

\[
\widetilde A_4=R\otimes I_{19}+
\operatorname{blockdiag}(G_4(p_L)^T,G_4(p_H)^T),
\]

with dimension 38. This is an exact unconditional formula for the explicitly
defined regime model, not for the original Moran process. A multi-branch
version would need a regime inheritance kernel at speciation and must be
specified separately.

## Conclusions for the manuscript

The original finite-(M) model should use the count-indexed sparse or
matrix-free operator for exact unconditional results. The prescribed-frequency
formula is appropriate for conditional analyses and for trajectory averaging.
The 38-state regime model may be used as a transparent sensitivity model, but
its results must not be described as an exact reduction of the Moran model.

No universal finite genealogical dual for the implemented frequency-dependent
coalescence rates was established. Constructing a different biological model
with a finite ancestral dual would require changing and revalidating the
simulator assumptions.
