# First complete symmetric four-taxon MSRC species-tree theorem

> Scope note: the constant-frequency `lambda` calculation in this document is
> a fixed-frequency `p(t)=1/2` special-case result. It is not the general
> Moran-consistent MSRC quartet law. `p=1/2` does not imply `T_L=T_R=0`;
> `T_L=T_R=0` is the shared-population limit, while `A_m` additionally fixes
> the shared-population frequency at `p(t)=1/2`. See
> `../moran_theory/DERIVATION.md` for the exact finite-M Moran-integrated
> shared-population kernel.

This derivation conditions on a fixed terminal background pattern
`Z=(z_1,z_2,z_3,z_4)` for the rooted species tree
`S=((1,2),(3,4))). It uses the symmetric backward process

\[
m_{01}=m_{10}=m,\qquad \lambda_0=\lambda_1=\lambda.
\]

The pattern probabilities `pi_z=P(Z=z)` are introduced only after the
conditional theorem. They are arbitrary forward-sorting weights.

The implementation in `quartet_theory.py` is independent of the old
conditional experiment but reproduces its `H_m` kernel. The simulator files
used as the model authority are `msrcsim/structured_coalescent.py`,
`msrcsim/analytic.py`, and `msrcsim/pattern_probabilities.py`.

## Theorem

Let `u_i` be the terminal single-lineage branch length for taxon `i`,
`T_L,T_R` the shared-population branch lengths for pairs `12` and `34`,
and `T_A` the duration of the final four-lineage structured ancestral
population. Define `ell_Z` and `r_Z` as the unnormalized distributions of
the surviving pairs at the ends of the two sister branches. Then

\[
\boxed{
\mathbf q_Z
=(1-s_Ls_R)\mathbf e_1+
(\ell_Z\otimes r_Z)H_m(T_A)
}
\]

where

\[
s_L=\ell_Z\mathbf 1,\qquad s_R=r_Z\mathbf 1,
\qquad \mathbf e_1=(1,0,0),
\]

and `H_m(T_A)` is the exact shared-population kernel with residual mass
resolved uniformly among the three topologies.

The first term is the probability that at least one sister pair coalesces on
its own species-tree branch. Either such event fixes `Q_1=12|34`, including
the case in which both sister pairs coalesce. The second term is the
probability-weighted shared ancestral calculation conditional on both pairs
surviving.

Because each row of `H_m` sums to one,

\[
\sum_j q_{Zj}=1-s_Ls_R+s_Ls_R=1.
\]

## Step 1: one-lineage switching

For one lineage,

\[
M=\begin{pmatrix}-m&m\\m&-m\end{pmatrix}.
\]

Diagonalization gives

\[
U(t)=e^{Mt}
=\frac12\begin{pmatrix}
1+e^{-2mt}&1-e^{-2mt}\\
1-e^{-2mt}&1+e^{-2mt}
\end{pmatrix}.
\]

The entering distribution of pair `12` is

\[
v_L(Z)=U(u_1)_{z_1,\cdot}\otimes U(u_2)_{z_2,\cdot},
\]

and pair `34` is

\[
v_R(Z)=U(u_3)_{z_3,\cdot}\otimes U(u_4)_{z_4,\cdot}.
\]

The tensor order is `00,01,10,11`.

## Step 2: one sister population

For a pair occupying one population, the transient states are
`00,01,10,11`. The row-vector generator is

\[
K=
\begin{pmatrix}
-(\lambda+2m)&m&m&0\\
m&-2m&0&m\\
m&0&-2m&m\\
0&m&m&-(\lambda+2m)
\end{pmatrix}.
\]

The loss from states `00` and `11` is coalescence. Thus the unnormalized
survivor vector is

\[
ell_Z=v_L(Z)e^{KT_L},\qquad r_Z=v_R(Z)e^{KT_R}.
\]

The survival probability is the row sum, and the branch coalescence
probability is

\[
P(\text{pair coalesces by }T)=1-v e^{KT}\mathbf1.
\]

### Closed form for e^{Kt}

Write

\[
S=x_{00}+x_{11},\quad M=x_{01}+x_{10},\quad
D=x_{00}-x_{11},\quad E=x_{01}-x_{10}.
\]

Then

\[
\frac{d}{dt}\binom{S}{M}
=
A\binom{S}{M},\qquad
A=\begin{pmatrix}-(\lambda+2m)&2m\\2m&-2m\end{pmatrix},
\]

while

\[
D(t)=e^{-(\lambda+2m)t}D(0),\qquad
E(t)=e^{-2mt}E(0).
\]

Let

\[
\Delta=\sqrt{\lambda^2+16m^2},\qquad
\mu_\pm=\frac{-(\lambda+4m)\pm\Delta}{2}.
\]

The symmetric block is

\[
e^{At}
=
\frac{e^{\mu_+t}(A-\mu_-I)-e^{\mu_-t}(A-\mu_+I)}{\Delta}.
\]

For each initial basis state, apply this two-by-two matrix to
`(S(0),M(0))`, apply the two scalar modes to `D(0),E(0)`, and reconstruct

\[
x_{00}=\frac{S+D}{2},\quad
x_{11}=\frac{S-D}{2},\quad
x_{01}=\frac{M+E}{2},\quad
x_{10}=\frac{M-E}{2}.
\]

This is the explicit closed form implemented in
`pair_kernel_closed_form`. It is algebraically equivalent to the
cosh/sinh eigendecomposition and is evaluated in spectral form for numerical
stability at large `m`. The unit tests compare it with
`scipy.linalg.expm(Kt)` to machine precision.

## Step 3: shared ancestral population and H_m

The exact shared-population kernel is formed from the same 16 background
states used by `msrcsim.analytic.compute_Hm`. For each background state
`z`, switching contributes rate `m` for each of its four lineages. Each
same-background pair contributes rate `lambda` to the absorbing topology
containing that pair.

If `Q` is the 16 by 16 transient generator and `R` the 16 by 3 absorption
matrix, then

\[
H_m(T_A)=
\left[e^{G T_A}\right]_{\rm transient,absorption}
+
\left[e^{G T_A}\right]_{\rm transient,transient}
\mathbf1\,(1/3,1/3,1/3),
\]

with

\[
G=\begin{pmatrix}Q&R\\0&0\end{pmatrix}.
\]

The second term resolves lineages that survive the structured interval in the
ordinary panmictic ancestral population. Consequently every row of `H_m`
sums to one.

## Relationship to the existing conditional H_m result.

Set `u_i=T_L=T_R=0`. Then `ell_Z` and `r_Z` are point masses and
`s_L=s_R=1`, so

\[
\mathbf q_Z=H_m(T_A)[z].
\]

For `Z=1010`, the symmetric calculation gives

\[
A_m(T_A)=
\frac{\lambda}{\lambda+2m}
\left[1-e^{-2(\lambda+2m)T_A}\right],
\]

and therefore

\[
H_m(1010)=
\left(
\frac{1-A_m}{3},
\frac{1+2A_m}{3},
\frac{1-A_m}{3}
\right).
\]

The new theorem is the species-tree extension: it propagates each lineage
through terminal branches, allows sister-pair coalescence on the two
species-tree branches, and only then invokes `H_m`.

## Closed-form symmetric quartet probabilities

For equal sister branches, the scalar coefficients and the four class formulas
are collected in `SCALAR_FORMULAS.md`. They follow by inserting the four
survivor vectors into the master identity. In particular, for the discordant
2:2 class,

\[
q_2-q_3=\alpha e^{-4mT},
\]

and the conditional topology-reversal condition is

\[
\alpha(d^2+e^2-2c^2)>1-s_D^2.
\]

The implementation checks these formulas against the current matrix evaluator
for every pattern over the requested numerical grid.

## All 16 patterns

Define

\[
F_L(ab)=
\left[U(u_1)_{a,\cdot}\otimes U(u_2)_{b,\cdot}\right]e^{KT_L},
\]

\[
F_R(cd)=
\left[U(u_3)_{c,\cdot}\otimes U(u_4)_{d,\cdot}\right]e^{KT_R},
\]

and

\[
\Phi(x,y)=(1-(x\mathbf1)(y\mathbf1))\mathbf e_1+
(x\otimes y)H_m(T_A).
\]

Then every row has the exact formula
`q_Z=Phi(F_L(z_1z_2),F_R(z_3z_4))`.

| Z | class | 2:2 split | exact conditional vector |
|---|---|---|---|
| 0000 | 4:0 | — | `Phi(F_L(00),F_R(00))` |
| 1111 | 4:0 | — | `Phi(F_L(11),F_R(11))` |
| 0001 | 3:1 | — | `Phi(F_L(00),F_R(01))` |
| 0010 | 3:1 | — | `Phi(F_L(00),F_R(10))` |
| 0100 | 3:1 | — | `Phi(F_L(01),F_R(00))` |
| 1000 | 3:1 | — | `Phi(F_L(10),F_R(00))` |
| 0111 | 3:1 | — | `Phi(F_L(01),F_R(11))` |
| 1011 | 3:1 | — | `Phi(F_L(10),F_R(11))` |
| 1101 | 3:1 | — | `Phi(F_L(11),F_R(01))` |
| 1110 | 3:1 | — | `Phi(F_L(11),F_R(10))` |
| 0011 | 2:2 | W1 | `Phi(F_L(00),F_R(11))` |
| 1100 | 2:2 | W1 | `Phi(F_L(11),F_R(00))` |
| 0101 | 2:2 | W2 | `Phi(F_L(01),F_R(01))` |
| 1010 | 2:2 | W2 | `Phi(F_L(10),F_R(10))` |
| 0110 | 2:2 | W3 | `Phi(F_L(01),F_R(10))` |
| 1001 | 2:2 | W3 | `Phi(F_L(10),F_R(01))` |

With equal terminal lengths, equal sister-branch lengths, and the same
symmetric rates everywhere, complementing all bits, swapping taxa within a
sister pair, and exchanging sister clades produce the expected equivalences.
The 4:0 rows form one class, all 3:1 rows form one class, and the 2:2 rows
form the three arrangement-split classes W1, W2, and W3. Without those equal
branch-length symmetries, the displayed formulas remain exact but rows that
are related only by a broken branch relabeling need not be numerically equal.

A 3:1 or 4:0 pattern is not generally `(1/3,1/3,1/3)`: the first term in
`Phi` favors `Q_1` whenever either sister pair has a nonzero chance to
coalesce before the shared population.

## MSC limit

For fast switching, the stationary background distribution is uniform. A
given pair is in the same background with probability one half, so its
effective coalescence rate is

\[
\frac12\lambda.
\]

The two sister branches contribute independently, hence

\[
s_Ls_R\longrightarrow
\exp\left[-\frac{\lambda}{2}(T_L+T_R)\right].
\]

The shared-population kernel becomes uniform across topologies, so with

\[
\tau=\frac{\lambda}{2}(T_L+T_R),
\]

\[
\boxed{
q_1\to1-\frac23e^{-\tau},\qquad
q_2\to\frac13e^{-\tau},\qquad
q_3\to\frac13e^{-\tau}.
}
\]

Under the simulator parameterization at `p=1/2`,
`lambda=1/(2Np)=1/N`, giving
`tau=(T_L+T_R)/(2N)`, the unrooted internal branch length in coalescent
units. This limit is independent of `Z`.

At `m=0`, the formula reduces to deterministic background persistence:
only same-background sister pairs can coalesce, and the shared `H_m`
kernel is evaluated from the unaltered entering pattern. At zero branch
lengths, the result is the `H_m` conditional formula. As sister branches
become long, `s_Ls_R\to0` and `q_Z\to e_1`. These cases are covered by
unit tests and direct numerical evaluation.

## Averaging over arbitrary forward sorting weights

Only after the conditional theorem is established, let

\[
\pi_z=P(Z=z),\qquad \sum_z\pi_z=1.
\]

Then

\[
\boxed{
\mathbf q_{\rm MSRC}=\sum_z\pi_z\mathbf q_z.
}
\]

For equal branch lengths, the conditional vector is constant within the five
symmetry classes 4:0, 3:1, W1, W2, and W3. Therefore the average depends on
the weights only through the total 4:0 mass, total 3:1 mass, and

\[
W_1=\pi_{0011}+\pi_{1100},\quad
W_2=\pi_{0101}+\pi_{1010},\quad
W_3=\pi_{0110}+\pi_{1001}.
\]

This is a proven reduction under the stated equal-branch symmetric model. For
unequal terminal or sister-branch durations, the reduction generally fails and
the full 16 weights must be retained. The weights can be obtained from the
forward Wright–Fisher dynamic program in
`msrc-sim/msrcsim/pattern_probabilities.py`; they are not assumed to be
uniform.

## Status of claims

* **Established analytically:** the one-lineage transition, the closed
  two-lineage kernel, the factorization through sister survival, probability
  conservation, the H_m special case, the five-class symmetry reduction under
  equal branch lengths, and the MSC limit.
* **Numerically validated:** closed kernel versus `scipy.linalg.expm`,
  all-16 probability conservation, H_m recovery, MSC convergence, and direct
  msrc-sim comparisons in the validation script.
* **Unproven simplification:** none is used in the theorem. Any scalar
  simplification beyond the five equal-branch symmetry classes must be checked
  against the matrix evaluator.
