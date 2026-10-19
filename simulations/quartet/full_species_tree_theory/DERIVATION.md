# Exact four-taxon MSRC species-tree quartet distribution

This note derives the quartet distribution for the balanced species tree
`S=((1,2),(3,4))`. The event rules are those implemented in
`msrc-sim/msrcsim/structured_coalescent.py`: a lineage switches background
with rates

\[
m_{01}=\rho p,\qquad m_{10}=\rho(1-p),
\]

and two lineages coalesce only when they are in the same population and have
the same background. Their coalescence rate is

\[
\lambda_0=\{2N(1-p)\}^{-1},\qquad
\lambda_1=\{2Np\}^{-1}.
\]

The terminal sampled arrangement is `Z=z=(z_1,z_2,z_3,z_4)`.

## State space and exact matrix expression

Before the first coalescence there are four surviving sampled lineages. A
transient state is

\[
x=(b_1,b_2,b_3,b_4; a_1,a_2,a_3,a_4),
\]

where `b_i` is the current species-tree population and `a_i` is the current
background. Within one constant-rate epoch, the 16 background vectors form the
transient state space because the population assignment is fixed by the
species-tree epoch. The three absorbing states are

\[
Q_1=12|34,\qquad Q_2=13|24,\qquad Q_3=14|23.
\]

For a transient row state `a`, the generator has the following entries.

* A switch of lineage `i` from 0 to 1 has rate `m01`; a switch from 1 to 0
  has rate `m10`.
* For every pair `(i,j)` in the same population and with `a_i=a_j`, add
  `lambda_{a_i}` to the absorbing column corresponding to the topology
  containing that pair.
* The diagonal is minus the sum of all switch and coalescence rates.

Write this block generator as

\[
G=\begin{pmatrix}Q&R\\0&0\end{pmatrix},
\]

where `Q` is 16 by 16 and `R` is 16 by 3. For an epoch of length `d`,

\[
(\alpha',q')=(\alpha,q)\exp(dG).
\]

At a species-tree boundary, every lineage in a child population is mapped to
its parent population, while its background is unchanged. This is a
deterministic map on transient states. For the balanced tree the epochs are:

1. terminal branches, with populations `(1,2,3,4)`;
2. the two sister branches, with populations `(12,12,34,34)`;
3. the shared root, with populations `(r,r,r,r)`.

Thus, for finite epoch lengths `d_1,d_2`,

\[
q(z)=\alpha_z K_{\rm terminal}K_{\rm sister}K_{\rm root},
\]

where each `K` is the augmented matrix exponential with the appropriate
deterministic boundary map. The root is infinite in the implementation. If
`Q_r` and `R_r` are its transient and absorbing blocks, its contribution is

\[
\alpha_r\int_0^\infty e^{Q_rt}R_r\,dt
=\alpha_r(-Q_r)^{-1}R_r.
\]

This is the exact matrix answer. A scalar expression is generally unavailable
because the boundary maps mix several 16-state background kernels.

## The four-taxon species-tree probabilities

For each of the 16 terminal patterns, the implementation evaluates

\[
P(Q_j\mid Z=z,S)=q_j(z).
\]

If the forward frequency-history calculation supplies `P(Z=z)`, then

\[
P(Q_j\mid S)=\sum_{z\in\{0,1\}^4}P(Z=z)q_j(z).
\]

The Wright–Fisher dynamic program in `msrc-sim/msrcsim/pattern_probabilities.py`
is the appropriate source for the weights. The biological bridge is therefore

\[
(N,s,\rho,\text{origin time},S)\longrightarrow P(Z)
\longrightarrow P(Q\mid Z,S)\longrightarrow P(Q\mid S).
\]

The quartet CTMC is conditional on the branchwise piecewise-constant values of
`p` and `N`; it does not replace the simulator's forward population process.

## All 16 patterns and symmetry classes

The table uses `n_1` for the number of 1 backgrounds and records the pattern's
species-tree pair relationship. Complementing every bit preserves the model
after exchanging labels 0 and 1 when `p` is replaced by `1-p`.

| pattern | class | n_1 | sister-pair relation |
|---|---:|---:|---|
| 0000 | 4:0 | 0 | both concordant pairs equal |
| 1111 | 4:0 | 4 | both concordant pairs equal |
| 0001 | 3:1 | 1 | 1/2 equal, 3/4 mixed |
| 0010 | 3:1 | 1 | 1/2 mixed, 3/4 equal |
| 0100 | 3:1 | 1 | 1/2 mixed, 3/4 equal |
| 1000 | 3:1 | 1 | 1/2 mixed, 3/4 equal |
| 0111 | 3:1 | 3 | 1/2 mixed, 3/4 equal |
| 1011 | 3:1 | 3 | 1/2 mixed, 3/4 equal |
| 1101 | 3:1 | 3 | 1/2 equal, 3/4 mixed |
| 1110 | 3:1 | 3 | 1/2 mixed, 3/4 equal |
| 0011 | 2:2 | 2 | both concordant pairs mixed |
| 1100 | 2:2 | 2 | both concordant pairs mixed |
| 0101 | 2:2 | 2 | each sister pair mixed |
| 1010 | 2:2 | 2 | each sister pair mixed |
| 0110 | 2:2 | 2 | both concordant pairs mixed |
| 1001 | 2:2 | 2 | both concordant pairs mixed |

The table is a bookkeeping classification, not an assertion that all rows in
a class have equal quartet probabilities. Equality requires the corresponding
taxon/background relabeling to preserve both the species tree and the rate
parameters. In particular, 4:0 and 3:1 patterns can still favor `Q_1` on the
full species tree because the `(1,2)` and `(3,4)` boundaries create a
concordant opportunity.

## Required limits

Let `tau` be the unrooted internal species-tree branch length in coalescent
units. In the rooted ultrametric representation used here, two equal sister
branches of length `T` correspond to `tau=2T/(2N)=T/N`. As
`rho` tends to infinity, each lineage rapidly reaches the stationary
background distribution `(1-p,p)`. For any pair, the averaged same-background
coalescence hazard is

\[
(1-p)^2\lambda_0+p^2\lambda_1
=\frac{1-p}{2N}+\frac{p}{2N}=\frac1{2N}.
\]

All pairs therefore have the same ordinary coalescent hazard after entering a
common population. No lineage switches during the limiting boundary layer, so
the structured model converges to the MSC:

\[
P(Q_1)=1-\frac23e^{-\tau},\qquad
P(Q_2)=P(Q_3)=\frac13e^{-\tau}.
\]

With one background, or with a fixed/lost rearrangement that removes switching
and background dependence, the same unstructured coalescent argument gives the
MSC directly. Probability conservation follows from the augmented generator's
row sums and the root absorption solve.

## Relationship to the previous `H_m` conditional formula

The earlier `compute_Hm` calculation starts after all four lineages have
already entered one shared structured population. It propagates the same
16-state background CTMC for a finite time and assigns the residual mass
uniformly to the three topologies. The helper `shared_population_kernel` in
`quartet_theory.py` reproduces that calculation exactly. For `z=1010`,
`m01=m10=m`, and `lambda0=lambda1=lambda`, it gives

\[
A_m(t)=\frac{\lambda}{\lambda+2m}
[1-e^{-2(\lambda+2m)t}],
\]

\[
P(13|24\mid Z=1010)=\frac{1+2A_m(t)}3,
\qquad
P(12|34\mid Z=1010)=P(14|23\mid Z=1010)=\frac{1-A_m(t)}3.
\]

This is a special case of the shared-population kernel. It is not the
species-tree distribution, because it omits terminal and sister-branch
propagation and their boundary opportunities for coalescence.

## Time-varying frequencies

For a time-varying `p(t)`, replace each constant-rate exponential by the
time-ordered solution of `dP/dt=P G(t)`. Piecewise-constant approximations
give products of ordinary exponentials, which is also the natural numerical
route for frequencies generated by the simulator's forward history.
