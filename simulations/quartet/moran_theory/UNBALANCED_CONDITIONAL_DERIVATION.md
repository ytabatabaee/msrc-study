# Exact arrangement-conditional quartet law for `(((1,2),3),4)`

## Scope

This derivation extends the finite-state hierarchical Moran recursion to the
terminal arrangement pattern

\[
Z=(Z_1,Z_2,Z_3,Z_4)\in\{0,1\}^4,
\qquad Q_1=12|34,;Q_2=13|24,;Q_3=14|23.
\]

It is exact for the model implemented by `msrc-sim`: a forward finite-\( M \)
Moran count process, binomial daughter initialization conditional on the
shared parental count, and a backward structured genealogy conditional on the
realized count history. It is not a general theorem for every biological
Moran process and does not assert a microscopic ancestral dual for the
forward Moran reproduction events.

## Pattern-specific terminal vectors

At a terminal branch with count \( k \), write \( p=k/M \). For a fixed
sampled arrangement bit \( z_j \), the young-boundary vector is

\[
u_{j,z_j}(k,b)=
1_{\{b=z_j\}}p^{z_j}(1-p)^{1-z_j}.
\]

The factor is the probability of sampling the terminal arrangement and the
indicator records the arrangement carried by the sampled lineage. Summing over
the two possible sampled bits gives the ordinary unconditional vector:

\[
u_{j,0}(k,b)+u_{j,1}(k,b)=1_{\{b=0\}}(1-p)+1_{\{b=1\}}p.
\]

Each vector is propagated through its terminal branch with the existing dual
one-lineage operator and then through the appropriate inheritance kernel.
For a pair of terminal bits \( z_1z_2 \), the pair vector is constructed by
conditional multiplication given the sister count, exactly as in the
unconditional recursion. This retains the shared parental count and the joint
sampling mass.

## Three-lineage stage and early mass

Let \( L_{12,z_1z_2}(k_I,s) \) be the six-state pair result at the younger
end of the intermediate branch and \( L_{3,z_3}(k_I,b) \) the one-lineage
result for taxon 3. The younger-boundary three-lineage vector is

\[
Y_{z}(k_I,Q_1)=
\big[L_{12,z_1z_2}(k_I,C0)+L_{12,z_1z_2}(k_I,C1)\big]
\sum_b L_{3,z_3}(k_I,b),
\]

and, for an unresolved state \( (b_1,b_2,b_3) \),

\[
Y_z(k_I,b_1b_2b_3)=
L_{12,z_1z_2}(k_I,b_1b_2)L_{3,z_3}(k_I,b_3).
\]

The first expression is essential: a sister-pair coalescence has already
fixed \( Q_1 \), but the terminal pattern still includes taxon 3. The
three-lineage vector is propagated with \( A_3 \) and inherited into the
root count. Let \( W_z(k_R,Q_i) \) denote its absorbing columns.

The early joint contribution must also retain taxon 4's sampling mass. If
\( L_{4,z_4}(k_R,b) \) is the root-count conditional one-lineage result,

\[
c_i(z)=\sum_{k_R}\pi_R(k_R)W_z(k_R,Q_i)
\sum_b L_{4,z_4}(k_R,b).
\]

Dropping this factor would incorrectly treat \( Z_4 \) as independent of
the root count.

## Root-entry vector and joint law

For the 16 unresolved root-background states,

\[
F_z(k_R,b_1b_2b_3b_4)=
W_z(k_R,b_1b_2b_3)L_{4,z_4}(k_R,b_4),
\]

and the three absorbed root states are zero. This is an unnormalized joint
probability mass. It must not be normalized separately for each pattern.
Using the existing four-lineage operator and residual vectors gives

\[
J[z,i]=c_i(z)+
(\pi_{\rm old}\otimes h_i)^T
e^{A_4T_A}F_z.
\]

The implementation evaluates all root patterns in one batched
`expm_multiply` call. It returns `joint[z,i]`, `early_topology_joint[z,i]`,
`root_topology_joint[z,i]`, and the unnormalized root-entry arrays.

## Probability identities

The pattern-specific terminal vectors partition the ordinary terminal vector.
Because inheritance kernels are row-stochastic and the conditional genealogy
operators partition probability into transient and absorbing states, summing
over patterns gives the unconditional recursion. Therefore

\[
\sum_{z,i}J[z,i]=1,
\qquad \sum_zJ[z,i]=P(Q_i),
\qquad \sum_iJ[z,i]=P(Z=z).
\]

For \( P(Z=z)>0 \),

\[
P(Q_i\mid Z=z)=J[z,i]/P(Z=z).
\]

The implementation returns `NaN` for patterns with zero probability. The test
suite verifies all identities numerically, including the equality of the
conditional marginal to the previously validated unconditional result.

## Conditional asymmetry and anomalies

Define

\[
\Delta(z)=P(Q_2\mid z)-P(Q_3\mid z),
\qquad
\Gamma(z)=\max(P(Q_2\mid z),P(Q_3\mid z))-P(Q_1\mid z).
\]

These are conditional quantities. A large value must be interpreted together
with \( P(Z=z) \). The compact analytical exploration found, for the base
condition \( M=10, p_0=0.5, \rho=0.17, T_S=0.3,T_I=0.2,T_A=0.5 \),

\[
P(Z=0101)=P(Z=0110)=P(Z=1001)=P(Z=1010)=0.0443121.
\]

For `0101`,

\[
(q_1,q_2,q_3)=(0.275862,0.449849,0.274289),
\]

so \(\Delta=0.175560\) and \(\Gamma=0.173987\). The complementary
patterns exchange \( q_2 \) and \( q_3 \). Thus the model can produce a
conditional anomalous quartet, although the corresponding unconditional law
for this condition remains \( q_1,q_2,q_3)=(0.356291,0.321855,0.321855) \).

This is a numerical result for the stated parameter regime. A sufficient
symmetry condition for \( q_2(z)=q_3(z) \) is that the conditional process
and terminal pattern are invariant under the taxon permutation exchanging the
two discordant resolutions. In particular, when all sampled lineages are in
one common background and switching does not create a distinction, the
ordinary coalescent exchangeability of the two discordant resolutions gives
\( q_2=q_3 \). The unbalanced tree does not make every pattern invariant;
crossed patterns can break this symmetry.

## Balanced comparison

The balanced comparison uses the existing conditional balanced recursion with
the same root count, branch counts, terminal counts, and root-to-tip depths.
The left sister split is assigned duration \( T_S \), the right sister
split duration \( T_S+T_I \), and terminal branches are zero. This is a
time-depth-matched comparison, not a claim that the two trees have identical
ancestral-stage histories.

Under the base condition and `Z=0101`, the balanced result is

\[
(q_1,q_2,q_3)=(0.287204,0.430903,0.281893),
\]

which also has conditional discordance and an anomalous margin, but with
smaller \(\Delta\) and \(\Gamma\) than the pectinate result. This comparison
supports a difference between the matched stage structures, while not by
itself isolating rooted shape from the differing branch decomposition.

## Independent validation and limitations

Fresh `msrc-sim` histories were simulated for the base, zero-switching, and
high-switching conditions, with 4,000 genealogies per condition. The joint
table has 48 cells per condition. The maximum absolute cell error was
0.00663, and 99.31% of cells were within three multinomial Monte Carlo
standard errors. Sparse-cell comparisons should be interpreted cautiously;
the validation is a targeted smoke check rather than a full parameter sweep.

The result establishes exact joint and conditional laws for the hierarchical
finite-state model. It does not establish that a conditional anomalous quartet
is a genome-wide species-tree prediction after averaging over unobserved
arrangement states, nor does it identify introgression or any other biological
mechanism beyond the model's specified genealogy process.
