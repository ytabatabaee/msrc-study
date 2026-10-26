# Analytical conditions for anomalous MSRC quartets

## Definitions

For the pectinate species tree `(((1,2),3),4)`, write

\[
Q_1=12|34,\qquad Q_2=13|24,\qquad Q_3=14|23.
\]

For an arrangement pattern \(z\) and parameters \(\theta\), define

\[
q_i(z;\theta)=P_\theta(Q_i\mid Z=z),\quad
\Delta(z;\theta)=q_2(z;\theta)-q_3(z;\theta),
\]

and

\[
\Gamma(z;\theta)=\max\{q_2(z;\theta),q_3(z;\theta)\}-q_1(z;\theta).
\]

The conditional quartet is anomalous precisely when \(\Gamma>0\). This is a
conditional statement. The unconditional quantities are obtained by summing
the joint masses over all patterns and need not have the same behavior.

## Exact operator inequality

The conditional implementation returns the joint mass

\[
J[z,i]=P(Z=z,Q_i)=c_i(z)+
(\pi_{\rm old}\otimes h_i)^T e^{A_4T_A}F_z,
\]

where \(F_z\) is not normalized by its pattern probability. If
\(P(Z=z)>0\), then

\[
q_2(z)>q_1(z)\iff J[z,2]-J[z,1]>0,
\]

with the exact difference

\[
J[z,2]-J[z,1]=c_2(z)-c_1(z)+
\{\pi_{\rm old}\otimes(h_2-h_1)\}^T e^{A_4T_A}F_z.
\]

The analogous expression with \(h_3-h_1\) decides \(q_3>q_1\). This is a
necessary and sufficient condition for the finite-\(M\) model. It is useful
because the early term and root term can be inspected separately, but the
matrix-exponential term generally has no scalar monotonicity in \(\rho\) or
the branch lengths.

## A proved zero-switching sufficient condition

Consider the fixed-frequency, zero-switching limit: the arrangement frequency
is constant at \(p\in(0,1)\), \(\rho=0\), and the chromosome count is \(M\).
For background 0 and 1, the same-background pair rates are

\[
c_0=\frac{1}{M(1-p)},\qquad c_1=\frac{1}{Mp}.
\]

For \(z=0101\), taxa 1 and 2 differ, so no sister-branch coalescence is
possible. On the intermediate branch, taxa 1 and 3 are the only matching pair
among taxa 1,2,3. Therefore

\[
P(Q_2\mid z)\ge 1-e^{-c_0T_I},
\qquad
P(Q_1\mid z)\le e^{-c_0T_I}.
\]

Consequently,

\[
\boxed{T_I>M(1-p)\log 2\quad\Longrightarrow\quad
P(Q_2\mid Z=0101)>P(Q_1\mid Z=0101).}
\]

This is a rigorous conservative sufficient condition. It is not necessary;
the exact operator can produce an anomaly at much shorter durations because
root-stage events and competing background configurations contribute additional
mass. The corresponding results are:

\[
\begin{array}{c|c|c}
z & \text{winning intermediate pair} & \text{sufficient condition}\\
\hline
0101,1010 & 13 & T_I>M(1-p)\log2\text{ for the matching background}\\
0110,1001 & 23 & T_I>M(1-p)\log2\text{ or }Mp\log2\text{ as appropriate}.
\end{array}
\]

More explicitly, the bound is \(M(1-p)\log2\) when the matching pair is in
background 0 and \(Mp\log2\) when it is in background 1. It is independent
of \(T_S\) and \(T_A\) because the event used in the lower bound occurs before
lineage 4 enters and fixes the quartet topology.

For `0000` and `1111` in this same fixed-background limit, all lineages occupy
one background and the process reduces to an ordinary structured coalescent on
the species tree. The two discordant topologies are exchangeable and the
species-tree topology is at least as probable, so no conditional anomaly occurs.
The remaining patterns, including `0011` and `1100`, require the exact
operator or additional case-specific bounds.

## Symmetry results

For neutral dynamics, complementing every arrangement bit and replacing
\(p\) by \(1-p\) exchanges background 0 and 1 while preserving all
coalescence and switching rates. Thus the conditional law is carried to the
complementary pattern with the same topology probabilities. With genic
selection, the corresponding label-swap transformation also changes the
fitness parameter from \(s\) to \(-s/(1+s)\), when \(s>-1\).

If a measure-preserving taxon permutation preserves the species-tree stages,
the branch parameters, and the specified pattern event while exchanging the
two discordant resolutions, then \(q_2=q_3\). The crossed patterns in the
pectinate tree do not satisfy this condition, which is why `0101` and its
relatives can have \(\Delta\ne0\). Complementary crossed patterns exchange
the sign of \(\Delta\).

Complement symmetry at \(p=1/2\) restores unconditional \(q_2=q_3\) when
patterns are summed in complementary pairs. It does not force unconditional
\(q_1\) to exceed the discordant probabilities without an additional
species-tree argument. The numerical calculations in this repository found no
unconditional anomaly in the tested grid; this is reported as an observation,
not as a general impossibility theorem.

## Short-branch and fast-switching limits

For any fixed pattern, every branch propagation is an analytic finite-matrix
function. Hence, for a branch operator \(A\),

\[
e^{AT}=I+AT+O(T^2),
\]

with a standard finite-dimensional remainder bound
\(\|R_2\|\le e^{\|A\|T}\|A\|^2T^2/2\). This gives an exact first derivative of
each joint-mass difference by replacing one branch exponential at a time by
its generator. It supplies a local boundary test but does not establish a
global threshold.

In the fixed-frequency model, as \(\rho\to\infty\), background states mix
rapidly to \((1-p,p)\). Averaging the pair hazards gives

\[
p^2c_1+(1-p)^2c_0=1/M,
\]

so the leading fast-switching limit is the ordinary MSC with pair rate
\(1/M\), and arrangement-conditioned asymmetry disappears. The finite-M
Moran model with a changing count process has an analogous timescale-separation
expectation, but a full singular-perturbation proof for the coupled operator is
outside this analysis.

## Thresholds and numerical findings

For the base condition \(M=10,p=1/2,\rho=0.17,T_S=0.3,T_I=0.2,T_A=0.5\),

\[
P(Q\mid Z=0101)=(0.275862,0.449849,0.274289),
\]

so \(\Delta=0.175560\) and \(\Gamma=0.173987\), with
\(P(Z=0101)=0.044312\). In the exact operator calculation, the numerical
root of \(q_2-q_1=0\) for this pattern is approximately \(\rho^*=2.75\) at
\(T_I=0\), decreasing to approximately 2.21 at \(T_I=0.8\). These roots are
reported as numerical crossings; uniqueness and global monotonicity were not
proved. The anomaly region is connected over the plotted \((\rho,T_I)\) grid,
but that is also a numerical observation.

The matched balanced comparison uses the same root-to-tip depths, with sister
durations \((T_S,T_S+T_I)\). For `0101` under the base condition it gives
\((0.287204,0.430903,0.281893)\), also anomalous but with smaller margins.
This does not isolate rooted shape from ancestral-stage decomposition because
the two histories are not identical after matching time depths.

## Scope

The proved sufficient condition and symmetry statements apply only under their
explicit fixed-frequency or neutral assumptions. The full finite-M
Moran-integrated model is governed by the exact joint operator inequality, and
its threshold surfaces may be nonmonotone or disconnected. A conditional
anomaly is an arrangement-stratified prediction; after marginalization it may
cancel, as it does for the tested symmetric base example. No result here says
that recombination suppression necessarily causes anomalous quartets or
identifies introgression.
