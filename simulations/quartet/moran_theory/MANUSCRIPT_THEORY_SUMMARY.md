# Manuscript theory summary: finite-M MSRC quartet laws

This summary is suitable for adapting into the theory section. It records the
claims that are established by the current operator derivations and separates
them from numerical observations.

## Result 1: Exact quartet probabilities

For a finite chromosome-copy count, the arrangement-count process is a finite
Moran chain. Conditional genealogical states are propagated by the backward
structured-coalescent operators, while the Moran chain is coupled through the
dual forward/backward operator. Matrix-exponential propagation through the
species-tree branches gives exact quartet probabilities for the balanced and
pectinate histories represented by the implementation. Exactness here means
exactness for the stated hierarchical finite-state model, including its
binomial inheritance kernels and arrangement-origin convention.

## Result 2: Arrangement-conditional laws

For the pectinate tree, the conditional recursion retains the terminal
arrangement pattern as a joint probability mass. It returns

$$J[z,i]=P(Z=z,Q_i),\qquad
P(Q_i\mid Z=z)=J[z,i]/\sum_jJ[z,j],$$

for patterns with positive mass. The pattern-specific root-entry vectors are
not normalized before propagation. Summing over patterns recovers the validated
unconditional law. Zero-probability patterns have undefined conditional laws.

## Result 3: Conditional anomalous quartets

The exact finite-state anomaly criterion is

$$\max\{J[z,2],J[z,3]\}>J[z,1].$$

A proved scalar sufficient condition is available only in the fixed-frequency,
zero-switching special case. If frequency $p$ is fixed through all relevant
pre-root branches, labels do not switch, and `Z=0101` places taxa 1 and 3 in
background 0 while taxon 2 is in background 1, then

$$T_I>M(1-p)\log2$$

ensures $P(Q_2\mid Z=0101)>P(Q_1\mid Z=0101)$. The proof uses the exact
intermediate-branch event that taxa 1 and 3 coalesce. This result is
conservative and does not apply to the full Moran-integrated model merely from
the terminal pattern, because conditioning can alter ancestral count histories.

## Result 4: Symmetry and limits

Neutral complementation of all arrangement labels with $p\mapsto1-p$ preserves
the topology labels and maps the law at $z$ to the law at $1-z$. It therefore
preserves, rather than reverses, the discordant difference
$\Delta=q_2-q_3$. Sign reversal requires an appropriate taxon permutation.
For the simulator's genic-selection convention, the corresponding label swap
uses $s\mapsto-s/(1+s)$. Rapid switching has a fixed-frequency averaged pair
hazard $1/M$ at leading order, but a complete conditional convergence theorem
for the coupled Moran model remains open because of terminal-state boundary
layers.

The numerical base example has a conditional anomaly at `Z=0101`, whereas the
tested symmetric settings have no unconditional anomaly after marginalizing
over patterns. These observations do not establish a universal unconditional
bound or a causal claim about recombination suppression.
