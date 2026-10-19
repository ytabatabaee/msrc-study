# Manuscript theorem and corollary

**Theorem (four-taxon MSRC quartet law).** Consider the ultrametric species
tree `S=((1,2),(3,4))`, one sampled lineage per taxon, and a piecewise-constant
MSRC history. On each population branch, background switching has rates
`m01=rho p` and `m10=rho(1-p)`, and same-background coalescence has rates
`lambda_0={2N(1-p)}^{-1}` and `lambda_1={2Np}^{-1}`. Conditional on the
terminal background pattern `Z=z`, the quartet probabilities are obtained by
propagating the 16-state transient background distribution through the
terminal and sister population epochs with augmented matrix exponentials,
applying deterministic species-boundary maps, and multiplying by the root
absorption matrix `(-Q_root)^(-1) R_root`. The resulting three-vector is
`P(Q | Z=z,S)` and sums to one.

**Corollary (unconditional quartet law).** If `w_z=P(Z=z)` is obtained from
the forward population-frequency process, then

`P(Q_j | S) = sum_z w_z P(Q_j | Z=z,S)`.

This separates forward arrangement prevalence from backward lineage sorting.
When `rho -> infinity`, or when only one background is available, the
effective pairwise coalescence rate is `1/(2N)` and the law reduces to the
ordinary MSC:

`P(Q_S)=1-(2/3) exp(-tau)` and `P(Q_D1)=P(Q_D2)=(1/3) exp(-tau)`,

where `tau` is the unrooted internal branch length in coalescent units.
