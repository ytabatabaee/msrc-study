# Manuscript theorem and corollaries

**Theorem (symmetric four-taxon MSRC quartet law).** For
`S=((1,2),(3,4))`, one sampled lineage per species, fixed terminal pattern
`Z=z`, symmetric switching `m01=m10=m`, and symmetric coalescence
`lambda0=lambda1=lambda`, let `ell_z` and `r_z` be the unnormalized
four-state survivor distributions of pairs `12` and `34` after their
species-tree population branches. Then

`q_z = (1-s_L s_R)(1,0,0) + (ell_z kron r_z) H_m(T_A)`,

where `s_L=ell_z 1`, `s_R=r_z 1`, and `H_m(T_A)` is the exact shared
structured-population kernel. The vector sums to one.

**Corollary (MSC limit).** As `m -> infinity`, the conditional distribution
becomes independent of `z` and

`q_1=1-(2/3) exp(-tau)`, `q_2=q_3=(1/3) exp(-tau)`,

with `tau=lambda(T_L+T_R)/2`. Under the simulator's `p=1/2`
parameterization, `lambda=1/N), so `tau=(T_L+T_R)/(2N)`.

**Corollary (arbitrary forward sorting).** For arbitrary forward pattern
weights `pi_z=P(Z=z)`,

`q_MSRC = sum_z pi_z q_z`.

With equal terminal and sister branch lengths, the conditional laws reduce to
the five classes 4:0, 3:1, W1, W2, and W3. The corresponding aggregate
weights are the total 4:0 mass, total 3:1 mass, and
`W1=pi_0011+pi_1100`, `W2=pi_0101+pi_1010`, and
`W3=pi_0110+pi_1001`.
