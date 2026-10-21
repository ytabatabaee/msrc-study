# Derivation: the shared-population Moran-integrated kernel

This result is the exact finite-
`M` kernel for the hierarchical model currently implemented by `msrc-sim`.
It averages a forward Moran frequency process and a backward structured
coalescent conditional on the realized frequency history. It is a shared-
population result and is not yet the full species-tree theorem.

## 1. Forward Moran process

Let `M=2 Ne` and let `K_t` count arrangement-1 chromosome copies. Biological
time runs from the older boundary to the younger boundary. On `K=0,...,M`,
the row generator has

```text
Q_M[k,k+1] = q_plus(k),   Q_M[k,k-1] = q_minus(k),
Q_M[k,k]   = -(q_plus(k)+q_minus(k)).
```

For neutrality, matching `msrcsim/moran.py`,

```text
q_plus(k) = q_minus(k) = k(M-k)/(2M).
```

The endpoints are absorbing. The optional selection parameter uses the same
genic Moran fitness transform as `msrcsim.moran_rates`.

## 2. Backward genealogy conditional on `K_t=k`

For a fixed count, put `p=k/M`. A lineage switches `0 -> 1` at `rho p` and
`1 -> 0` at `rho(1-p)`, where `rho` is the effective background-switching
intensity passed to `structured_coalescent.py`.

The simulator's finite floor gives

```text
c1(k) = 1/max(k,1),   c0(k) = 1/max(M-k,1).
```

The 16 transient states use the same `itertools.product((0,1), repeat=4)`
ordering as `msrcsim/analytic.py`. For every background state `z`, `S(k)`
contains switching and the negative total switching/coalescence loss on its
diagonal. The absorption matrix has

```text
R[z,1] = 1[z1=z2] c_z1 + 1[z3=z4] c_z3,
R[z,2] = 1[z1=z3] c_z1 + 1[z2=z4] c_z2,
R[z,3] = 1[z1=z4] c_z1 + 1[z2=z3] c_z2.
```

The augmented generator is `G4(k)=[[S(k),R(k)],[0,0]]`. Every row sums to
zero. The three absorbing columns represent `12|34`, `13|24`, and `14|23`.

## 3. Forward/backward ordering and the dual operator

Suppose a forward path visits `k1,...,kn` over durations `Delta1,...,Deltan`.
A young-end genealogy sees the reverse chronological product

```text
U_back = exp(G4(kn) Deltan) ... exp(G4(k1) Delta1).
```

Therefore its transpose acts in forward chronological order:

```text
U_back^T = exp(G4(k1)^T Delta1) ... exp(G4(kn)^T Deltan).
```

The averaged implementation uses the sparse operator

```text
A4 = Q_M kron I_19
     + blockdiag(G4(0)^T, ..., G4(M)^T).
```

This is a dual/Feynman--Kac operator used with row-vector propagation. It
should not be described as a literal forward biological genealogy CTMC.

## 4. Terminal sampling, residual mass, and probabilities

At the young boundary, conditional on `K_T=k`, the four sampled arrangements
are independent. For pattern `z`, with `n1(z)=sum(z)`,

```text
beta_z(k) = (k/M)^n1(z) (1-k/M)^(4-n1(z)).
```

The terminal column `b_z` contains this value in the joint state `(k,z)` and
zero elsewhere. `b_all=sum_z b_z`.

At the older boundary, unresolved mass is assigned uniformly by the ordinary
panmictic ancestor. For topology `r`, `h_r` is one third on all 16 transient
states, one on absorbed state `Q_r`, and zero on the other absorbed states.
For the older count distribution `pi_old`, set `a_r=pi_old kron h_r`.

The exact joint probabilities are

```text
J[z,r] = a_r exp(A4 T) b_z.
```

Thus `q_r=sum_z J[z,r]`, `P(Z=z)=sum_r J[z,r]`, and
`P(Q_r|Z=z)=J[z,r]/P(Z=z)` when the denominator is positive.

## 5. Checks and the fixed-frequency limit

The implementation checks probability conservation, the direct terminal
Moran count law, zero time, zero switching, and explicit forward/backward
operator ordering. With Moran evolution disabled and `K=M/2` fixed,

```text
m = rho/2,   lambda = 2/M = 1/Ne,
A_m(T) = lambda/(lambda+2m) [1-exp(-2(lambda+2m)T)].
```

For `Z=1010`, this gives

```text
P(Q2|1010) = (1+2 A_m)/3,
P(Q1|1010) = P(Q3|1010) = (1-A_m)/3.
```

This establishes that the old `A_m` formula is the fixed-`p=1/2` special
case. It is not the general finite-M Moran-integrated law. In particular,
`p=1/2` does not imply `T_L=T_R=0`; the latter is the shared-population
limit, and the old formula additionally holds the shared-population frequency
fixed at `p(t)=1/2`.

The calculation is exact for the stated hierarchical model. It does not claim
that this operator has been derived as the microscopic ancestral dual of the
same Moran reproduction events.

