# Global conditional Moran simplex coverage: resolution status

## Conclusion

**FULL COVERAGE UNRESOLVED.**

For the original finite-`M` Moran MSRC model, the repository currently proves
vertex closure and documents exact conditional asymmetry witnesses, but it
does not prove or disprove

\[
\overline{\mathcal Q_{\mathrm{MSRC}\mid Z}}=\Delta_2.
\]

No numerical optimizer, local Jacobian, vertex construction, or external
mixture is being promoted to a global theorem.

## What is established

The exact conditional law is computed from the joint finite-M recursion,

\[
q_i(z;\theta)=J(z,i;\theta)/\sum_jJ(z,j;\theta),
\]

with the stochastic Moran count process, original switching/coalescence
generators, and binomial daughter-count inheritance. For both balanced and
pectinate rooted trees, paired patterns can approach the three simplex
vertices in valid parameter limits with positive conditioning probability at
each finite approximation.

For `Z=0101`, the exact pectinate solver also gives valid cases on both sides
of `r=1/2`, where

\[
r=q_1/(q_1+q_3).
\]

Therefore a simple universal constraint such as `q1 >= q3` is false. This does
not imply that all `r` values are available at fixed `d=q2`.

## Why the global proof is still missing

The required sector theorem would need to establish, for every feasible
`d`, that the level set

\[
\{r:q_2(z;\theta)=d\}
\]

contains the complete interval required by the corresponding dominance
sector. The existing finite-M operator gives this set through coupled matrix
exponentials and conditional normalization, but no monotonicity, connectedness,
or endpoint theorem has been proved for it.

The model is not known to be convex under parameter changes. An average of
three separately parameterized vertex constructions is an external mixture,
not a single Moran process. Likewise, the observed nonzero local Jacobians do
not establish global degree, boundary winding, or surjectivity.

No universal inequality excluding an open subset has been proved either. Thus
the available evidence does not justify either `FULL COVERAGE PROVED` or
`FULL COVERAGE DISPROVED`.

## Scope of the unresolved claim

| Question | Current status |
|---|---|
| One fixed arrangement pattern covers a dominance sector | Unresolved |
| Union over patterns covers the interior | Unresolved |
| All six strict quartet orderings occur | Unresolved globally |
| Exact finite-parameter interior coverage | Unresolved |
| Density in the closed simplex | Unresolved |
| Balanced and pectinate trees both satisfy global coverage | Unresolved |
| Conditional quartet-frequency overlap with introgression | Only conditional: it follows from a future coverage theorem; no current global result |

The union over arrangement patterns is permitted by the target definition, but
there is no proof that it is necessary or sufficient. No single-pattern
construction has yet been shown to cover a sector.

## Limiting constructions

The known limits use zero lower-branch durations, zero switching, and long
root duration. They prove vertex closure only. Extending them to edges would
require controlling the relative mass of the two discordant root resolutions;
the existing recursion provides no established independent control of that
mass. A conditioning event may remain positive along a finite approximating
sequence, but that fact alone does not determine its limiting quartet ratio.

## Introgression implication

There is no new exact or approximate introgression-indistinguishability
theorem. If a future proof established exact interior coverage, every
interior quartet vector from a specified four-taxon introgression model could
be matched conditionally at the quartet-frequency level. A density theorem
would give only approximate matching. Neither statement has been established
for the original Moran MSRC model here, and neither would imply equality of
full genomic or genealogical distributions.

## Precise next lemma

The most direct missing result is a global level-set lemma for one pattern:

> For every feasible `d` in the `Q2` dominance sector, the exact finite-M
> Moran conditional image at `Z=0101` contains every admissible `r` in the
> sector interval, or has closure containing that interval.

Proving this lemma, or proving a universal inequality that contradicts it,
would settle the global question. The current repository contains neither.
