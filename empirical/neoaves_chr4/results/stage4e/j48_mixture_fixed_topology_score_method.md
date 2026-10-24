# Fixed-topology quartet-score method

For fixed candidate trees T0 and NONCHR4, define the diagnostic quartet
universe U as every four-taxon subset on which the two candidates induce
different resolved unrooted topologies. For each input gene-tree record g and
quartet u in U, let I(g,u,T)=1 when g contains all four taxa, resolves u, and
matches candidate T; otherwise the observation is excluded. Repeated records
are repeated observations. Then `Q(T) = sum_g sum_u I(g,u,T)` and
`Delta-Q = Q(T0)-Q(T_NONCHR4)`.

The normalized difference is `Delta-Q / N`, where N is the number of resolved,
eligible gene-tree-by-diagnostic-quartet observations. The same N is used for
both candidates. Missing taxa and unresolved/polytomous quartets are excluded
symmetrically. This is an exact cancellation-aware fixed-candidate quartet
concordance score: quartets on which the candidates agree add equally to both
scores and cannot change Delta-Q.

The repository also contains `score_gene_records_bruteforce`, an intentionally
slow all-four-taxon reference implementation used on small fixtures. The
optimized implementation was checked against it exactly. For the 48-taxon
candidates, the diagnostic universe contains 3,597 quartets; a full production
score should be benchmarked on Bridges-2 before processing all mixtures.

This score is distinct from ASTRAL localPP, branch-local q1/q2/q3, the ASTRAL
topology optimization objective, and the published PNAS Delta-Q unless those
definitions are independently shown to coincide. The executable scorer does
not synthesize scores from completed ASTRAL annotations. ASTRAL/ASTER documents
`-C -c` as a scoring/annotation mode for a fully resolved candidate tree, but
the pinned Bridges-2 executable still needs a tiny known-answer validation
before its output can be used as this score.
