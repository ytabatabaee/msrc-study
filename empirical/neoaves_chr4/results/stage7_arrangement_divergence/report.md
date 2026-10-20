# Neoaves chromosome 4 Stage 7 report

## Scope

This stage audits whether the independently defined structural-state tables support an arrangement-informed divergence contrast analogous to A(w) in the other empirical systems. The local trees are ML trees with substitution-unit branch lengths. They are not calibrated coalescence-time observations.

## Findings

- The structural labels are region-specific, multistate synteny/anchor classes and are independent of the local ML-tree topology used for the QQS tracks.
- Exact structural-species labels were mapped to six local-tree tips. The 47-taxon Jarvis mapping is retained as a separate provenance reference; it does not contain Ciconia maguari.
- Of 11 canonical regions, the numerically eligible exploratory regions are CR_01, CR_10, CR_11. The other regions have only singleton or otherwise insufficient state replication.
- The exploratory contrast is `A_distance = mean(different-state pairs) - mean(same-state pairs)` with a matched region-specific background `B_region`; `A_relative = A_distance/B_region`.
- These tracks compare deep species-level states and lack within-state population replication. They therefore do not support a primary population-level A(w) claim and should not be interpreted as direct evidence of recombination suppression.

## Figure recommendation

The topology overview is the primary Stage 7 figure: it keeps Columbea, N61, and N62 as separate source-defined quartet designs and shows both alternative contrasts. Region-specific distance tracks are shown separately as exploratory sensitivity results. The Stage 6 fixed-taxon global B(w) remains the appropriate arrangement-blind reference; no binary qH or universal A(w) was manufactured.
