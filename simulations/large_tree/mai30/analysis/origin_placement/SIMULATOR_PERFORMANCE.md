# Optimized simulator performance

The authoritative backward implementation now retains the original function as
`simulate_genealogy_reference` and uses `simulate_genealogy` as a grouped exact
implementation. It aggregates switching and same-background coalescence
hazards by `(population, background)` rather than rebuilding a lineage-level
pair list at every event. It makes one frequency-history lookup per occupied
population and samples a lineage or unordered pair uniformly within the chosen
hazard group.

This changes random-number consumption but preserves the target rates,
species-boundary transitions, structured-origin boundary, and ordinary
panmictic fallback. It is not an approximation or a replacement population
process.

On Mai replicate 007, one 30-taxon genealogy at rho=0.02 took 0.276 seconds
with the reference implementation and 0.128 seconds with the grouped
implementation (10-replicate benchmark; approximately 2.15x faster). The
reference and optimized implementations produced 29 coalescences per
30-taxon genealogy and comparable event counts. A 1,500-replicate four-taxon
topology comparison passed the distributional regression threshold.

The remaining bottleneck for the origin-placement experiment is forward
Wright--Fisher history construction for long branches and doubled root
extensions, followed by very high backward event counts at rho=0.1--1.0.
Those cases were not silently approximated. The 1,000-history all-tree
persistence run was stopped before completion and is recorded as unexecuted.
