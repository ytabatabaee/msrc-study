# Stage 13 Phlag benchmark methods

Phlag was kept at upstream commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6` (version 0.2.0) with its default model settings. The finite bundled-example check uses the previously recorded reduced-iteration reproduction (`-L 1 -l 5 --n-replicates 10`); a fresh full-default toy run was not completed within the practical runtime window. For Atlantic cod, the unchanged Stage 11 runs were copied into `raw_results/` and reclassified from their diagnostics. For the other systems, the benchmark stops at the input audit when the local-tree and reference-tree taxon sets differ.

Phlag uses a shared DendroPy taxon namespace for the species tree and local gene trees. Collapsing individuals, haplotypes, or pseudo-t-haplotype tips to population labels would change the input data and is therefore not used. A numerical failure is never scored as a biological no-anomaly result.

The Phlag objective detects genealogical deviations/model violations. It is not itself a recombination-suppression or inversion-boundary estimator.
