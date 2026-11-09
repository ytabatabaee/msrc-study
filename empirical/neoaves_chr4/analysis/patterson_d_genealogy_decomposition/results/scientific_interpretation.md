# Neoaves N61 genealogy-based Patterson's D decomposition

## Scope and validation
The analysis used the existing checksum-verified Stiller alignment archive without redownloading it. Exact chr4 alignment-member mappings were used; 6247 alignment members were read and callable sites were calculated for the four taxa in each design. The named gene-tree file was parsed for 6247 loci and each design was classified by a root-independent four-point induced quartet.

The induced quartet labels are q1=(P1,P2)|(P3,O), q2=(P1,P3)|(P2,O), q3=(P1,O)|(P2,P3). Under outgroup polarization, ABBA is q3 and BABA is q2. Internal lengths are induced quartet internal edges in ML substitution units, not coalescence times.

## Main quantitative result
For N61, the original observed locus-level shift is ΔD_site = -0.177615. The genealogy predictions are ΔD_topology = -0.878071 and ΔD_tree = -0.889805. The exact values, exposures, and denominators are in `observed_vs_predicted_D.tsv`.
The symmetric decomposition attributes -0.642899 of the predicted shift to topology/exposure composition and -0.246906 to conditional internal-branch lengths; their sum is -0.889805, equal to the predicted D_tree shift -0.889805 up to floating-point precision.
The unweighted resolved-topology contrast shifts by -0.873833; callable-site exposure changes the topology-only contrast by -0.004238 beyond that unweighted comparison. Thus the large predicted shift is primarily a change in discordant-topology composition, with a secondary contribution from conditional internal branches; callable exposure is a small additional effect in this summary.

By interval, `(D_site, D_topology, D_tree)` were PNAS_1 `(-0.2584, -0.9261, -0.9709)`, PNAS_2 `(-0.3305, -0.8734, -0.9657)`, and PNAS_3 `(-0.2533, -0.9482, -0.9839)`. The same sign occurs in all three intervals, but the gene-tree predictions remain substantially more extreme than nucleotide D.

The alternative galliform outgroup gives N61 `(D_site, D_topology, D_tree) = (-0.2529, -0.9364, -0.9784)` inside the combined PNAS mask, confirming robustness of the qualitative pattern.

## Interpretation

The local trees recover the sign but substantially overpredict the observed magnitude. Agreement in magnitude is not required by the infinite-sites approximation: observed ABBA/BABA counts also depend on finite-sites homoplasy, ancestral-state errors, alignment filtering, model misspecification, and tree-estimation dependence. The named trees were estimated from the same alignments, so this is not independent validation.

A finite-sites D prediction was not fit here. A JC69 or richer model would require an explicit likelihood over the full quartet branch configuration, ascertainment, and empirical base composition; applying a fitted correction without independent calibration would be misleading. Transversion-only sensitivity is retained in the observed tables and robustness output.

The PNAS intervals are published GalGal4 predicates applied to the Stiller GalGal6 axis; no exact lift was introduced. Results are therefore regional compatibility analyses, not nucleotide-resolution boundary validation.

The analysis is compatible with a shared recombination-partitioned genealogy regime, but it does not discriminate MSRC from introgression, incomplete lineage sorting, mutation-rate variation, or gene-tree error. Both the predicted genealogies and observed D derive from the same underlying sequence data.
