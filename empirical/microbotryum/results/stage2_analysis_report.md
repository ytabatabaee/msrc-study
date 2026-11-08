# Microbotryum Stage 2 analysis report

## Recommendation

**CONDITIONAL GO — strong but taxon-pair-specific stratum-level evidence, with spatial validation incomplete.** Stage 2 recovered no defensible chromosome-scale coordinates or orthogroup-to-tree mapping. It therefore does not claim a spatial change point or generate chromosome tracks.

## Assembly and source validation

The paper identifies the phased biological assemblies as M. superbum strain 1065 and M. shykoffianum strain c212. The local `MPbarcode50.fa` is described only as a haploid Flye ONT assembly with anonymous contigs, and its GTF uses an `Msup_tig...` namespace. The paper separately describes a submitted M. superbum strain 6P annotation. The tree tips are Mv-sup-1065 and Mv-cat-C212. Because no local file links MPbarcode50 to 1065, 6P, or a mating type, the GTF match is provisional.

The source paper assigns strain C212 to M. shykoffianum based on barcode genes, despite the deposited tree prefix `Mv-cat`. Stage 2 therefore labels the C212 taxon as M. shykoffianum (C212) and retains `Mv-cat-C212` as the literal sample identifier. The `Mv-cat` prefix is not treated as an independent species assignment.

The authors' paper documents a PR chromosome fused with part of the ancestral HD chromosome, small PARs at chromosome edges, and M. lagerheimii as a proxy for ancestral gene order. It also describes black as shared/older and turquoise/ruby as young species-specific strata. These are independent biological annotations, but the submitted files currently available here do not provide a chromosome-coordinate transform.

## Complete tree-to-coordinate audit

The one-row-per-tree mapping contains 260 records: 260 provisional contig-local matches, 0 unmapped records, 0 chromosome-scale coordinates, 0 recovered orthogroup IDs, and 260 deposited stratum labels. The numeric suffix match to MPbarcode50 is not treated as proof of biological orthology.

The six taxa are retained as Mv-cat-C212, Mv-lag-1253, and Mv-sup-1065, each with A1/A2 labels. The published biological assignment used in this report is C212 = M. shykoffianum, 1065 = M. superbum, and the `Mv-lag` sample as M. lagerheimii; sample prefixes remain preserved in all machine-readable tables.

## Cross-species-pair results

The three complete species-pair quartets were analyzed from the existing Stage 1 assignments. q1 is species clustering and q2 is mating-type clustering, with fixed labels across strata.

| stratum | species pair | n | q1 | q2 | q3 | D | B | internal |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| HD_RR | cat|lag | 9 | 1.000 | 0.000 | 0.000 | -1.000 | 0.13121 | 0.19582 |
| HD_RR | cat|sup | 9 | 1.000 | 0.000 | 0.000 | -1.000 | 0.03056 | 0.04360 |
| HD_RR | lag|sup | 9 | 1.000 | 0.000 | 0.000 | -1.000 | 0.12157 | 0.17920 |
| PARs | cat|lag | 29 | 1.000 | 0.000 | 0.000 | -1.000 | 0.24247 | 0.36311 |
| PARs | cat|sup | 29 | 1.000 | 0.000 | 0.000 | -1.000 | 0.07030 | 0.10437 |
| PARs | lag|sup | 29 | 1.000 | 0.000 | 0.000 | -1.000 | 0.22953 | 0.34267 |
| black | cat|lag | 134 | 1.000 | 0.000 | 0.000 | -1.000 | 0.19286 | 0.13676 |
| black | cat|sup | 134 | 0.119 | 0.873 | 0.007 | 0.754 | 0.14770 | 0.05926 |
| black | lag|sup | 134 | 1.000 | 0.000 | 0.000 | -1.000 | 0.15520 | 0.13411 |
| blue | cat|lag | 17 | 1.000 | 0.000 | 0.000 | -1.000 | 0.40828 | 0.31990 |
| blue | cat|sup | 17 | 1.000 | 0.000 | 0.000 | -1.000 | 0.28713 | 0.15573 |
| blue | lag|sup | 17 | 1.000 | 0.000 | 0.000 | -1.000 | 0.16459 | 0.21025 |
| dark_ruby | cat|lag | 25 | 1.000 | 0.000 | 0.000 | -1.000 | 0.34134 | 0.47030 |
| dark_ruby | cat|sup | 25 | 0.600 | 0.400 | 0.000 | -0.200 | 0.25593 | 0.31392 |
| dark_ruby | lag|sup | 25 | 1.000 | 0.000 | 0.000 | -1.000 | 0.15375 | 0.17560 |
| in_HD | cat|lag | 2 | 1.000 | 0.000 | 0.000 | -1.000 | 0.13806 | 0.12426 |
| in_HD | cat|sup | 2 | 0.000 | 1.000 | 0.000 | 1.000 | 0.08562 | 0.04156 |
| in_HD | lag|sup | 2 | 1.000 | 0.000 | 0.000 | -1.000 | 0.12947 | 0.12426 |
| light_ruby | cat|lag | 13 | 1.000 | 0.000 | 0.000 | -1.000 | 0.42767 | 0.63820 |
| light_ruby | cat|sup | 13 | 1.000 | 0.000 | 0.000 | -1.000 | 0.26830 | 0.37704 |
| light_ruby | lag|sup | 13 | 1.000 | 0.000 | 0.000 | -1.000 | 0.24329 | 0.33726 |
| orange | cat|lag | 4 | 1.000 | 0.000 | 0.000 | -1.000 | 0.14045 | 0.11147 |
| orange | cat|sup | 4 | 0.250 | 0.750 | 0.000 | 0.500 | 0.10283 | 0.04266 |
| orange | lag|sup | 4 | 1.000 | 0.000 | 0.000 | -1.000 | 0.12327 | 0.10810 |
| purple | cat|lag | 3 | 1.000 | 0.000 | 0.000 | -1.000 | 0.17516 | 0.04899 |
| purple | cat|sup | 3 | 0.000 | 1.000 | 0.000 | 1.000 | 0.15853 | 0.09364 |
| purple | lag|sup | 3 | 1.000 | 0.000 | 0.000 | -1.000 | 0.11899 | 0.04899 |
| red | cat|lag | 2 | 1.000 | 0.000 | 0.000 | -1.000 | 2.71194 | 2.31042 |
| red | cat|sup | 2 | 1.000 | 0.000 | 0.000 | -1.000 | 1.77099 | 0.84642 |
| red | lag|sup | 2 | 1.000 | 0.000 | 0.000 | -1.000 | 3.09172 | 1.08613 |
| turquoise | cat|lag | 22 | 1.000 | 0.000 | 0.000 | -1.000 | 0.13484 | 0.17541 |
| turquoise | cat|sup | 22 | 0.955 | 0.045 | 0.000 | -0.909 | 0.04997 | 0.02626 |
| turquoise | lag|sup | 22 | 1.000 | 0.000 | 0.000 | -1.000 | 0.12557 | 0.16461 |

For the primary C212–1065 pair, q2 changes from 0.000 in PAR to 0.873 in black and 0.045 in turquoise. This reversal is **not shared** by the other two pairings: C212–lagerheimii and lagerheimii–1065 remain q1 = 1.000 in black. The black signal is therefore taxon-pair-specific in these trees, not evidence of a uniform three-species reversal.

This restriction weakens a broad ancestral-process interpretation and makes introgression, species-specific history, or taxon-specific gene-tree effects important alternatives. It remains a stratum contrast rather than a boundary-localized test.

## Boundary and age analysis

PAR–black and black–turquoise effects are reported in `stage2_boundary_transition_summary.tsv`, but every row is marked `not_spatially_testable`: no validated physical or gene-order boundary was recovered. These deltas must not be interpreted as topology-derived boundary discovery.

The age pattern is compatible with the published chronology: black is older/shared and turquoise is younger/species-specific; the paper reports turquoise and ruby as lower-dS young strata relative to black and identifies light/dark ruby sub-strata. The gene-tree data show the strongest mating-type support in black only for the C212–1065 pair, not a monotonic or universal age trend. Small strata such as orange, purple, and red are too sparsely sampled for strong age-trend claims.

## Spatial analysis status

No q1/q2/q3, D, A, or B chromosome tracks were generated. The available GTF matches span multiple anonymous contigs and the tree prefix cannot be linked to an orthogroup or published ancestral gene-order rank. Contig-local plotting would not establish adjacency of strata, so the spatial analysis remains blocked.

## Biological interpretation and limitations

The C212–1065 black-stratum reversal is compatible with an older mating-type-associated genealogy maintained under recombination suppression, but the lack of the same reversal in the other two species-pair quartets is a major limitation. Selection, gene conversion, introgression, branch-length heterogeneity, alignment quality, and gene-tree estimation error remain alternatives. In particular, the fungal mating-type chromosome history includes chromosome fusion, rearrangement, selfing/intra-tetrad mating, and loss of HD mating function, which differs from an idealized two-background model.

## Files and next blocker

- `stage2_coordinate_mapping.tsv`: complete one-row-per-tree mapping audit.
- `stage2_stratum_order_metadata.tsv`: published regime and age annotations without invented coordinates.
- `cross_species_pair_summary.tsv` and `cross_species_pair_bootstrap.tsv`: all three complete species-pair quartet comparisons.
- `stage2_boundary_transition_summary.tsv`: descriptive stratum contrasts explicitly marked non-spatial.
- `microbotryum_stage2_cross_species_quartets.*`: cross-species quartet figure.
- `microbotryum_stage2_genomic_organization_status.*`: evidence-status organization diagram, not a physical map.

The concrete next step is to obtain the authors' source-data/table containing the M. lagerheimii ancestral-order ranks and the exact tree-prefix-to-orthogroup conversion. With those files, the existing assignments can be projected onto gene-order bins without reconstructing trees.

References: [Lucotte et al. 2025](https://doi.org/10.1038/s41467-025-60222-5), [Dryad](https://doi.org/10.5061/dryad.0k6djhb89), [authors' repository](https://github.com/eliselucotte/MicrobotryumOnDianthus).
