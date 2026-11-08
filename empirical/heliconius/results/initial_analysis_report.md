# Initial Heliconius MSRC pilot report

## Status

This is an exploratory pilot. Dryad’s landing page verifies the four candidate tree tables and `eratoGenomeInfo.zip`; local contents are **available and parsed**. No result is reported for an inversion interval whose assembly-specific numeric boundaries have not been independently validated.

## Sources and design

The targeted source is Edelman et al. (2019), Science, DOI 10.1126/science.aaw2090, Dryad DOI 10.5061/dryad.b7bj832. Primary coordinates are H. erato demophoon v1; Hmel2.5 remains separate. The primary analysis enumerates all four-taxon subsets of each source-specific seven-taxon universe (35 quartets for complete trees). Per-quartet baselines are frozen from source-specific genome-wide backgrounds excluding chromosomes 2 and 15; alternatives are assigned deterministically. The earlier `HeraRef, Hhim, Htel, Hhsa` result is retained only as a legacy focal-quartet comparison.

## Audit and normalization

The pipeline found 57488 usable local-tree rows and 57488 rows with complete terminal branch lengths. A single tree contributes one topology assignment per available quartet; aggregate q1/q2/q3 are fractions across quartets within a tree and are not independent-window probabilities. The 50-kb table is summarized into explicit physical bins with `n_trees`; overlapping sliding windows are not treated as independent replicates.

## Inversions and recombination

The validated H. erato 2b interval and the separate Hmel2.5 2b/15b intervals are used only in their documented assembly coordinate systems. H. erato 15b remains unresolved in the available assembly mapping, so no Hera 15b inside/outside result is claimed. `eratoGenomeInfo.zip` supplies quantitative LDhelmet-derived recombination estimates; these are population-scaled/source-scaled values, not direct crossover measurements in inversion heterozygotes. Structural evidence of an inversion is not equated with measured suppression.

## Findings

The all-quartet results and figures are descriptive spatial summaries. The strongest transitions are genome-wide exploratory bin boundaries and do not independently validate MSRC or inversion effects; the documented 2b comparison does not show a clean alternative-topology replacement.

## Interpretation and limitations

Any topology shift could reflect introgression, ancestral polymorphism/ILS, mutation-rate or alignment variation, or tree-estimation error; the data cannot by themselves establish MSRC as the mechanism. Substitution-scale patristic distances are not coalescence times.

## Missing resources / next analysis

The public local trees and recombination metadata have been downloaded, parsed, analyzed, and plotted. Next steps are block-aware regional comparisons, boundary-sensitivity analyses, and independent structural validation for the unresolved H. erato 15b interval.

## Actual parsed results

The normalized dataset contains 57488 local trees: 32870 10-kb Hera trees, 13743 Hera 50-kb trees, and 10875 Hmel2.5 50-kb trees. Primary Hera coverage includes chromosomes 1–21; chromosome 2 has 1563 trees and chromosome 15 has 1762 trees.

46174 primary H. erato local-tree windows have quantitative LDhelmet-derived values, with mean 2.603355 cM/Mb and range 0.120808–8.671181 cM/Mb. These are population-scaled/source-scaled recombination estimates from the published LDhelmet workflow, not direct crossover measurements in the inversion heterozygotes.

| chromosome | region | n | q1 | q2 | q3 | D | internal branch mean (n) | mean pairwise patristic (n) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| chr2 | inside | 19 | 1.000000 | 0.000000 | 0.000000 | -1.000000 | 0.018260 (19) | 0.124377 (19) |
| chr2 | overlap | 3 | 1.000000 | 0.000000 | 0.000000 | -1.000000 | 0.021227 (3) | 0.109516 (3) |
| chr2 | outside | 1541 | 0.993511 | 0.002596 | 0.003894 | -0.990915 | 0.029865 (1541) | 0.155574 (1541) |
| chr15 | unclassified | 1762 | 0.995460 | 0.001703 | 0.002838 | -0.993757 | 0.034886 (1762) | 0.148071 (1762) |

The chr2 exploratory quartet is overwhelmingly baseline inside (q1=1.000, n=19) and outside (q1=0.9935, n=1541); q2 is 0 inside and 0.0026 outside. This does not support a localized alternative-topology replacement in the tested quartet.

Figures generated: `heliconius_chr2_tracks.*`, `heliconius_chr15_tracks.*`, `heliconius_figure3_chr2b_inside_outside.*`, and `heliconius_figure4_quartet_simplex.*`.

## All-quartet analysis

The generalized analysis enumerates all four-taxon subsets of each seven-taxon source: 35 quartets per local tree when all taxa are present. The Hera 10-kb and Hera 50-kb sources use the same taxon universe but remain separate spatial datasets; Hmel2.5 uses its own seven-taxon universe and assembly. The baseline for each quartet is frozen from its source-specific genome-wide background excluding chromosomes 2 and 15. The two alternatives use a deterministic canonical ordering, with no selection based on focal-region enrichment. Mapping details are in `data/heliconius/processed/all_quartet_reference_baselines.tsv`.

The aggregate q1/q2/q3 values below are the fractions of resolved quartets within each individual local tree, then spatially averaged; they are not fixed-quartet frequencies across independent windows. Individual assignments are retained in compressed `data/heliconius/processed/all_quartet_assignments.tsv.gz`, and quartet-specific region summaries are in `data/heliconius/processed/all_quartet_quartet_summary.tsv`.

| source | local trees | quartets per tree | mean q1 | mean q2 | mean q3 | mean D |
|---|---:|---:|---:|---:|---:|---:|
| erato_10kb | 32870 | 35 | 0.8173 | 0.1141 | 0.0686 | -0.7032 |
| erato_50kb_hera | 13743 | 35 | 0.8628 | 0.0995 | 0.0377 | -0.7634 |
| erato_50kb_hmel | 10875 | 35 | 0.8906 | 0.0898 | 0.0196 | -0.8008 |

The strongest genome-wide adjacent 50-kb-bin transitions in the deterministic summary were:

| source | chromosome | boundary between bins (bp) | D jump |
|---|---|---:|---:|
| erato_50kb_hmel | chr11 | 899999–900000 | -1.286 |
| erato_50kb_hmel | chr11 | 799999–800000 | 1.271 |
| erato_50kb_hmel | chr15 | 5249999–5250000 | -1.171 |
| erato_50kb_hmel | chr13 | 12049999–12050000 | -1.029 |
| erato_50kb_hmel | chr13 | 11949999–11950000 | 1.029 |

These are exploratory spatial transitions, not independently validated inversion boundaries. In particular, the 2b interval is independently documented only for H. erato chromosome 2 and Hmel2.5 scaffold coordinates as listed above; H. erato 15b remains unresolved in the available assembly mapping. The all-quartet chromosome-2 summaries show broad baseline/alternative competition rather than a clean, inversion-specific alternative-topology replacement, and no MSRC mechanism is claimed.

All-quartet figures generated include `heliconius_genomewide_all_quartet_tracks.*`, source-specific chromosome-2 and chromosome-15 tracks for 10-kb Hera, 50-kb Hera, and 50-kb Hmel2.5, plus the earlier focal-quartet figures.
