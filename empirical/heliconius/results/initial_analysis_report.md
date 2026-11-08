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

32712 primary H. erato local-tree windows have quantitative LDhelmet-derived values, with mean 2.633240 cM/Mb and range 0.120808–7.751381 cM/Mb. These are population-scaled/source-scaled recombination estimates from the published LDhelmet workflow, not direct crossover measurements in the inversion heterozygotes.

| chromosome | region | n | q1 | q2 | q3 | D | internal branch mean (n) | mean pairwise patristic (n) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| chr2 | inside | 19 | 1.000000 | 0.000000 | 0.000000 | -1.000000 | 0.018260 (19) | 0.124377 (19) |
| chr2 | overlap | 3 | 1.000000 | 0.000000 | 0.000000 | -1.000000 | 0.021227 (3) | 0.109516 (3) |
| chr2 | outside | 1541 | 0.993511 | 0.002596 | 0.003894 | -0.990915 | 0.029865 (1541) | 0.155574 (1541) |
| chr15 | unclassified | 1762 | 0.995460 | 0.001703 | 0.002838 | -0.993757 | 0.034886 (1762) | 0.148071 (1762) |

The chr2 exploratory quartet is overwhelmingly baseline inside (q1=1.000, n=19) and outside (q1=0.9935, n=1541); q2 is 0 inside and 0.0026 outside. This does not support a localized alternative-topology replacement in the tested quartet.

Figures generated: `heliconius_chr2_tracks.*`, `heliconius_chr15_tracks.*`, `heliconius_figure3_chr2b_inside_outside.*`, and `heliconius_figure4_quartet_simplex.*`.

## All-quartet analysis

The generalized analysis enumerates all four-taxon subsets of each seven-taxon source: 35 quartets per local tree when all taxa are present. The Hera 10-kb and Hera 50-kb sources use the same taxon universe but remain separate spatial datasets; Hmel2.5 uses its own seven-taxon universe and assembly. The baseline for each quartet is the fixed published species-history resolution documented in `config/published_species_trees.tsv`; the source-specific genome-wide background is retained as an agreement comparison. The two alternatives use a deterministic canonical ordering, with no selection based on focal-region enrichment. Mapping details are in `data/heliconius/processed/all_quartet_reference_baselines.tsv`.

The aggregate q1/q2/q3 values below are the fractions of resolved quartets within each individual local tree, then spatially averaged; they are not fixed-quartet frequencies across independent windows. Individual assignments are retained in compressed `data/heliconius/processed/all_quartet_assignments.tsv.gz`, and quartet-specific region summaries are in `data/heliconius/processed/all_quartet_quartet_summary.tsv`.

| source | local trees | quartets per tree | mean q1 | mean q2 | mean q3 | mean D |
|---|---:|---:|---:|---:|---:|---:|
| erato_10kb | 32870 | 35 | 0.6784 | 0.1486 | 0.1730 | -0.5297 |
| erato_50kb_hera | 13743 | 35 | 0.6736 | 0.1363 | 0.1900 | -0.5373 |
| erato_50kb_hmel | 10875 | 35 | 0.6453 | 0.2542 | 0.1006 | -0.3911 |

The strongest genome-wide adjacent 50-kb-bin transitions in the deterministic summary were:

| source | chromosome | boundary between bins (bp) | D jump |
|---|---|---:|---:|
| erato_10kb | chr3 | 949999–1100000 | 0.886 |
| erato_10kb | chr17 | 9049999–9050000 | 0.880 |
| erato_10kb | chr15 | 15099999–15100000 | 0.829 |
| erato_10kb | chr3 | 1149999–1150000 | -0.762 |
| erato_50kb_hmel | chr15 | 5199999–5200000 | 0.743 |

These are exploratory spatial transitions, not independently validated inversion boundaries. The source-matched 2b definitions and the unresolved H. erato 15b mapping are detailed in the critical reassessment section. The all-quartet chromosome-2 summaries show broad baseline/alternative competition rather than a clean, inversion-specific alternative-topology replacement, and no MSRC mechanism is claimed.

All-quartet figures generated include `heliconius_genomewide_all_quartet_tracks.*`, source-specific chromosome-2 and chromosome-15 tracks for 10-kb Hera, 50-kb Hera, and Hmel2.5, plus the critical quartet heatmaps and source-matched 2b zoom figures.

## Critical reassessment (published topology, coordinate-matched regions, sustained regimes)

This report supplements `initial_analysis_report.md`. It uses the downloaded Edelman et al. local trees, published species-history resolutions, and source-specific physical coordinates. All run/segment summaries use non-overlapping 50-kb bins; sliding-window trees remain descriptive linked observations.

## Coordinate resolution

Thawornwattana et al. Supplementary Table 2 contains two distinct 2b entries. The `etales-8spp` entry pairs Hmel202001o:2,147,755–2,230,606 with Herato0206:880,283–1,059,137, while the `hmelv25-res` entry uses Hmel202001o:4,398,947–6,357,046 and has no H. erato coordinate. They are therefore retained as different dataset-specific definitions, not equated by label. The available Hmel-coordinate erato-clade trees cover both positions, but their source is not the `hmelv25-res` dataset; only the etales-8spp interval is used for their region classification. H. erato 15b remains unresolved in the independent table.

## Published topology comparison

Per-quartet q1 now uses the fixed published species-history topology from Edelman et al. Fig. 1A/2 Tree 2, represented in `config/published_species_trees.tsv`; q2 and q3 are the other two deterministic resolutions. The prior background-majority topology is retained in the baseline mapping with an agreement flag. This makes disagreements visible rather than allowing region-specific relabeling.

## Aggregate dominance

| source | chromosome | n windows | q2>q1 | q3>q1 | q2+q3>q1 | longest run (bp) |
|---|---:|---:|---:|---:|---:|---:|
| erato_10kb | 1 | 2051 | 0.002 | 0.001 | 0.012 | 300000 |
| erato_10kb | 10 | 2300 | 0.003 | 0.003 | 0.017 | 250000 |
| erato_10kb | 11 | 1359 | 0.000 | 0.006 | 0.023 | 250000 |
| erato_10kb | 12 | 1993 | 0.002 | 0.003 | 0.012 | 150000 |
| erato_10kb | 13 | 2147 | 0.000 | 0.002 | 0.010 | 250000 |
| erato_10kb | 14 | 1114 | 0.004 | 0.004 | 0.036 | 200000 |
| erato_10kb | 15 | 1251 | 0.004 | 0.006 | 0.031 | 700000 |
| erato_10kb | 16 | 1245 | 0.005 | 0.005 | 0.033 | 150000 |
| erato_10kb | 17 | 1755 | 0.002 | 0.006 | 0.025 | 250000 |
| erato_10kb | 18 | 2021 | 0.002 | 0.002 | 0.018 | 300000 |
| erato_10kb | 19 | 2070 | 0.002 | 0.004 | 0.018 | 150000 |
| erato_10kb | 2 | 1108 | 0.005 | 0.008 | 0.041 | 1100000 |
| erato_10kb | 20 | 1759 | 0.002 | 0.005 | 0.018 | 250000 |
| erato_10kb | 21 | 1577 | 0.004 | 0.001 | 0.008 | 200000 |
| erato_10kb | 3 | 1248 | 0.003 | 0.006 | 0.025 | 250000 |
| erato_10kb | 4 | 1127 | 0.000 | 0.005 | 0.030 | 200000 |
| erato_10kb | 5 | 1176 | 0.002 | 0.002 | 0.030 | 150000 |
| erato_10kb | 6 | 1638 | 0.002 | 0.001 | 0.012 | 200000 |
| erato_10kb | 7 | 1762 | 0.001 | 0.005 | 0.015 | 250000 |
| erato_10kb | 8 | 1115 | 0.002 | 0.006 | 0.027 | 200000 |
| erato_10kb | 9 | 1054 | 0.007 | 0.003 | 0.031 | 150000 |
| erato_50kb_hera | 1 | 854 | 0.000 | 0.000 | 0.000 | 400000 |
| erato_50kb_hera | 10 | 940 | 0.001 | 0.002 | 0.004 | 500000 |
| erato_50kb_hera | 11 | 559 | 0.000 | 0.000 | 0.000 | 350000 |
| erato_50kb_hera | 12 | 810 | 0.000 | 0.000 | 0.000 | 400000 |
| erato_50kb_hera | 13 | 900 | 0.000 | 0.000 | 0.000 | 750000 |
| erato_50kb_hera | 14 | 481 | 0.000 | 0.000 | 0.006 | 300000 |
| erato_50kb_hera | 15 | 511 | 0.000 | 0.000 | 0.000 | 750000 |
| erato_50kb_hera | 16 | 523 | 0.000 | 0.002 | 0.002 | 550000 |
| erato_50kb_hera | 17 | 723 | 0.000 | 0.000 | 0.000 | 650000 |
| erato_50kb_hera | 18 | 824 | 0.001 | 0.000 | 0.005 | 550000 |
| erato_50kb_hera | 19 | 866 | 0.000 | 0.000 | 0.000 | 400000 |
| erato_50kb_hera | 2 | 455 | 0.000 | 0.000 | 0.000 | 1200000 |
| erato_50kb_hera | 20 | 746 | 0.000 | 0.000 | 0.000 | 600000 |
| erato_50kb_hera | 21 | 644 | 0.000 | 0.000 | 0.000 | 400000 |
| erato_50kb_hera | 3 | 534 | 0.000 | 0.000 | 0.000 | 850000 |
| erato_50kb_hera | 4 | 483 | 0.000 | 0.000 | 0.000 | 350000 |
| erato_50kb_hera | 5 | 502 | 0.000 | 0.000 | 0.000 | 200000 |
| erato_50kb_hera | 6 | 702 | 0.000 | 0.000 | 0.000 | 350000 |
| erato_50kb_hera | 7 | 765 | 0.000 | 0.000 | 0.000 | 500000 |
| erato_50kb_hera | 8 | 463 | 0.000 | 0.000 | 0.000 | 400000 |
| erato_50kb_hera | 9 | 458 | 0.000 | 0.000 | 0.000 | 300000 |
| erato_50kb_hmel | 0 | 100 | 0.000 | 0.000 | 0.010 | 100000 |
| erato_50kb_hmel | 1 | 684 | 0.000 | 0.000 | 0.003 | 250000 |
| erato_50kb_hmel | 10 | 717 | 0.001 | 0.001 | 0.001 | 250000 |
| erato_50kb_hmel | 11 | 467 | 0.009 | 0.004 | 0.024 | 500000 |
| erato_50kb_hmel | 12 | 647 | 0.000 | 0.000 | 0.008 | 400000 |
| erato_50kb_hmel | 13 | 719 | 0.000 | 0.003 | 0.008 | 550000 |
| erato_50kb_hmel | 14 | 362 | 0.000 | 0.003 | 0.017 | 150000 |
| erato_50kb_hmel | 15 | 402 | 0.005 | 0.005 | 0.007 | 500000 |
| erato_50kb_hmel | 16 | 396 | 0.000 | 0.000 | 0.008 | 350000 |
| erato_50kb_hmel | 17 | 585 | 0.000 | 0.003 | 0.005 | 600000 |
| erato_50kb_hmel | 18 | 669 | 0.000 | 0.000 | 0.006 | 400000 |
| erato_50kb_hmel | 19 | 651 | 0.000 | 0.000 | 0.006 | 250000 |
| erato_50kb_hmel | 2 | 352 | 0.000 | 0.000 | 0.011 | 2000000 |
| erato_50kb_hmel | 20 | 589 | 0.003 | 0.003 | 0.014 | 1050000 |
| erato_50kb_hmel | 21 | 532 | 0.000 | 0.000 | 0.000 | 200000 |
| erato_50kb_hmel | 3 | 414 | 0.002 | 0.005 | 0.007 |  |
| erato_50kb_hmel | 4 | 375 | 0.011 | 0.011 | 0.027 | 250000 |
| erato_50kb_hmel | 5 | 385 | 0.000 | 0.003 | 0.003 | 200000 |
| erato_50kb_hmel | 6 | 550 | 0.004 | 0.002 | 0.007 | 450000 |
| erato_50kb_hmel | 7 | 570 | 0.000 | 0.000 | 0.007 | 200000 |
| erato_50kb_hmel | 8 | 365 | 0.005 | 0.000 | 0.008 | 250000 |
| erato_50kb_hmel | 9 | 344 | 0.003 | 0.003 | 0.009 | 350000 |

The complete per-quartet dominance table is `all_quartet_dominance_summary.tsv`; aggregate dominance is not required for a localized MSRC-compatible signal because an effect can involve only a subset of quartets.

## Sustained regimes

| source | chromosome | start | end | span | bins | trees | mean alternative-minus-baseline | direction |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| erato_50kb_hmel | 2 | 4350000 | 6349999 | 2000000 | 40 | 79 | -0.376 | baseline |
| erato_50kb_hera | 2 | 5600000 | 6799999 | 1200000 | 24 | 47 | -0.393 | baseline |
| erato_10kb | 2 | 3750000 | 4849999 | 1100000 | 22 | 101 | -0.315 | baseline |
| erato_10kb | 2 | 5750000 | 6799999 | 1050000 | 21 | 102 | -0.353 | baseline |
| erato_50kb_hera | 2 | 3800000 | 4849999 | 1050000 | 21 | 40 | -0.400 | baseline |
| erato_50kb_hmel | 20 | 7000000 | 8049999 | 1050000 | 21 | 42 | 0.114 | alternative |
| erato_50kb_hera | 3 | 10500000 | 11349999 | 850000 | 17 | 34 | 0.114 | alternative |
| erato_50kb_hera | 13 | 4000000 | 4749999 | 750000 | 15 | 30 | 0.114 | alternative |
| erato_50kb_hera | 15 | 2250000 | 2999999 | 750000 | 15 | 30 | -0.328 | baseline |
| erato_10kb | 15 | 2300000 | 2999999 | 700000 | 14 | 70 | -0.310 | baseline |

The strongest runs are not automatically inversion effects; they are genome-wide candidates defined before looking at the annotated intervals. Missing physical bins are not bridged.

## Quartet-specific chr2 effects

The complete table is `quartet_region_effects_chr2.tsv`. It contains every quartet with both inside and outside observations for each source; the ranking below identifies contributors by absolute change in alternative support without discarding the remaining quartets.

| source | quartet | n inside | n outside | delta alternative | delta internal | delta B |
|---|---|---:|---:|---:|---:|---:|
| erato_50kb_hera | Hdem|HeraRef|Hhsa|Htel | 6 | 447 | -0.591 | -0.0088 | -0.0219 |
| erato_50kb_hera | Hdem|Hhim|Hhsa|Htel | 6 | 447 | -0.591 | -0.0088 | -0.0219 |
| erato_50kb_hera | HeraRef|Hhsa|Hsar|Htel | 6 | 447 | -0.591 | -0.0088 | -0.0219 |
| erato_50kb_hera | Hhim|Hhsa|Hsar|Htel | 6 | 447 | -0.591 | -0.0088 | -0.0219 |
| erato_50kb_hera | HeraRef|Hhsa|HmelRef|Htel | 6 | 447 | 0.564 | -0.0041 | -0.0219 |
| erato_50kb_hera | Hhim|Hhsa|HmelRef|Htel | 6 | 447 | 0.564 | -0.0041 | -0.0219 |
| erato_50kb_hmel | Etal|HeraRef|Hhsa|Htel | 2 | 346 | -0.500 | -0.0100 | -0.0174 |
| erato_50kb_hmel | Etal|Hhim|Hhsa|Htel | 2 | 346 | -0.500 | -0.0100 | -0.0174 |
| erato_50kb_hmel | Hdem|HeraRef|Hhsa|Htel | 2 | 346 | 0.451 | -0.0064 | -0.0174 |
| erato_50kb_hmel | Hdem|Hhim|Hhsa|Htel | 2 | 346 | 0.451 | -0.0064 | -0.0174 |
| erato_50kb_hmel | HeraRef|Hhsa|Hsar|Htel | 2 | 346 | 0.451 | -0.0064 | -0.0174 |
| erato_50kb_hmel | Hhim|Hhsa|Hsar|Htel | 2 | 346 | 0.451 | -0.0064 | -0.0174 |

## Branch-length-only candidates

| source | chromosome | bin | B robust z | topology effect | candidate |
|---|---:|---:|---:|---:|---|
| erato_50kb_hmel | 4 | 3900000–3949999 | 713.52 | 0.086 | true |
| erato_50kb_hmel | 4 | 4100000–4149999 | 56.73 | -0.086 | true |
| erato_50kb_hmel | 3 | 1400000–1449999 | 30.40 | 0.029 | true |
| erato_10kb | 4 | 6200000–6249999 | 25.49 | 0.029 | true |
| erato_50kb_hmel | 19 | 14150000–14199999 | 20.59 | 0.029 | true |
| erato_50kb_hmel | 2 | 2700000–2749999 | 15.08 | -0.086 | true |
| erato_10kb | 16 | 15250000–15299999 | 14.60 | -0.086 | true |
| erato_10kb | 6 | 17850000–17899999 | 14.33 | 0.086 | true |
| erato_10kb | 15 | 4350000–4399999 | 13.61 | -0.086 | true |
| erato_50kb_hmel | 12 | 6650000–6699999 | 11.23 | 0.029 | true |

These are descriptive robust deviations from the chromosome-specific median; patristic distances remain source-tree substitution units, not coalescence times. A relative-rate or alignment-quality analysis would be needed to distinguish mutation-rate effects from genealogical-time effects.

## Recombination and validated regions

The H. erato LDhelmet track is compared continuously by chromosome-specific recombination quantiles in `recombination_quantile_summary.tsv`. The current data do not provide a direct crossover measurement in inversion heterozygotes. `validated_region_comparison.tsv` contains topology, B, internal branch, and recombination summaries for the source-matched chr2 intervals. No H. erato 15b comparison is reported because the independent mapping is n/a.

## Decision

Recommendation: retain Heliconius as exploratory and do not add it to the MSRC manuscript empirical section. The dataset contains extensive genome-wide discordance and candidate spatial regimes, but the available coordinate-matched comparisons do not yet demonstrate a coherent, independently supported topology-plus-branch-length regime coincident with a verified suppressed-recombination interval. Introgression, ancestral polymorphism/ILS, mutation-rate variation, and tree-estimation error remain strong alternatives.

Figures: `heliconius_*_quartet_heatmap.*` and `heliconius_genomewide_sustained_regimes.*`.
