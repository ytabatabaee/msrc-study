# Initial Heliconius MSRC pilot report

## Status

This is an exploratory pilot. Dryad’s landing page verifies the four candidate tree tables and `eratoGenomeInfo.zip`; local contents are **available and parsed**. No result is reported for an inversion interval whose assembly-specific numeric boundaries have not been independently validated.

## Sources and design

The targeted source is Edelman et al. (2019), Science, DOI 10.1126/science.aaw2090, Dryad DOI 10.5061/dryad.b7bj832. Primary coordinates are H. erato demophoon v1; Hmel2.5 remains separate. The frozen exploratory quartet is `HeraRef, Hhim, Htel, Hhsa`. Baseline is `HeraRef,Hhim|Hhsa,Htel` and alternative is `HeraRef,Htel|Hhim,Hhsa`; this is topology-based because arrangement assignments are not established for the four tips.

## Audit and normalization

The pipeline found 57488 usable local-tree rows and 57488 rows with complete terminal branch lengths. A single tree contributes a 0/1 topology indicator. The 50-kb table is summarized into explicit physical bins with `n_trees`; overlapping sliding windows are not treated as independent replicates.

## Inversions and recombination

The validated H. erato 2b interval is 9,114,112–9,292,966 (1-based inclusive), converted from the published Herato0206 scaffold interval using the archived scaffold-to-chromosome transitions. The published Hmel2.5 2b and 15b scaffold intervals are retained separately. H. erato 15b remains unresolved because the independent coordinate table reports n/a for that assembly. The parsed `eratoGenomeInfo.zip` supplies quantitative LDhelmet-derived recombination values; these are population-scaled estimates, not direct crossover measurements in inversion heterozygotes. Structural evidence of an inversion is not equated with a measured recombination rate.

## Findings

The descriptive tables and figures should be interpreted as topology/distance summaries, not independent-window significance tests. The chr2 2b inside/outside comparison is descriptive and uses the validated interval; chr15 15b inside/outside comparison remains blocked for the Hera tracks because its Hera coordinate is unresolved.

## Interpretation and limitations

Any topology shift could reflect introgression, ancestral polymorphism/ILS, mutation-rate or alignment variation, or tree-estimation error; the data cannot by themselves establish MSRC as the mechanism. Substitution-scale patristic distances are not coalescence times.

## Missing resources / next analysis

Download the small Dryad files, inspect the audit, validate inversion boundaries against a primary structural source or explicit assembly mapping, parse recombination metadata, then rerun normalization, statistics, and figures. Only after that should block-aware inside/outside comparisons and boundary sensitivity analyses be added.

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
