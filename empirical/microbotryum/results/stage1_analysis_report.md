# Microbotryum Stage 1 empirical analysis

## Recommendation

**CONDITIONAL GO — strong stratum-level signal, but incomplete spatial and assembly validation.** The deposited trees provide a compelling descriptive transition between the PAR control, the older black stratum, and the younger turquoise stratum. A chromosome-coordinate MSRC analysis is not yet justified because the MPbarcode50 annotation-to-tree sample/assembly correspondence and chromosome-scale contig order remain unvalidated.

## Data and coordinate status

The analysis uses the local Dryad `gene_trees` file and 260 parseable records. Each record contains six taxa and complete branch lengths. The 260 trees span 11 deposited stratum labels. The focal quartet is `Mv-cat-C212-A1/A2` plus `Mv-sup-1065-A1/A2`; the third species pair is `Mv-lag-1253-A1/A2`.

The separate `coordinate_validation_report.md` is the authoritative mapping assessment. All 260 tree gene IDs map numerically to `MPbarcode50_tsebra.gtf`, but the matching ID namespace does not establish the same sample or assembly. Contig-local coordinates are retained; contigs are not concatenated into an invented chromosome. Consequently no spatial track is generated.

## Primary quartet

Topology labels were fixed independently of strata: q1 is species clustering, q2 is same-mating-type clustering across species, and q3 is the third resolution. Frequencies are across genes within a stratum, not probabilities from one tree.

| stratum | n | q1 | q2 | q3 | B | A (opposite - same mating type) | internal |
|---|---:|---:|---:|---:|---:|---:|---:|
| HD_RR | 9 | 1.000 [1.000,1.000] | 0.000 [0.000,0.000] | 0.000 [0.000,0.000] | 0.03056 | 0.00000 | 0.04360 |
| PARs | 29 | 1.000 [1.000,1.000] | 0.000 [0.000,0.000] | 0.000 [0.000,0.000] | 0.07030 | 0.00000 | 0.10437 |
| black | 134 | 0.119 [0.067,0.179] | 0.873 [0.813,0.933] | 0.007 [0.000,0.022] | 0.14770 | 0.05470 | 0.05926 |
| blue | 17 | 1.000 [1.000,1.000] | 0.000 [0.000,0.000] | 0.000 [0.000,0.000] | 0.28713 | 0.00000 | 0.15573 |
| dark_ruby | 25 | 0.600 [0.400,0.800] | 0.400 [0.200,0.600] | 0.000 [0.000,0.000] | 0.25593 | 0.01131 | 0.31392 |
| in_HD | 2 | 0.000 [0.000,0.000] | 1.000 [1.000,1.000] | 0.000 [0.000,0.000] | 0.08562 | 0.04156 | 0.04156 |
| light_ruby | 13 | 1.000 [1.000,1.000] | 0.000 [0.000,0.000] | 0.000 [0.000,0.000] | 0.26830 | 0.00000 | 0.37704 |
| orange | 4 | 0.250 [0.000,0.750] | 0.750 [0.250,1.000] | 0.000 [0.000,0.000] | 0.10283 | 0.03672 | 0.04266 |
| purple | 3 | 0.000 [0.000,0.000] | 1.000 [1.000,1.000] | 0.000 [0.000,0.000] | 0.15853 | 0.09364 | 0.09364 |
| red | 2 | 1.000 [1.000,1.000] | 0.000 [0.000,0.000] | 0.000 [0.000,0.000] | 1.77099 | 0.00000 | 0.84642 |
| turquoise | 22 | 0.955 [0.864,1.000] | 0.045 [0.000,0.136] | 0.000 [0.000,0.000] | 0.04997 | 0.00127 | 0.02626 |

The black stratum has q2 = 0.873 (117/134) and q1 = 0.119; the PAR control has q1 = 1.000 (29/29) and q2 = 0; turquoise has q1 = 0.955 and q2 = 0.045. The black-to-turquoise contrast is therefore not a generic effect of every nonrecombining label: the younger turquoise stratum retains mostly species clustering.

The corrected A statistic is the mean distance for opposite mating-type cross-species pairs minus the mean distance for same-mating-type cross-species pairs. B and all branch lengths are in the source tree's substitution-distance units, not calibrated coalescence times. The gene-level bootstrap intervals are descriptive; linked genes do not provide independent physical replicates.

## All 15 quartets

All 15 possible four-taxon subsets are enumerated for the complete six-taxon trees. 3 quartet identities (the three species-pair comparisons) have a biologically meaningful species-versus-mating-type labeling; the remaining 12 quartets retain fixed generic unrooted resolution labels and are not assigned an arrangement interpretation. The full assignment and summary tables preserve every quartet and every gene.

The all-quartet heatmap uses the fixed second resolution for each quartet, not a post hoc choice of the largest local signal. The focal cat–sup quartet is the strongest and most directly interpretable contrast; the other species-pair quartets should be assessed as replication only after their species and mating-type metadata are independently confirmed.

## Recombination regimes and alternatives

The PARs are a categorical recombining control, not a quantitative crossover-rate track. The paper supports an older ancestral nonrecombining black region with trans-specific polymorphism and a younger turquoise species-specific stratum; it also reports additional colored strata, but this audit does not force an age ordering onto those labels without recovered physical boundaries. The observed black signal is compatible with long-lived mating-type-associated genealogies, but mutation-rate differences, alignment quality, gene-tree estimation error, residual gene conversion, selection, and introgression remain alternative explanations. The data do not make MSRC a unique causal explanation.

## Spatial feasibility and blockers

No chromosome-wide tracks, gene-order windows, or boundary tests are reported. A valid next step is to recover the authors' ortholog/contig placement table or an assembly-to-chromosome map, then repeat the analysis with block-aware windows and matched local controls. Until then, stratum-level comparisons are the defensible unit.

## Files

- `all_quartet_assignments.tsv` and `all_quartet_stratum_summary.tsv` — complete six-taxon quartet results.
- `primary_stratum_bootstrap.tsv` — deterministic descriptive gene-level intervals.
- `coordinate_validation_report.md` and `stratum_metadata.tsv` — mapping and biological regime documentation.
- `microbotryum_stage1_topology_by_stratum.*`, `microbotryum_stage1_distance_by_stratum.*`, `microbotryum_stage1_all_quartet_heatmap.*`, and `microbotryum_stage1_representative_trees.*` — Stage 1 figures.

References: [Lucotte et al. 2025](https://doi.org/10.1038/s41467-025-60222-5), [Dryad](https://doi.org/10.5061/dryad.0k6djhb89), [authors' analysis repository](https://github.com/eliselucotte/MicrobotryumOnDianthus).
