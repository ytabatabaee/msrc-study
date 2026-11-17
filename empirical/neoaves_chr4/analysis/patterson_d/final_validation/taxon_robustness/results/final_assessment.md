# Final assessment: N61 alternative-taxon validation

## Recommendation: CONDITIONAL GO

The N61 inside/outside nucleotide D contrast is manuscript-ready as a carefully delimited descriptive result for the tested taxa. It is robust in direction to the ten prespecified ingroup representative alternatives: 10/10 retain negative ΔD and 10/10 have 1-Mb block-jackknife 95% intervals below zero. Transversion-only ΔD is negative with its 95% interval below zero in 10/10 alternatives. This supports a regional site-pattern contrast in this chr4 dataset. It does not distinguish introgression, ancestral population structure, ILS, deep polymorphism, branch-length/substitution effects, homoplasy, or taxon-specific sampling, and it does not establish the MSRC mechanism.

The whole-chromosome nucleotide–genealogy correlation is positive for the original and all tested alternatives (1-Mb r range 0.837–0.929); the original r=0.914 falls within that range. After inside/outside mean differences are removed, correlations remain positive (range 0.540–0.761), but they are descriptive spatial associations, not independent tests of mechanism. Bootstrap intervals, region-specific correlations, and all resolutions are in `nucleotide_topology_robustness.tsv`.

## Alternative taxon panel and null topology

The original quartet remains primary: Tauraco_erythrolophus / Columba_livia / Ciconia_maguari / Gallus_gallus. The frozen panel includes two P1 replacements (Cuculus, Ardeotis), two P2 replacements (Pterocles, Mesitornis), three P3 alternatives/controls (Balaeniceps, Nipponia, and a Telluraves Aquila control), three combined ingroup substitutions, and one outgroup-only Meleagris sensitivity. All ten ingroup alternatives were selected from taxonomic/literature groupings and archive availability before D calculation; no quartet was selected by its D value. The design table gives each rationale.

The exact induced quartet in the Stiller 363-taxon published reference is q1 for all twelve designs, including the Gallus/Meleagris pair. All resolved local gene trees use exact taxon names and locus IDs; unresolved quartet trees are counted separately (0 across the exact quartet-bearing chr4 subset here). The reference tree comes from the same broad Stiller study as the alignments and is not an independent species-history replicate. Jarvis et al.'s 48-taxon TENT supports broad historical context but lacks *Ciconia maguari*, so an exact independent four-taxon species-tree check is unavailable.

Taxonomic clarification: the frozen P3 taxon *Ciconia maguari* is classified in Aequornithes/Ciconiiformes, not Telluraves. Therefore, P3 is described by its position on the tested quartet's non-P1/P2 side. *Balaeniceps* and *Nipponia* are waterbird-side alternatives; *Aquila* is a labeled broad Telluraves control. Their q1 placement is verified in the Stiller reference, but they are not assumed to be interchangeable members of one narrow monophyletic order.

## Complete regional D results

**Prespecified criterion:** an alternative retains the primary regional direction if all-sites ΔD<0; interval support is summarized by the 1-Mb paired physical-block delete-one jackknife 95% interval, with ‘excludes zero’ meaning its upper endpoint is below 0. This rule is stated before interpreting the results and does not count loci or quartets as independent replicates.

All estimates pool integer ABBA/BABA counts before calculating D. Block intervals use genomic blocks; informative sites are not independent replicates. The table includes the original, every ingroup alternative, and the alternate-outgroup sensitivity. Locus/site counts are inside/outside the union of the three source-coordinate PNAS intervals. Detailed estimates for each PNAS interval, outside-union-only estimates, and 500-kb/2-Mb block results are in `taxon_robustness_block_summary.tsv`.

| Frozen design | D inside | D outside | ΔD (95% 1-Mb block interval) | Loci in/out | Informative sites in/out | Transversion-only ΔD (95% interval) |
|---|---:|---:|---:|---:|---:|---:|
| N61_primary | -0.258 | -0.081 | -0.178 [-0.198, -0.157] | 1221/4427 | 24683/73729 | -0.318 [-0.365, -0.271] |
| N61_P1_Cuculus | -0.123 | 0.056 | -0.179 [-0.205, -0.154] | 1206/4290 | 25786/78254 | -0.338 [-0.387, -0.289] |
| N61_P1_Ardeotis | -0.262 | -0.070 | -0.191 [-0.217, -0.166] | 1222/4371 | 24234/72350 | -0.345 [-0.383, -0.306] |
| N61_P2_Pterocles | -0.234 | -0.015 | -0.218 [-0.240, -0.197] | 1186/4446 | 23066/70522 | -0.371 [-0.420, -0.323] |
| N61_P2_Mesitornis | -0.263 | -0.062 | -0.201 [-0.220, -0.181] | 1229/4417 | 24456/72610 | -0.340 [-0.385, -0.295] |
| N61_P3_Balaeniceps | -0.258 | -0.084 | -0.174 [-0.197, -0.152] | 1049/3864 | 21082/63592 | -0.311 [-0.355, -0.266] |
| N61_P3_Nipponia | -0.269 | -0.073 | -0.196 [-0.214, -0.178] | 1216/4418 | 24927/74527 | -0.330 [-0.370, -0.290] |
| N61_P3_Aquila | -0.255 | -0.071 | -0.184 [-0.207, -0.161] | 1219/4427 | 24850/74134 | -0.329 [-0.371, -0.288] |
| N61_combo_Cuculus_Pterocles | -0.086 | 0.117 | -0.204 [-0.226, -0.181] | 1171/4306 | 24457/75819 | -0.377 [-0.428, -0.327] |
| N61_combo_Ardeotis_Mesitornis | -0.262 | -0.052 | -0.210 [-0.231, -0.188] | 1233/4360 | 24097/71285 | -0.366 [-0.407, -0.324] |
| N61_combo_Cuculus_Pterocles_Balaeniceps | -0.103 | 0.111 | -0.214 [-0.238, -0.189] | 1011/3751 | 20922/65447 | -0.387 [-0.435, -0.339] |
| N61_primary_Meleagris | -0.253 | -0.086 | -0.167 [-0.186, -0.147] | 1206/4376 | 24243/71749 | -0.307 [-0.354, -0.260] |

## Exact genealogy comparison

`q1/q2/q3` are observed proportions across the exact resolved chr4 quartet trees; D_topology=q3−q2. Spatial correlations pair each nucleotide locus to its exact tree by identical archive locus ID, then pool into 250-kb, 500-kb, and 1-Mb bins. Block-bootstrap intervals resample 2-Mb physical blocks (1,000 replicates, seed 20261118). The result table includes whole chromosome, inside, outside, and within-region-centered correlations.

| Design | Exact trees | Unresolved | Mean q1/q2/q3 | Mean D_topology | 1-Mb r whole (95% block CI) | r after region-centering (95% block CI) |
|---|---:|---:|---:|---:|---:|---:|
| N61_primary | 5648 | 0 | 0.284/0.473/0.243 | -0.230 | 0.914 [0.864, 0.945] | 0.732 [0.586, 0.834] |
| N61_P1_Cuculus | 5496 | 0 | 0.298/0.450/0.252 | -0.198 | 0.837 [0.707, 0.909] | 0.540 [0.216, 0.737] |
| N61_P1_Ardeotis | 5593 | 0 | 0.293/0.460/0.247 | -0.212 | 0.876 [0.804, 0.919] | 0.608 [0.382, 0.748] |
| N61_P2_Pterocles | 5632 | 0 | 0.269/0.473/0.258 | -0.215 | 0.892 [0.832, 0.923] | 0.612 [0.244, 0.740] |
| N61_P2_Mesitornis | 5646 | 0 | 0.274/0.472/0.254 | -0.217 | 0.871 [0.787, 0.924] | 0.572 [0.224, 0.729] |
| N61_P3_Balaeniceps | 4913 | 0 | 0.279/0.472/0.249 | -0.223 | 0.881 [0.797, 0.928] | 0.671 [0.440, 0.800] |
| N61_P3_Nipponia | 5634 | 0 | 0.285/0.471/0.244 | -0.226 | 0.929 [0.885, 0.954] | 0.761 [0.499, 0.866] |
| N61_P3_Aquila | 5646 | 0 | 0.285/0.456/0.260 | -0.196 | 0.886 [0.825, 0.922] | 0.604 [0.407, 0.729] |
| N61_combo_Cuculus_Pterocles | 5477 | 0 | 0.286/0.450/0.264 | -0.187 | 0.901 [0.829, 0.940] | 0.646 [0.424, 0.780] |
| N61_combo_Ardeotis_Mesitornis | 5593 | 0 | 0.281/0.463/0.256 | -0.207 | 0.874 [0.786, 0.926] | 0.592 [0.291, 0.746] |
| N61_combo_Cuculus_Pterocles_Balaeniceps | 4762 | 0 | 0.285/0.454/0.261 | -0.192 | 0.887 [0.801, 0.927] | 0.624 [0.395, 0.759] |
| N61_primary_Meleagris | 5582 | 0 | 0.284/0.473/0.243 | -0.230 | 0.900 [0.845, 0.938] | 0.686 [0.519, 0.808] |

Among ingroup alternatives, the full-chromosome r range is 0.837–0.929; the original is not unusually high. The centered correlation is also positive throughout this panel, but it is lower than whole-chromosome r because regional mean differences explain part of the association. Shared loci and overlapping taxa make these sensitivity designs correlated; do not count them as independent biological replicates.

## Independent raw-alignment validation

A separate parser and direct pattern counter independently decoded 46 deterministic, stratified raw alignments across PNAS_1, PNAS_2, PNAS_3, outside high/low-information loci, and loci with ambiguous/gap data. It compared the three original primary quartets and their predeclared Meleagris outgroup designs against saved all-site and transversion-only counts: **544 exact ABBA/BABA matches, 0 discrepancies, 4 unavailable rows**. The unavailable N62 rows are the two sampled loci missing *Pterocles gutturalis* under each of the two N62 outgroups. There were 0 sampled loci with duplicate sequence IDs and 0 with unequal full-alignment sequence lengths. No production filters were changed. See `raw_alignment_validation_report.md` and row-level `raw_alignment_validation.tsv`.

No estimator error was found, so the frozen primary results did not require reanalysis. Regression tests cover parsing, duplicate IDs, case normalization, gaps/ambiguity, multistate filtering, outgroup polarization, ABBA/BABA orientation, and transition/transversion classification.

## Answers to the requested scientific questions

1. **Alternative ingroup quartets:** 10 biologically motivated sensitivity designs were tested, in addition to the original primary and an outgroup-only sensitivity.
2. **Original regional direction:** 10/10 alternatives retain negative all-sites ΔD.
3. **Block uncertainty excludes zero:** 10/10 alternatives have a 1-Mb 95% paired block-jackknife interval wholly below zero. Transversion-only direction and interval are also negative in 10/10.
4. **Dependence on P1/P2/P3:** the contrast direction is not dependent on a single representative. The absolute inside/outside D background is taxon-sensitive: e.g., replacing P1 with *Cuculus* and replacing P1+P2 together changes the outside estimate's sign while leaving ΔD negative.
5. **D_site–D_topology correspondence:** whole-chromosome correspondence remains positive and of similar magnitude across alternatives; the original r is typical. Within-region-centered values are smaller but positive across the panel. These are block-aware descriptive correlations.
6. **Independent raw alignment decode:** exact integer ABBA/BABA agreement for every available sampled design/filter row; 0 discrepancies.
7. **Columbea/N62:** no new full ingroup sensitivity was run for those hypotheses. Their previous near-zero regional contrasts and interpretation remain unchanged; the raw validation checks the original saved counts but does not re-estimate their full archives.
8. **Remaining limitations:** exact species-history quartets are not independently verified for every alternative; the tree reference and alignments share the Stiller source; Gallus and Meleagris are related Galliformes; raw decoding is stratified sampling rather than a full independent reanalysis of all 6,247 alignments; block bootstrap/jackknife is descriptive; alternative designs share linked loci; no demographic/finite-sites model, branch-length-time conversion, or spatial-null test; D/topology cannot resolve introgression versus ILS, ancestral structure, deep polymorphism, homoplasy, or mechanism.
9. **Manuscript readiness:** **CONDITIONAL GO** for the narrow N61 regional site-pattern result and its descriptive genealogy concordance, with the nonzero outside baseline and taxon sensitivity stated. Do not present these results as unique evidence for introgression, an inversion, or MSRC.

## Figures

- Figure E: `figures/figure_E_alternative_quartet_delta_D` (PNG, PDF, SVG)
- Figure F: `figures/figure_F_site_vs_topology` (PNG, PDF, SVG)
- Figure G: `figures/figure_G_sensitivity_matrix` (PNG, PDF, SVG)
