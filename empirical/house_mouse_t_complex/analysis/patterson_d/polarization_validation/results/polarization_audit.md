# Correct ancestral-state polarization audit

The primary estimator requires all 8 outgroup genotypes complete and fixed to one nucleotide, inferred from GT calls. Fixed-ALT outgroup sites are polarized by complementing each population ALT frequency (derived REF frequency = 1−ALT frequency). Polymorphic/missing outgroups are excluded. The missing class also records sites with no complete call in at least one ingroup population; partial missing individuals are omitted from that population frequency. Only biallelic one-base A/C/G/T records from the indexed PASS VCF are used. D is computed from pooled site-pattern contributions.

| Design | fixed-ALT additional informative SNPs | old full-chr D | corrected full-chr D | change | old informative | corrected informative |
|---|---:|---:|---:|---:|---:|---:|
| D1_Mmd_t_vs_standard | 141,605 | -0.13038 | +0.19955 | +0.32993 | 139,857 | 281,462 |
| D2_Mmm_t_vs_standard | 239,222 | +0.02569 | +0.13186 | +0.10616 | 124,559 | 363,781 |
| D3_Mmc_t_vs_standard | 305,377 | +0.05725 | +0.08359 | +0.02635 | 360,362 | 665,739 |

## Inversion-union contrasts

| Design | corrected D inside union | full-chromosome complement D | inside−complement | matched-range inside D | matched-range outside D |
|---|---:|---:|---:|---:|---:|
| D1_Mmd_t_vs_standard | +0.39396 | +0.04378 | +0.35019 | +0.35280 | +0.46256 |
| D2_Mmm_t_vs_standard | +0.24055 | +0.02010 | +0.22045 | +0.19674 | +0.34702 |
| D3_Mmc_t_vs_standard | +0.26206 | +0.01446 | +0.24760 | +0.24164 | +0.21648 |

## Polarization class counts and contributions

| Design | class | accepted records | informative sites | ABBA | BABA | D |
|---|---|---:|---:|---:|---:|---:|
| D1_Mmd_t_vs_standard | fixed_ref | 2,784,305 | 139,857 | 10256.253 | 13331.650 | -0.13038 |
| D1_Mmd_t_vs_standard | fixed_alt | 957,554 | 141,605 | 23597.740 | 9258.945 | +0.43640 |
| D1_Mmd_t_vs_standard | polymorphic | 785,558 | 0 | 0.000 | 0.000 | +nan |
| D1_Mmd_t_vs_standard | missing | 558,737 | 0 | 0.000 | 0.000 | +nan |
| D2_Mmm_t_vs_standard | fixed_ref | 2,783,260 | 124,559 | 5309.081 | 5043.089 | +0.02569 |
| D2_Mmm_t_vs_standard | fixed_alt | 956,900 | 239,222 | 33522.124 | 24740.733 | +0.15072 |
| D2_Mmm_t_vs_standard | polymorphic | 785,558 | 0 | 0.000 | 0.000 | +nan |
| D2_Mmm_t_vs_standard | missing | 560,436 | 0 | 0.000 | 0.000 | +nan |
| D3_Mmc_t_vs_standard | fixed_ref | 2,779,117 | 360,362 | 39663.870 | 35368.423 | +0.05725 |
| D3_Mmc_t_vs_standard | fixed_alt | 955,542 | 305,377 | 44760.387 | 36030.090 | +0.10806 |
| D3_Mmc_t_vs_standard | polymorphic | 785,558 | 0 | 0.000 | 0.000 | +nan |
| D3_Mmc_t_vs_standard | missing | 565,937 | 0 | 0.000 | 0.000 | +nan |

Fixed-REF-only and fixed-ALT-only D and contribution counts are in `polarization_class_contributions.tsv`. The corrected full-contig tracks are in `corrected_multiresolution_D.tsv`; the paired corrected/original 5-kb table is `original_vs_corrected_D.tsv`. The full-chromosome complement is descriptive and is not a neutral control. All inversion/tree contrasts in the exact-quartet report use the verified 5–37.055 Mb matched tree region.

Synthetic validation passed for REF/ALT recoding invariance, pure ABBA, pure BABA, balanced site patterns, fixed REF, fixed ALT, polymorphic and missing outgroups, heterozygous and missing ingroup calls, and multiallelic/multibase exclusion. See tests under `../polarization_validation/tests/` and exact quartet scorer tests under `../exact_quartet_comparison/tests/`.
