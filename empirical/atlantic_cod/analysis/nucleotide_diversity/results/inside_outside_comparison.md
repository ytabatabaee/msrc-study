# Inside versus outside comparison

The primary arrangement-specific table continues to report pi_X, pi_Y, dXY, and FST only within fully contained inversion windows. This additional comparison uses the same per-LG ancestral-associated (17 individuals) and derived-associated (7 individuals) cohorts throughout the chromosome. Inside suppressed regions, those cohorts correspond to homozygous arrangement backgrounds based on the published study selection. Outside, values are population-grouped summaries only; recombination means their alleles cannot be assigned to the inversion arrangement. The outside results therefore describe group differentiation, not continued X/Y chromosomal backgrounds.

| LG | inside windows | outside windows | inside πX | outside πX group | inside πY | outside πY group | inside dXY | outside dXY group | inside FST | outside FST group |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LG01 | 66 | 43 | 0.07114 | 0.05490 | 0.03963 | 0.05638 | 0.07644 | 0.05608 | 0.27111 | 0.00825 |
| LG02 | 21 | 70 | 0.06013 | 0.06110 | 0.03833 | 0.05957 | 0.12996 | 0.06121 | 0.57843 | 0.01362 |
| LG07 | 36 | 81 | 0.08629 | 0.05798 | 0.03167 | 0.05643 | 0.14575 | 0.05785 | 0.59157 | 0.00819 |
| LG12 | 53 | 51 | 0.04323 | 0.06112 | 0.04315 | 0.05927 | 0.06233 | 0.06210 | 0.28774 | 0.02179 |

All estimates use the same biallelic PASS/dot, QUAL≥30, GQ≥20, DP≥8 filters and at least four observed haploid copies in each group at a site. The same eligible SNP set defines within-group pi and between-group dXY at each window. Values are SNP-conditional because the source VCFs are variant-only and no callable invariant mask is available. Boundary-crossing windows remain NA. Summary values are equal-window averages; do not interpret the contrast causally because population structure and ecotype are confounded with the background cohorts outside inversions. Numeric outside-minus-inside differences are in inside_outside_cohort_contrast.tsv. Full window values are in inside_outside_cohort_group_summary.tsv and inside_outside_cohort_group_windows.tsv.
