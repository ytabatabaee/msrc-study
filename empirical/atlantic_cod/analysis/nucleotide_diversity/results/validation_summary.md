# Atlantic cod arrangement-specific diversity: validation summary

## Identifiability conclusion

The VCFs contain diploid per-individual GT, DP, and GQ calls for 29 samples per LG (24 Atlantic cod plus five outgroups). This corrects the prior audit claim that only population summaries or unavailable VCF headers were present. Matschiner et al. (2022) report that 22 additional WGS cod were selected as homokaryotypic at all four supergenes; the gadMor2/gadMor_Stat assembly specimens are the complementary homozygous references. The published cohort mapping supports two classes per LG. X is ancestral; Y is derived. The 24 Atlantic cod samples are mapped by source population/ecotype prefix. PCA/GMM is diagnostic and does not assign classes. The cohort was selected to omit heterokaryotypes, so no XY class or 3-cluster test can be validated here.

Population/ecotype-to-arrangement assignments reuse the source paper convention frozen in population_arrangements.tsv and orientation evidence in arrangement_orientation.tsv. LG01, LG02, and LG07 orientations have direct outgroup/assembly evidence; LG12 derived orientation is demographic inference and therefore moderate-confidence. Individual calls do not establish a third arrangement class. See individual_genotype_classification.tsv and individual_arrangement_assignments.tsv.

## VCF QC and PCA

All 24 Atlantic cod individuals were used for PCA after PASS/QUAL≥30, GQ≥20, DP≥8, site call-rate ≥90%, and MAF ≥0.02; missing genotypes at retained SNPs were mean-imputed for PCA. PCA uses no LD pruning and is descriptive. For all LGs, GMM BIC over PC1–PC5 was lowest at k=4; this is not interpreted as four arrangement states and likely reflects within-background population structure.

| LG | VCF samples | cod samples | SNP records | mean GQ | mean DP | genotype missing after QC | PCA SNPs | PC1 variance | PC2 variance | best exploratory k |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LG01 | 29 | 24 | 11474 | 48.16 | 16.85 | 0.0252 | 3630 | 0.280 | 0.067 | 4 |
| LG02 | 29 | 24 | 3123 | 47.64 | 16.55 | 0.0271 | 1060 | 0.354 | 0.053 | 4 |
| LG07 | 29 | 24 | 10411 | 48.68 | 17.06 | 0.0253 | 3930 | 0.421 | 0.053 | 4 |
| LG12 | 29 | 24 | 10339 | 48.24 | 16.83 | 0.0251 | 2731 | 0.200 | 0.068 | 4 |

The GMM minimum BIC did not identify three clusters. The expected arrangement classes in these selected samples are XX and YY only; ordinary SNP heterozygosity is not arrangement heterozygosity. PCA plots are colored by population and annotate individual IDs.

## Window estimates and missingness

Filters: biallelic SNPs, FILTER PASS or ., QUAL≥30; called diploid genotypes require GQ≥20 and DP≥8. The common eligible panel requires at least four observed haploid copies in both backgrounds (at least two diploid individuals per class). The identical sites enter πX, πY, and dXY. Within-background π uses the unbiased finite-sample correction; dXY uses cross-class pairwise differences. FST is 1−(πX+πY)/(2dXY), not truncated, and undefined at dXY=0. Values are SNP-conditional, not per-base; no callable invariant mask exists.

| LG | XX (ancestral) | YY (derived) | fully inside windows with estimates | boundary-overlap windows (NA) | fully outside windows (NA) | mean πX SNP | mean πY SNP | mean dXY SNP | mean Hudson FST | minimum eligible SNPs/window |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LG01 | 17 | 7 | 66/66 | 1 | 43 | 0.07114 | 0.03963 | 0.07644 | 0.27111 | 38 |
| LG02 | 17 | 7 | 21/21 | 1 | 70 | 0.06013 | 0.03833 | 0.12996 | 0.57843 | 33 |
| LG07 | 17 | 7 | 36/36 | 2 | 81 | 0.08629 | 0.03167 | 0.14575 | 0.59157 | 118 |
| LG12 | 17 | 7 | 53/53 | 1 | 51 | 0.04323 | 0.04315 | 0.06233 | 0.28774 | 15 |

In total, 176/176 fully contained suppressed-recombination windows have estimates: LG01 66, LG02 21, LG07 36, LG12 53. There are 245 fully outside and five boundary-overlap windows; all remain NA because ancestral/derived background labels are not extended through recombining sequence or across mixed windows. No inside window was zero-SNP or fell below the two-individual-per-class threshold. Missingness is observed-call missingness after GQ/DP filtering; it is not a callable-base denominator.

## Comparison with existing tracks

Correlations use only fully inside windows with estimates, against the existing stage4a arrangement-quartet fraction and mean opposite-arrangement stage7 pairwise MRCA-time shift. These are exploratory spatial associations, not independent tests; strong autocorrelation makes nominal p-values descriptive. The LG01 quartet-fraction track is constant over eligible windows, so its correlations are undefined. Most coefficients are modest; the clearest positive associations are on LG12, where πX, πY, dXY, and FST correlate positively with both existing tracks (Spearman rho about 0.32–0.54). On LG02, dXY/FST track branch-time shifts positively (rho 0.66/0.37); LG07 dXY has a modest positive association with branch-time shift (rho 0.39). Full coefficients and window counts are in comparison_with_quartet_branch_tracks.tsv.

## Limitations

1. VCF records are variant-only; no callable invariant-site mask or per-base π/dXY is available. 2. The 24-individual cohort was selected for homozygous arrangement backgrounds and cannot reveal arrangement heterozygotes or estimate XY diversity. 3. Arrangement states are mapped from documented sampling ecotype and published cohort selection; exact external specimen accession-to-VCF-alias mapping is not separately deposited locally. 4. LG12 orientation is demographic inference, and published Bornholm samples carry a roughly 275-kb tract of mixed local ancestry. 5. SNP-conditional tracks describe the selected sample backgrounds; population structure is not separated from arrangement effects. 6. Estimates outside suppressed-recombination regions and across boundaries are intentionally not identifiable.

## Inside versus outside cohort comparison

The original arrangement-specific tracks remain limited to fully inside windows. For a direct regional comparison, the same 17 ancestral-associated and seven derived-associated samples per LG are also grouped across the flanking/outside SNP VCF. Those outside statistics are population-group summaries, not arrangement-linked backgrounds because recombination decouples them from the inversion. This comparison is descriptive and confounded with population structure/ecotype.

| LG | scope | valid windows | mean pi_X group | mean pi_Y group | mean dXY group | mean Hudson FST |
|---|---|---:|---:|---:|---:|---:|
| LG01 | inside | 66 | 0.07114 | 0.03963 | 0.07644 | 0.27111 |
| LG01 | outside | 43 | 0.05490 | 0.05638 | 0.05608 | 0.00825 |
| LG02 | inside | 21 | 0.06013 | 0.03833 | 0.12996 | 0.57843 |
| LG02 | outside | 70 | 0.06110 | 0.05957 | 0.06121 | 0.01362 |
| LG07 | inside | 36 | 0.08629 | 0.03167 | 0.14575 | 0.59157 |
| LG07 | outside | 81 | 0.05798 | 0.05643 | 0.05785 | 0.00819 |
| LG12 | inside | 53 | 0.04323 | 0.04315 | 0.06233 | 0.28774 |
| LG12 | outside | 51 | 0.06112 | 0.05927 | 0.06210 | 0.02179 |

Full window values and interpretation labels are in inside_outside_cohort_group_windows.tsv; the separate inside arrangement-specific tracks retain NA outside.
