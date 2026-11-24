# Individual genotype and arrangement identifiability audit

The inversion VCFs contain diploid individual genotype calls (GT), genotype quality (GQ), read depth (DP), and site-level QUAL; they are not population summaries. This corrects the earlier local audit, which had not inspected the available headers. Arrangement classes are assigned from the published homokaryotype-selected cohort and documented population/ecotype mapping; PCA/GMM do not assign classes.

## LG01

- Samples: 29 total VCF columns; 24 Gadmor cod individuals used in PCA.
- Biallelic SNP records in inversion VCF: 11474; PASS/QUAL≥30: 11474; PCA SNPs after ≥90% call rate and MAF≥0.02 with GQ≥20, DP≥8: 3630.
- Raw missing genotype fraction: 0.0; after genotype filters: 0.02522644900314354; mean QUAL: 2609.629097960607; mean GQ: 48.16006202929562; mean DP: 16.84881260781497.
- PC variance: 0.2801506073646739, 0.06697764414704765; exploratory GMM BIC (PC1–PC5): {1: 1063.02117341888, 2: 782.5431639013267, 3: 677.489848339767, 4: 608.9255815304068}; lowest BIC k=4.
- Interpretation: PCA/clustering describes SNP structure. All 24 cod individuals are part of the published cohort: 22 additional individuals were selected as homokaryotypic for all four supergenes, and the two assembly specimens were selected/identified as complementary homozygotes. We therefore map their VCF population/ecotype IDs to published arrangement classes; PCA is not used to define those classes. No XY class is expected because the cohort was ascertained to exclude arrangement heterokaryotypes.

## LG02

- Samples: 29 total VCF columns; 24 Gadmor cod individuals used in PCA.
- Biallelic SNP records in inversion VCF: 3123; PASS/QUAL≥30: 3123; PCA SNPs after ≥90% call rate and MAF≥0.02 with GQ≥20, DP≥8: 1060.
- Raw missing genotype fraction: 0.0; after genotype filters: 0.02714012830280345; mean QUAL: 2679.5451809157867; mean GQ: 47.63665573553281; mean DP: 16.55365640906732.
- PC variance: 0.3537052724219268, 0.05319249953642574; exploratory GMM BIC (PC1–PC5): {1: 913.1039949181885, 2: 656.4955791581577, 3: 629.5507904489191, 4: 595.521238837666}; lowest BIC k=4.
- Interpretation: PCA/clustering describes SNP structure. All 24 cod individuals are part of the published cohort: 22 additional individuals were selected as homokaryotypic for all four supergenes, and the two assembly specimens were selected/identified as complementary homozygotes. We therefore map their VCF population/ecotype IDs to published arrangement classes; PCA is not used to define those classes. No XY class is expected because the cohort was ascertained to exclude arrangement heterokaryotypes.

## LG07

- Samples: 29 total VCF columns; 24 Gadmor cod individuals used in PCA.
- Biallelic SNP records in inversion VCF: 10411; PASS/QUAL≥30: 10411; PCA SNPs after ≥90% call rate and MAF≥0.02 with GQ≥20, DP≥8: 3930.
- Raw missing genotype fraction: 0.0; after genotype filters: 0.025251805947952927; mean QUAL: 3154.6562943041013; mean GQ: 48.68177557556828; mean DP: 17.063305058641557.
- PC variance: 0.42080773216134626, 0.05349677417804503; exploratory GMM BIC (PC1–PC5): {1: 1057.7513317480439, 2: 716.2212949127974, 3: 646.7193627066481, 4: 613.32466855671}; lowest BIC k=4.
- Interpretation: PCA/clustering describes SNP structure. All 24 cod individuals are part of the published cohort: 22 additional individuals were selected as homokaryotypic for all four supergenes, and the two assembly specimens were selected/identified as complementary homozygotes. We therefore map their VCF population/ecotype IDs to published arrangement classes; PCA is not used to define those classes. No XY class is expected because the cohort was ascertained to exclude arrangement heterokaryotypes.

## LG12

- Samples: 29 total VCF columns; 24 Gadmor cod individuals used in PCA.
- Biallelic SNP records in inversion VCF: 10339; PASS/QUAL≥30: 10339; PCA SNPs after ≥90% call rate and MAF≥0.02 with GQ≥20, DP≥8: 2731.
- Raw missing genotype fraction: 0.0; after genotype filters: 0.02508746593914572; mean QUAL: 2487.51689428378; mean GQ: 48.241502713195096; mean DP: 16.827489485743637.
- PC variance: 0.2003779949131612, 0.06798375860151987; exploratory GMM BIC (PC1–PC5): {1: 1028.7860646544725, 2: 836.1989038333833, 3: 740.324885195052, 4: 672.6237380100362}; lowest BIC k=4.
- Interpretation: PCA/clustering describes SNP structure. All 24 cod individuals are part of the published cohort: 22 additional individuals were selected as homokaryotypic for all four supergenes, and the two assembly specimens were selected/identified as complementary homozygotes. We therefore map their VCF population/ecotype IDs to published arrangement classes; PCA is not used to define those classes. No XY class is expected because the cohort was ascertained to exclude arrangement heterokaryotypes.

All 24 cod samples map to the published cohort: 22 selected homokaryotypes plus the two complementary assembly specimens. The cohort has no expected XY class because heterokaryotypes were excluded by study design. The PCA does not recover three clusters; its lowest GMM BIC among k=1–4 is k=4 for each LG, likely reflecting population structure, so it is not interpreted as arrangement count. SNP-conditional tracks are estimated separately in estimate_tracks.py; no callable invariant mask exists, and windows outside or crossing suppressed-region boundaries remain unassigned.
