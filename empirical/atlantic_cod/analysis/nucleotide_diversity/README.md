# Atlantic cod arrangement-specific nucleotide diversity

Run from repository root:

    python3 empirical/atlantic_cod/analysis/nucleotide_diversity/scripts/audit_and_report.py
    python3 -m unittest discover -s empirical/atlantic_cod/analysis/nucleotide_diversity/tests -v

## Identifiability and classification

The four inversion-region SNAPP VCFs contain individual diploid GT, DP, GQ, and site QUAL fields: 29 sample columns per file, including 24 Atlantic cod and five outgroups. These are not population summaries. The earlier audit incorrectly said local VCF headers were unavailable; direct inspection corrected that.

Matschiner et al. (2022) report that their 22 additional WGS cod were selected after preliminary analysis showed each was homokaryotypic at all four supergene regions. The gadMor2 and gadMor_Stat assembly specimens are the complementary homozygous references. The local VCFs contain the 24 Atlantic cod individuals from this cohort; their VCF names map to source population/ecotype prefixes. We map each sample to ancestral/derived state using data/atlantic_cod/processed/population_arrangements.tsv and orientation evidence in data/atlantic_cod/metadata/arrangement_orientation.tsv. This is literature-supported cohort assignment, not an assumption that arbitrary population labels are karyotypes. X denotes ancestral and Y derived. individual_genotype_classification.tsv records per-sample call rate, heterozygosity, and putative XX/YY. individual_arrangement_assignments.tsv records evidence and confidence per LG.

The study selected homozygotes, so this dataset cannot estimate an XY class or validate three clusters. PCA plots show individual genotype structure colored by population; GMM BIC over PC1–PC5 is descriptive, and k=4 has the lowest BIC for every LG. That does not indicate four arrangements: PCA is not used to assign genotypes. LG12 orientation confidence is moderate because ancestral/derived orientation was inferred demographically. The study also reports two Bornholm individuals with a roughly 275-kb tract of heterozygous haplotype-informative sites in LG12; local ancestry is mixed there despite their overall ancestral chromosome-arrangement assignment.

## Estimator and scope

The analysis filters to biallelic SNPs with FILTER PASS or ., QUAL ≥30, and individual calls with GQ ≥20 and DP ≥8. A site enters the common eligible panel when each background contributes at least four called haploid copies (at least two diploid individuals). This same site set is used for πX, πY, and dXY within each window. Within-background π uses 2k(n-k)/(n(n-1)), equivalent to the finite-sample correction. dXY is the cross-background pairwise allele-difference estimate. Hudson FST is 1 - (πX + πY)/(2 dXY); negative values are retained, and FST is NA when dXY=0.

The outputs are SNP-conditional because the VCFs contain variant records only and no callable invariant-site mask was found. Estimates are never labelled per-base. Only 250-kb windows fully inside frozen suppressed-recombination intervals receive estimates. Windows outside or crossing boundaries remain NA; arrangement identity is not extrapolated through recombining sequence. Within the LG12 Bornholm ancestry tract, arrangement assignment does not imply every nucleotide has that arrangement's local ancestry.

## Outputs

- results/atlantic_cod_pi_windows.tsv and results/atlantic_cod_pi_dxy_fst_windows.tsv: every existing window, sample sizes, SNP denominators, missingness, SNP-conditional πX, πY, dXY, and Hudson FST.
- results/genotype_qc.tsv, results/individual_genotype_classification.tsv, results/individual_arrangement_assignments.tsv: genotype QC and per-sample classes/evidence.
- results/genotype_pca_lg_summary.tsv, figures/LG01_pca.png through LG12_pca.png: PCA and descriptive cluster diagnostics.
- figures/LG01_pi_dxy_fst_tracks.png/.pdf through LG12_pi_dxy_fst_tracks.png/.pdf: chromosome tracks; outside estimates remain NA.
- results/comparison_with_quartet_branch_tracks.tsv: exploratory within-inversion correlations with quartet arrangement fraction and opposite-arrangement pairwise MRCA-time shifts.
- results/inside_outside_cohort_group_windows.tsv and inside_outside_cohort_group_summary.tsv: secondary comparison using the same 17 ancestral-associated and seven derived-associated samples. Inside windows are homozygous arrangement backgrounds; outside windows are population-group summaries only because arrangement linkage is not assumed. The contrast is descriptive and potentially confounded by population/ecotype.
- results/inside_outside_comparison.md and figures/atlantic_cod_inside_outside_cohort_comparison.png/.pdf: inside/outside summary table and figure.
- results/track_validation_report.md: assignment evidence and estimability summary.

Correlations are descriptive because the spatial windows are not independent. The quartet/branch correlations use only inside windows with valid diversity. The separate inside/outside cohort comparison estimates outside population-group summaries, not arrangement-specific values. VCF genotype filters approximate callability at recorded SNPs; they do not provide callable-base denominators.

Reference: Matschiner et al. (2022), “Supergene origin and maintenance in Atlantic cod,” Nature Ecology & Evolution 6:469–481, https://doi.org/10.1038/s41559-022-01661-x. The paper describes the homokaryotype selection (Results/Methods), the gadMor2/gadMor_Stat arrangement orientations, and the LG12 Bornholm ancestry tract (Fig. 5).
