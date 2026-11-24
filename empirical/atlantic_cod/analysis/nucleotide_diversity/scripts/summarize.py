#!/usr/bin/env python3
"""Build final validation reports from generated TSVs."""
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[5]
BASE=ROOT/'empirical/atlantic_cod/analysis/nucleotide_diversity'
R=BASE/'results'
def read(name):
 with (R/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def mean(rows,col):
 a=[float(x[col]) for x in rows if x[col] not in ('NA','nan','')]
 return sum(a)/len(a) if a else float('nan')
rows=read('atlantic_cod_pi_dxy_fst_windows.tsv')
qc=read('genotype_qc.tsv')
pca=read('genotype_pca_lg_summary.tsv')
lines=['# Atlantic cod arrangement-specific diversity: validation summary','',
'## Identifiability conclusion','',
'The VCFs contain diploid per-individual GT, DP, and GQ calls for 29 samples per LG (24 Atlantic cod plus five outgroups). This corrects the prior audit claim that only population summaries or unavailable VCF headers were present. Matschiner et al. (2022) report that 22 additional WGS cod were selected as homokaryotypic at all four supergenes; the gadMor2/gadMor_Stat assembly specimens are the complementary homozygous references. The published cohort mapping supports two classes per LG. X is ancestral; Y is derived. The 24 Atlantic cod samples are mapped by source population/ecotype prefix. PCA/GMM is diagnostic and does not assign classes. The cohort was selected to omit heterokaryotypes, so no XY class or 3-cluster test can be validated here.','',
'Population/ecotype-to-arrangement assignments reuse the source paper convention frozen in population_arrangements.tsv and orientation evidence in arrangement_orientation.tsv. LG01, LG02, and LG07 orientations have direct outgroup/assembly evidence; LG12 derived orientation is demographic inference and therefore moderate-confidence. Individual calls do not establish a third arrangement class. See individual_genotype_classification.tsv and individual_arrangement_assignments.tsv.','',
'## VCF QC and PCA','',
'All 24 Atlantic cod individuals were used for PCA after PASS/QUAL≥30, GQ≥20, DP≥8, site call-rate ≥90%, and MAF ≥0.02; missing genotypes at retained SNPs were mean-imputed for PCA. PCA uses no LD pruning and is descriptive. For all LGs, GMM BIC over PC1–PC5 was lowest at k=4; this is not interpreted as four arrangement states and likely reflects within-background population structure.','',
'| LG | VCF samples | cod samples | SNP records | mean GQ | mean DP | genotype missing after QC | PCA SNPs | PC1 variance | PC2 variance | best exploratory k |','|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
for q,p in zip(qc,pca):
 lines.append(f"| {q['lg']} | {q['n_vcf_samples']} | {q['n_codon_samples_assigned']} | {q['PASS_biallelic_SNPs_QUAL30']} | {float(q['mean_GQ']):.2f} | {float(q['mean_DP']):.2f} | {float(q['GQ_DP_filtered_missing_fraction']):.4f} | {p['pca_snps']} | {float(p['pc1_variance']):.3f} | {float(p['pc2_variance']):.3f} | {p['best_exploratory_k']} |")
lines+=['','The GMM minimum BIC did not identify three clusters. The expected arrangement classes in these selected samples are XX and YY only; ordinary SNP heterozygosity is not arrangement heterozygosity. PCA plots are colored by population and annotate individual IDs.','',
'## Window estimates and missingness','',
'Filters: biallelic SNPs, FILTER PASS or ., QUAL≥30; called diploid genotypes require GQ≥20 and DP≥8. The common eligible panel requires at least four observed haploid copies in both backgrounds (at least two diploid individuals per class). The identical sites enter πX, πY, and dXY. Within-background π uses the unbiased finite-sample correction; dXY uses cross-class pairwise differences. FST is 1−(πX+πY)/(2dXY), not truncated, and undefined at dXY=0. Values are SNP-conditional, not per-base; no callable invariant mask exists.','',
'| LG | XX (ancestral) | YY (derived) | fully inside windows with estimates | boundary-overlap windows (NA) | fully outside windows (NA) | mean πX SNP | mean πY SNP | mean dXY SNP | mean Hudson FST | minimum eligible SNPs/window |','|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
for q in qc:
 lg=q['lg'];g=[x for x in rows if x['chromosome']==lg]
 inside=[x for x in g if x['inside_suppressed_recombination']=='true']
 valid=[x for x in inside if x['pi_X_SNP']!='NA']
 boundary=[x for x in g if x['overlaps_suppressed_recombination']=='true' and x['inside_suppressed_recombination']=='false']
 outside=[x for x in g if x['overlaps_suppressed_recombination']=='false']
 lines.append(f"| {lg} | {q['n_X_individuals']} | {q['n_Y_individuals']} | {len(valid)}/{len(inside)} | {len(boundary)} | {len(outside)} | {mean(valid,'pi_X_SNP'):.5f} | {mean(valid,'pi_Y_SNP'):.5f} | {mean(valid,'dXY_SNP'):.5f} | {mean(valid,'FST_Hudson'):.5f} | {min(int(x['joint_eligible_SNPs']) for x in valid)} |")
lines+=['','In total, 176/176 fully contained suppressed-recombination windows have estimates: LG01 66, LG02 21, LG07 36, LG12 53. There are 245 fully outside and five boundary-overlap windows; all remain NA because ancestral/derived background labels are not extended through recombining sequence or across mixed windows. No inside window was zero-SNP or fell below the two-individual-per-class threshold. Missingness is observed-call missingness after GQ/DP filtering; it is not a callable-base denominator.','',
'## Comparison with existing tracks','',
'Correlations use only fully inside windows with estimates, against the existing stage4a arrangement-quartet fraction and mean opposite-arrangement stage7 pairwise MRCA-time shift. These are exploratory spatial associations, not independent tests; strong autocorrelation makes nominal p-values descriptive. The LG01 quartet-fraction track is constant over eligible windows, so its correlations are undefined. Most coefficients are modest; the clearest positive associations are on LG12, where πX, πY, dXY, and FST correlate positively with both existing tracks (Spearman rho about 0.32–0.54). On LG02, dXY/FST track branch-time shifts positively (rho 0.66/0.37); LG07 dXY has a modest positive association with branch-time shift (rho 0.39). Full coefficients and window counts are in comparison_with_quartet_branch_tracks.tsv.','',
'## Limitations','',
'1. VCF records are variant-only; no callable invariant-site mask or per-base π/dXY is available. 2. The 24-individual cohort was selected for homozygous arrangement backgrounds and cannot reveal arrangement heterozygotes or estimate XY diversity. 3. Arrangement states are mapped from documented sampling ecotype and published cohort selection; exact external specimen accession-to-VCF-alias mapping is not separately deposited locally. 4. LG12 orientation is demographic inference, and published Bornholm samples carry a roughly 275-kb tract of mixed local ancestry. 5. SNP-conditional tracks describe the selected sample backgrounds; population structure is not separated from arrangement effects. 6. Estimates outside suppressed-recombination regions and across boundaries are intentionally not identifiable.','']
for name in ['track_validation_report.md','validation_summary.md']:
 (R/name).write_text('\n'.join(lines))
comparison=R/'inside_outside_cohort_group_summary.tsv'
if comparison.exists():
 group_rows=[]
 with comparison.open() as f:group_rows=list(csv.DictReader(f,delimiter='\t'))
 with (R/'validation_summary.md').open('a') as f:
  f.write('\n## Inside versus outside cohort comparison\n\n')
  f.write('The original arrangement-specific tracks remain limited to fully inside windows. For a direct regional comparison, the same 17 ancestral-associated and seven derived-associated samples per LG are also grouped across the flanking/outside SNP VCF. Those outside statistics are population-group summaries, not arrangement-linked backgrounds because recombination decouples them from the inversion. This comparison is descriptive and confounded with population structure/ecotype.\n\n')
  f.write('| LG | scope | valid windows | mean pi_X group | mean pi_Y group | mean dXY group | mean Hudson FST |\n|---|---|---:|---:|---:|---:|---:|\n')
  for x in group_rows:f.write(f"| {x['lg']} | {x['scope']} | {x['valid_windows']} | {float(x['mean_pi_X_group_SNP']):.5f} | {float(x['mean_pi_Y_group_SNP']):.5f} | {float(x['mean_dXY_group_SNP']):.5f} | {float(x['mean_FST_group_Hudson']):.5f} |\n")
  f.write('\nFull window values and interpretation labels are in inside_outside_cohort_group_windows.tsv; the separate inside arrangement-specific tracks retain NA outside.\n')
