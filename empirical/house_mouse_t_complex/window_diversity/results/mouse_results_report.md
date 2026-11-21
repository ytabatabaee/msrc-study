# House mouse t-complex diversity and differentiation (phase-aware SNP-panel analysis)

Date: 2026-11-20. This is a sitewise reconstructed t-allele analysis on a PASS SNP panel, not a per-callable-base estimate and not a chromosome-scale phased reconstruction. The later callable-site audit confirmed that public individual gVCFs exist, but indexed retrieval failed before a complete chr17 mask could be built; the results below remain exploratory SNP-conditioned estimates.

## Sources, samples, and coordinate system

The sequence source is Harr et al.'s publicly indexed AllMouse PASS VCF (67 samples), accessed by bounded range queries through `https://wwwuser.gwdguser.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz`; the matching local TBI was used. A direct chr17 query verified sample IDs, contig length (94,987,271 bp), and coordinate compatibility. This PASS VCF is variant-only and contains no invariant-call denominator. The tree-matched scan covered 4,046 non-overlapping 5-kb local-tree windows from mm10/GRCm38 chr17 5,000,000 through 37,055,000 and retained 1,044,870 biallelic SNP records. A subsequent indexed scan covered the entire VCF contig, [0, 94,987,271) in 0-based half-open coordinates, across 18,998 5-kb bins and 5,086,154 biallelic PASS SNP records. The full scan includes 37.055–42.05 Mb and extends the genotype tracks through the chromosome end; local gene trees still stop at 37.055 Mb.

The frozen sample design in `empirical/house_mouse_t_complex/analysis/patterson_d/config/frozen_quartet_members.tsv` defines separate within-subspecies comparisons: Mmd/domesticus (4 t carriers, 12 standard controls), Mmm/musculus (8, 13), and Mmc/castaneus (3, 7). No subspecies were pooled. Standard genotypes are represented by both diploid allele copies. Carriers are not treated as pure t samples: a homozygous carrier contributes its observed allele as the putative t allele; at a heterozygous carrier site, the t allele is inferred only if exactly one of its two alleles is absent from all same-subspecies standard chromosomes. Missing carrier calls and heterozygotes for which neither or both alleles are unique are masked. All standard calls must be present at the SNP; at least two reconstructed carrier alleles are required to estimate within-t diversity. This is a conservative, per-site adaptation of published pseudo-t reconstruction logic and does not phase alleles across loci.

Four inversion intervals are shaded using the existing validated GRCm38 coordinate table in `empirical/house_mouse_t_complex/results/t_complex_inversion_intervals.tsv`: In(17)1 5,747,331–6,957,919; In(17)2 7,435,308–13,016,510; In(17)3 14,068,229–18,051,507; In(17)4 19,401,000–40,049,964. The archived tree windows end at 37,055,000, so only the portion of In(17)4 within the tree-data span is analyzed.

## Estimators and missing data

At a biallelic SNP, `pi = n/(n-1) * (1 - sum(p_k^2))`, using `n` reconstructed haploid t alleles for pi_t and `2*n_standard` diploid standard chromosomes for pi_standard. `dXY = 1 - sum(p_t,k*p_standard,k)`. Hudson FST is `sum[dXY - (pi_t + pi_standard)/2] / sum[dXY]` over retained SNP records in each window. pi and dXY are arithmetic means across the same retained polymorphic records. The previous denominator was observed PASS SNP count because that was the site set exposed by the indexed input. This yields values per observed PASS SNP, not per callable base; it does not estimate missing positions as invariant.

Per-window mean reconstructed t-carrier counts are 3.91/4 for Mmd, 7.70/8 for Mmm, and 2.90/3 for Mmc. Among t-carrier call opportunities at sites where all standard calls were available, missing or phase-ambiguous carrier genotypes accounted for 4.9% (Mmd), 6.3% (Mmm), and 8.8% (Mmc). These excluded alleles are not imputed. No BAM-level DP/GQ filters were applied to the PASS VCF calls. The public `vcf/GVCF/` archive has individual reference-block gVCFs; a query confirmed reference-block GT, DP, GQ and MIN_DP fields and mm10 chr17 length 94,987,271. But a low-concurrency multi-sample range pilot failed on `CR46.all.vcf.gz` at BGZF offset 8,170,095,228 with a truncated-file read and close timeout. No complete mask was produced, so no callable-base estimate is reported. See `results/mouse_callable_site_audit.md` for the input inventory and recovery attempt. Geographic sampling varies within the Mmd and Mmm strata, so within-subspecies estimates still carry residual population structure.

Inside/outside contrasts use a 500-kb physical-block percentile bootstrap with 10,000 resamples. Windows overlapping any of the four intervals are inside; all other retained tree windows are outside. Each contrast contains 3,824 inside and 222 outside windows; one Mmm FST window has a zero summed dXY denominator and is undefined, so its FST contrast uses 3,823 inside windows. The block resampling uses 56 inside versus only 4 outside blocks. Intervals are broad and descriptive; no independent-window tests are used.

## Numerical results

All values are conditional per observed SNP. Difference is inside minus outside.

| Subspecies design | Statistic | Inside | Outside | Difference | 95% 500-kb block interval |
|---|---|---:|---:|---:|---:|
| Mmd/domesticus | Hudson FST | 0.5570 | 0.3557 | +0.2013 | [−0.3176, +0.4936] |
|  | dXY | 0.0993 | 0.0812 | +0.0181 | [−0.0814, +0.0681] |
|  | pi t pseudo-haplotype | 0.0238 | 0.0169 | +0.0068 | [−0.0029, +0.0152] |
|  | pi standard | 0.0328 | 0.0356 | −0.0027 | [−0.0220, +0.0130] |
| Mmm/musculus | Hudson FST | 0.6821 | 0.3192 | +0.3628 | [−0.2039, +0.6557] |
|  | dXY | 0.1328 | 0.0885 | +0.0443 | [−0.0421, +0.0841] |
|  | pi t pseudo-haplotype | 0.0235 | 0.0313 | −0.0078 | [−0.0175, +0.0090] |
|  | pi standard | 0.0419 | 0.0422 | −0.0003 | [−0.0123, +0.0182] |
| Mmc/castaneus | Hudson FST | 0.5413 | 0.2857 | +0.2556 | [−0.1478, +0.4804] |
|  | dXY | 0.1366 | 0.1250 | +0.0116 | [−0.0390, +0.0352] |
|  | pi t pseudo-haplotype | 0.0263 | 0.0630 | −0.0367 | [−0.0587, +0.0052] |
|  | pi standard | 0.0854 | 0.0957 | −0.0102 | [−0.0274, +0.0129] |

Mean FST and dXY are directionally higher inside the intervals in all three strata, a direction compatible with reduced exchange between t and standard backgrounds. However, all block intervals include zero, and the outside contrast uses only four physical blocks. The pattern is not precise evidence for an inside/outside effect. The pi_t direction is mixed: higher inside in Mmd, lower inside in Mmm and Mmc. No simple universal pi prediction is supported by these SNP-panel summaries.

## Relationship to existing genealogy tracks

Correlations are descriptive and use the same local windows (4,046 except the Mmm FST row, which has one undefined window). Hudson FST has Pearson r=0.507, 0.465, and 0.577 with existing A_distance in Mmd, Mmm, and Mmc, respectively; dXY has r=0.684, 0.674, and 0.697. The B-like `mean_pair_distance_no_spretus` relation is weak or negative (FST r=−0.214, −0.025, −0.014; dXY r=−0.202, −0.064, −0.115). Correlations with quartet `q_t_alt` and Patterson's D are near zero for the diversity measures (absolute Pearson r≤0.082 for FST/D and ≤0.047 for dXY/D). Per-design Pearson and Spearman summaries are in `results/mouse_t_pseudohap_tree_track_correlations.tsv`.

## Outputs and interpretation

- Reproducible range-query and reconstruction script: `scripts/analyze_t_pseudohap_snp_panel.py`.
- Per-window statistics plus A/B/quartet tracks: `results/mouse_t_pseudohap_snp_panel_window_stats.tsv`.
- Inside/outside block-bootstrap effects: `results/mouse_t_pseudohap_inside_outside.tsv`.
- Track correlations and source/reconstruction audit: `results/mouse_t_pseudohap_tree_track_correlations.tsv`, `results/mouse_source_and_reconstruction_audit.json`.
- Three separate Mmd/Mmm/Mmc four-track figures: `figures/D1_Mmd_t_vs_standard_snp_panel_tracks.{pdf,png}`, `D2_Mmm_t_vs_standard_snp_panel_tracks.{pdf,png}`, and `D3_Mmc_t_vs_standard_snp_panel_tracks.{pdf,png}`.
- `results/mouse_scan_checkpoint.json` retains completed range-block accumulators for transparent resumption/reproduction.
- Full chr17 5-kb window statistics and inversion-union contrast: `results/mouse_t_pseudohap_chr17_full_window_stats.tsv`, `results/mouse_t_pseudohap_chr17_full_summary.tsv`; chromosome-wide four-track figure: `figures/mouse_t_pseudohap_chr17_full_tracks.{pdf,png}`. `has_local_tree` flags the 4,046 matched 5-kb windows. The run checkpoint is `results/mouse_full_chr17_checkpoint_5Mb.json` (completed-range metadata includes the smaller resumed ranges after the original remote 5-Mb query truncated around 75 Mb).
- Full-chromosome reproduction script: `scripts/analyze_full_chr17_snp_panel.py`.
- Inversion-inside versus outside average bars with 95% 500-kb block-bootstrap intervals: `figures/mouse_t_pseudohap_chr17_inversion_barplots.{pdf,png}`; plotted means and intervals are in `results/mouse_t_pseudohap_chr17_inversion_barplot_summary.tsv`, generated by `scripts/plot_chr17_inversion_barplots.py`.
- Exact current-input, local-file, public-gVCF, and failed-retrieval audit: `results/mouse_callable_site_audit.md`. The full-chromosome tracks and bars remain explicitly SNP-conditioned exploratory outputs.

## Distal genotype windows beyond the tree tracks

The 37.055–42.050 Mb extension uses the same within-subspecies designs, sitewise carrier reconstruction, Hudson estimator and SNP-only denominators described above. It compares the distal In(17)4 portion (37.055–40.049964 Mb) with the immediately right flank (40.050–42.050 Mb), using 500-kb physical-block bootstrap intervals. The interval provides 7 inversion blocks and 5 flank blocks; its small number of spatial units makes the contrast exploratory. Estimates are means over windows with observed SNPs and valid estimates.

This narrower distal contrast was the initial extension beyond the tree data; it is now a subset of the full-chromosome scan below. Use the full-chromosome figure and table for current visualization and chromosome-wide summaries.

| Design | Measure | Distal In(17)4 | Right flank | Difference | 95% 500-kb block interval |
|---|---|---:|---:|---:|---:|
| Mmd/domesticus | Hudson FST | 0.7172 | 0.0781 | +0.6391 | [+0.2869, +0.6848] |
|  | dXY per observed SNP | 0.1509 | 0.0447 | +0.1062 | [+0.0045, +0.1133] |
|  | pi t pseudo-haplotype | 0.0215 | 0.0309 | −0.0094 | [−0.0527, −0.0013] |
|  | pi standard | 0.0307 | 0.0471 | −0.0164 | [−0.0620, −0.0080] |
| Mmm/musculus | Hudson FST | 0.5210 | 0.0970 | +0.4240 | [+0.1771, +0.5641] |
|  | dXY per observed SNP | 0.1679 | 0.0570 | +0.1108 | [+0.0710, +0.1331] |
|  | pi t pseudo-haplotype | 0.0377 | 0.0334 | +0.0044 | [−0.0172, +0.0374] |
|  | pi standard | 0.0950 | 0.0641 | +0.0309 | [+0.0069, +0.0572] |
| Mmc/castaneus | Hudson FST | 0.6183 | 0.1072 | +0.5111 | [+0.3568, +0.5988] |
|  | dXY per observed SNP | 0.1609 | 0.0964 | +0.0646 | [+0.0231, +0.0894] |
|  | pi t pseudo-haplotype | 0.0297 | 0.0776 | −0.0479 | [−0.0605, −0.0272] |
|  | pi standard | 0.0694 | 0.0925 | −0.0231 | [−0.0443, −0.0069] |

All three strata show higher SNP-panel FST and dXY in distal In(17)4 than in the right flank, directionally consistent with differentiation across the rearranged background; pi responses vary by group. These are not per-callable-base estimates and the comparison is not a full inversion-versus-genome contrast. The apparent precision of the block intervals should be read in light of only five right-flank blocks, the SNP ascertainment, and unresolved geographic/read-level effects. These distal files are retained as the initial subset: `results/mouse_t_pseudohap_distal_flank_window_stats.tsv`, `results/mouse_t_pseudohap_distal_flank_summary.tsv`, and `figures/mouse_t_pseudohap_distal_flank_tracks.{pdf,png}`.

## Full-chromosome genotype tracks

The all-sites-position indexed query completed from chr17 0 through the VCF contig end (94,987,271 bp) in 5-kb windows. The count is 5,086,154 biallelic A/C/G/T PASS SNP records. Each subspecies design contains 18,998 bins; valid estimates occur in 18,044 Mmd, 18,019 Mmm and 18,150 Mmc bins. The figure shades the four validated GRCm38 inversion intervals and marks the local-tree span limits (5.0 and 37.055 Mb) with dashed lines. `has_local_tree` in the per-window table explicitly distinguishes the 4,046 matched windows per design from all added genotype-only windows. No A/B/quartet correlations are inferred outside the tree span.

For a chromosome-scale summary, windows are classified by their midpoint as inside the union of In(17)1–4 or outside that union; each 500-kb genomic block is the resampling unit (10,000 bootstrap replicates). Values below are means of valid 5-kb window estimates, and differences are inversion union minus all outside windows. FST is the windowwise Hudson ratio; pi and dXY remain averages over retained observed SNPs, not invariant-inclusive per-base measures.

| Subspecies | Measure | Inversion union | Outside | Difference | 95% 500-kb block interval |
|---|---|---:|---:|---:|---:|
| Mmd/domesticus | Hudson FST | 0.5551 | 0.0734 | +0.4817 | [+0.4106, +0.5374] |
|  | dXY per observed SNP | 0.1003 | 0.0451 | +0.0552 | [+0.0409, +0.0709] |
|  | pi t pseudo-haplotype | 0.0223 | 0.0350 | −0.0127 | [−0.0182, −0.0073] |
|  | pi standard | 0.0346 | 0.0404 | −0.0058 | [−0.0112, −0.0015] |
| Mmm/musculus | Hudson FST | 0.6284 | 0.0576 | +0.5708 | [+0.4900, +0.6075] |
|  | dXY per observed SNP | 0.1321 | 0.0452 | +0.0869 | [+0.0744, +0.0974] |
|  | pi t pseudo-haplotype | 0.0262 | 0.0355 | −0.0093 | [−0.0126, −0.0021] |
|  | pi standard | 0.0491 | 0.0438 | +0.0053 | [−0.0007, +0.0129] |
| Mmc/castaneus | Hudson FST | 0.5307 | 0.0906 | +0.4401 | [+0.4016, +0.4669] |
|  | dXY per observed SNP | 0.1347 | 0.1040 | +0.0308 | [+0.0216, +0.0410] |
|  | pi t pseudo-haplotype | 0.0279 | 0.0811 | −0.0532 | [−0.0547, −0.0486] |
|  | pi standard | 0.0814 | 0.1028 | −0.0214 | [−0.0261, −0.0169] |

There are 68 inversion and 124 outside 500-kb blocks. FST and dXY are higher in inversion windows for all three within-subspecies designs; reconstructed t diversity is lower in all three, while standard diversity is mixed. This broad pattern is compatible with persistent differentiation across suppressed-recombination backgrounds. However, the outside group pools all other chr17 positions, the panel is SNP-ascertained, local-tree correspondence stops at 37.055 Mb, and the estimates are not per callable nucleotide. These chromosome-wide contrasts are descriptive support for an MSRC-consistent pattern, not a mechanistic test. Existing tree-track correlations remain restricted to the matched 5–37.055 Mb span.

For a direct visual comparison of average levels, the companion four-panel bar figure plots the arithmetic mean of valid 5-kb window estimates inside and outside the inversion union for each subspecies. Error bars are 95% percentile intervals from 10,000 resamples of 500-kb block means; they describe spatial uncertainty, not independent-window standard errors. Inversion windows show substantially higher FST and dXY in every subspecies. πt is lower inside in all three; πstandard is lower for Mmd and Mmc and slightly higher for Mmm. Because the bars inherit SNP-only ascertainment and the outside group pools the rest of chr17, interpret them as descriptive chromosome-wide contrasts.

The mean differentiation tracks are directionally compatible with t-versus-standard divergence persisting across the rearranged region, and FST/dXY covary strongly with the existing A_distance genealogy contrast within the tree-matched span. Yet uncertainty in that broad inside/outside analysis is large, the source panel is SNP-ascertained, t phase is reconstructed site by site, and geographic structure persists within subspecies. The distal extension has no local genealogy tracks. These results are exploratory MSRC-consistent directions, not callable-base estimates or conclusive tests of the mechanism. The next mouse step is to map the individual gVCF reference blocks for the exact frozen samples, apply a common callability/depth rule, recompute per-callable-base windows, and test geographic sensitivity.
