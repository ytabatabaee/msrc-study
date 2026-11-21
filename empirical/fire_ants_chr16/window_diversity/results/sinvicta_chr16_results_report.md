# Fire-ant SB/Sb diversity and differentiation (SNP-panel analysis)

Date: 2026-11-20. This preliminary analysis is restricted to *Solenopsis invicta* males from Argentina in the existing local Stolle-cohort VCF. It is a variant-panel analysis; the VCF has no invariant records or validated callable-base denominator, so values are per retained observed SNP, not per callable nucleotide.

## Data, sample selection, and coordinate matching

The genotype source is `data/fire_ants_chr16/raw/genotypes/gt.vcf.gz` (VCF reference label `gng20170922wFex.fa`, FORMAT GT and DP; 368 sample columns in its header). Direct IDs and assignments were intersected across `empirical/fire_ants_chr16/analysis/patterson_d/config/frozen_sample_membership.tsv` and `results/vcf_sample_inventory.tsv`, retaining only Q1 primary-Argentina *S. invicta* males: 78 SB and 30 Sb. All retained individuals are reported males. Calls were treated as haploid only when the diploid-encoded GT was homozygous; heterozygous and missing calls were excluded at a site. Per-sample global missingness is already audited in the source inventory.

The 96 chr16 local-gene-tree windows (window indices 118–213) are joined through the existing exact BUSCO-gene mapping in `analysis/patterson_d/results/window_gene_coordinate_mapping.tsv`. That table links each gene's coordinates on the VCF contig to a local tree window; estimates aggregate observed VCF SNPs over the mapped gene-body intervals, not all bases in the physical tree-window envelope. The local tree axis uses the Stolle / `gng20170922wFex` coordinates. The upstream pipeline and local workflow README identify `chr16nr` as the nonrecombining supergene and `chr16r` as recombining chromosome-16 background. The mapped comparison contains 52 `chr16nr` windows and 44 `chr16r` windows. The author-defined BUSCO-window supergene envelope (chr16:11,680,438–27,917,498) is shaded and bounded in the tracks; this is not a new whole-chromosome assembly lift or an exact inversion breakpoint set.

## Estimators and filters

Only biallelic A/C/G/T SNPs in mapped BUSCO intervals were considered. A haploid SB/Sb sample call is retained when GT is homozygous and DP≥3; each group must have at least 80% of its samples called at that SNP. The mean called sample counts are provided per window. The current VCF lacks GQ and read alignments, so no mapping-quality or allele-balance check is possible.

For the retained haploid sample alleles, within-group diversity is unbiased gene diversity `pi = n/(n-1) * (1 - sum(p_k^2))`; `dXY = 1 - sum(p_SB,k * p_Sb,k)`; Hudson FST is the ratio of summed sitewise `dXY - (pi_SB + pi_Sb)/2` to summed sitewise dXY. Window pi and dXY are arithmetic means across retained observed variable SNPs. They are conditional on the source's SNP ascertainment and call filter, never per-base values. Inside/outside intervals are 95% percentile intervals from 10,000 resamples of 1-Mb blocks after averaging windows within blocks. The contrasts use 13 inside and 9 outside physical 1-Mb blocks; the resulting percentile intervals are spatial descriptive intervals, not formal hypothesis tests.

## Results

The relaxed, explicit site rule retains 141,185 SNP observations across all 96 windows (all windows have at least one); effective mean calls per SNP are about 76/78 SB and 29/30 Sb. The corrected suppressed-recombination group contains 52 `chr16nr` windows; the recombining background contains 44 `chr16r` windows.

| Statistic (per observed SNP) | chr16nr supergene mean | chr16r background mean | Difference (inside − outside) | 95% 1-Mb block interval |
|---|---:|---:|---:|---:|
| Hudson FST | 0.54264 | 0.00062 | +0.54202 | [+0.41136, +0.59319] |
| dXY | 0.05299 | 0.03092 | +0.02207 | [+0.01648, +0.03191] |
| pi SB | 0.02740 | 0.03102 | −0.00362 | [−0.00472, +0.00469] |
| pi Sb | 0.01742 | 0.03076 | −0.01334 | [−0.01430, −0.00510] |

The supergene has higher SB/Sb differentiation and lower Sb diversity in this SNP panel, directionally consistent with a simple MSRC prediction. The pi_SB interval includes zero. The block bootstrap uses 13 inside and 9 outside 1-Mb blocks; intervals are descriptive rather than formal hypothesis tests. These remain local BUSCO-gene, SNP-ascertained comparisons without callable-base lengths or read-level quality validation, and should not be generalized to genome-wide nucleotide diversity.

Within `chr16nr`, dXY has modest correlation with the B-like mean focal-pair distance (Pearson r=0.166, rho=0.329) and with A_distance (r=0.224, rho=0.546). Hudson FST correlates with A_distance (r=0.198, rho=0.476) and q_H (r=0.645, rho=0.395); q_3 associations are weak. Full correlations with A-distance, A-relative, B-like mean focal pair distances, and quartet-support frequencies are retained in the track tables and are descriptive. The first report draft reversed the `chr16nr`/`chr16r` interpretation; the upstream README identifies `chr16nr` as the nonrecombining supergene, and all region contrasts and boundary shading here have been corrected from the saved window estimates.

## Outputs

- Script: `scripts/analyze_sinvicta_snp_panel.py`.
- Window table and effective per-window sample counts: `results/sinvicta_chr16_snp_panel_window_stats.tsv`.
- Spatial contrasts: `results/sinvicta_chr16_inside_outside.tsv`.
- Join to existing genealogy and quartet tracks: `results/sinvicta_chr16_tree_diversity_join.tsv`.
- Correlations: `results/sinvicta_chr16_tree_track_correlations.tsv`.
- Track figure with the mapped chr16nr supergene envelope: `figures/sinvicta_chr16_snp_panel_tracks.pdf` and `.png`.

The requested Dryad 16-male whole-genome variant VCF from [Dryad 10.5061/dryad.js509](https://datadryad.org/dataset/doi:10.5061/dryad.js509) was not retrieved because direct Dryad file-stream downloads returned HTTP 403 in this environment. The local genotype cohort supplies an alternative and supports a transparent SNP-panel analysis, but it does not replace all-sites callability or the published independent 8+8 male sample design.
