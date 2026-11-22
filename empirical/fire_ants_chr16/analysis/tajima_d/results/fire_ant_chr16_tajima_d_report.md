# Fire-ant chr16 Tajima's D — exploratory observed-panel analysis

**Status: EXPLORATORY.** This analysis uses the existing FreeBayes variant-only VCF and is not a validated unbiased site-frequency spectrum. The callset's completeness for all segregating variants in the sampled males is not established; missing invariant sites alone would not invalidate D, but filtered or uncalled segregating sites can bias it. The VCF covers selected BUSCO gene intervals mapped to 96 local chr16 gene-tree windows, not continuous callable chromosome sequence.

## Samples, ploidy, and coordinates

The frozen Q1 primary-Argentina sample membership intersected with the VCF inventory yields 78 SB and 30 Sb *S. invicta* males. The calculation evaluates SB (n=78 haplotypes), Sb (n=30), and pooled (n=108) groups. Male diploid-encoded GTs are treated as haploid only when homozygous; heterozygous calls are excluded, and each included genotype must have DP>=3. Statistics use only sites completely called for every chromosome in the named group; therefore n is fixed within each window. The 50-kb/100-kb concept is not applicable to these noncontiguous gene-body target intervals; results follow existing mapped gene-tree windows, local chr16 coordinate (Stolle/gng20170922wFex axis), and `chr16nr` versus `chr16r` annotations. The author-defined supergene envelope is chr16:11,680,438–27,917,498 in that local coordinate system; no chromosome-scale assembly transform is assumed.

The scan traversed 146,149 records in the mapped target intervals and produced 96 window rows. Summed across window assignments, segregating-site counts S were 3,760 pooled, 4,570 SB, and 6,776 Sb. The per-window output separately records all completely called biallelic VCF records examined for each group (including records monomorphic within that group) and the subset segregating in that group. Monomorphic records contribute zero to pi and do not enter S, so omitting them from the classical D vector does not change D. Each group had defined D in all 96 windows. Neither record count is a callable-base denominator; per-base pi/theta are unavailable. Across candidate biallelic records in the mapped target intervals, missing/partial GT fractions were 1.0147% (SB) and 0.6921% (Sb); DP<3 fractions were 1.9342% and 1.2331%, respectively. Heterozygous male calls were absent (0%). Heterozygous GT fractions in these mapped target records were 0. The per-window table retains these audit fractions, as well as counts of excluded non-biallelic sites. These missingness values are calculated over candidate VCF genotypes in target intervals; the complete-case Tajima estimates use only sites complete across the entire named group.

## Results

In the 52 `chr16nr` windows versus 44 `chr16r` windows, mean observed-panel D was:

| Statistic | Inside | Outside | Inside − outside | 1-Mb block interval |
|---|---:|---:|---:|---:|
| Pooled (n=108) | -2.125 | -2.119 | -0.006 | -0.185 to +0.073 |
| SB (n=78) | -2.028 | -2.139 | +0.111 | +0.001 to +0.228 |
| Sb (n=30) | -1.223 | -1.677 | +0.453 | +0.194 to +0.565 |
| Matched SB (shared sites) | -1.999 | -2.081 | +0.082 | -0.094 to +0.141 |
| Matched Sb (shared sites) | -1.063 | -1.856 | +0.793 | +0.341 to +0.954 |
| Matched pooled − unweighted background mean | -0.594 | -0.150 | -0.443 | -0.495 to -0.271 |
| Matched pooled − sample-size-weighted background mean | -0.386 | -0.100 | -0.285 | -0.305 to -0.186 |

The pooled-minus-background contrast is also computed on the **same sites complete across all 108 males** in each window (mean 203.3 matched variant records per window; 19,518 window-record assignments total). On that matched set, SB D averaged -1.999 inside vs -2.081 outside; Sb D averaged -1.063 vs -1.856; the unweighted pooled-minus-background contrast averaged -0.594 vs -0.150 (inside-minus-outside difference -0.443; 1-Mb block interval -0.495 to -0.271). The sample-size-weighted matched contrast was -0.386 inside vs -0.100 outside (difference -0.285; interval -0.305 to -0.186). The block intervals resample 1-Mb blocks along the local tree coordinate; they are descriptive spatial uncertainty intervals, not neutrality-test p-values. Adjacent windows and mapped genes are linked, and the local annotation map includes overlapping BUSCO intervals.

## Interpretation and limits

Pooling does **not** elevate D above the within-background average in this panel: the pooled-minus-background contrast is negative inside and outside, with a more negative mean in the supergene. This contrasts with the Neurospora pilot's positive pooled/background contrast. The pattern does not contradict MSRC by itself: Tajima's D depends on the full sampled genealogy and frequency spectrum, and its sign is not fixed by suppressed recombination. Population mixture, demography, selection, source filters, and severe complete-case ascertainment are viable alternatives. The current estimate is useful as an exploratory description only; no mechanism or neutrality-test significance is claimed.

The joined table maps these windows to the existing FST, dXY, pi, quartet, and tree-diversity tracks using exact `window_index`; correlations are descriptive and not independent replications. The common-site pooled/background contrast avoids comparing group-specific call masks, but remains conditional on the variant-only, gene-ascertained source. No Patterson's D value is relabeled as Tajima's D.

Reproduce from repository root with `python3 empirical/fire_ants_chr16/analysis/tajima_d/scripts/run_fire_ants.py` followed by `python3 empirical/fire_ants_chr16/analysis/tajima_d/scripts/summarize_fire_ants.py`. Outputs: `fire_ant_chr16_tajima_d_windows.tsv`, `fire_ant_chr16_tajima_d_inside_outside.tsv`, joined tracks, JSON audit, and PDF/PNG figure in this directory. The regional table distinguishes S records used in D from all complete biallelic records examined (including group-monomorphic records); the latter still do not form a callable-genome denominator. Reproduce by running both scripts in order.
