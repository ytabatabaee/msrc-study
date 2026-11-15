# Exact matched-quartet comparison

For each gene tree, the quartet scorer used the exact frozen sample IDs. A four-point distance test classified each induced quartet once as q1=(P1,P2)|(P3,O), q2=(P1,P3)|(P2,O), or q3=(P2,P3)|(P1,O); tied sums are unresolved. q frequencies are normalized over resolved induced quartets within a tree. Across loci, each tree receives equal weight; induced quartets are never treated as independent loci. D_topology=q3-q2. D_site uses fixed-REF plus fixed-ALT ancestral polarization and pooled ABBA/BABA counts.

## Exact topology frequencies and correspondence

- D1_Mmd_t_vs_standard: equal-locus q1=0.6514, q2=0.1350, q3=0.2136; D_topology=+0.0785; 20,133,386 summed induced quartet scores are diagnostics only across 4,043 loci.
- D2_Mmm_t_vs_standard: equal-locus q1=0.3357, q2=0.4128, q3=0.2515; D_topology=-0.1613; 40,251,014 summed induced quartet scores are diagnostics only across 4,042 loci.
- D3_Mmc_t_vs_standard: equal-locus q1=0.2668, q2=0.3313, q3=0.4019; D_topology=+0.0705; 8,816,571 summed induced quartet scores are diagnostics only across 4,045 loci.

## Spatial correspondence

Pearson r (Spearman rho) for corrected D_site versus D_topology at each pooled resolution; CIs resample matched 500-kb physical blocks and the circular-shift null rotates paired 500-kb block summaries. These are descriptive single-chromosome spatial tests, not independent-locus p-values.

| Design | Resolution | windows | Pearson r (95% block CI) | Spearman rho | sign agreement | circular-shift p |
|---|---:|---:|---:|---:|---:|---:|
| D1_Mmd_t_vs_standard | 50 kb | 467 | +0.615 [+0.490, +0.713] | +0.511 | 0.585 | 0.017 |
| D1_Mmd_t_vs_standard | 100 kb | 248 | +0.688 [+0.560, +0.776] | +0.622 | 0.669 | 0.017 |
| D1_Mmd_t_vs_standard | 250 kb | 108 | +0.768 [+0.673, +0.839] | +0.720 | 0.694 | 0.017 |
| D1_Mmd_t_vs_standard | 500 kb | 58 | +0.806 [+0.714, +0.879] | +0.784 | 0.707 | 0.017 |
| D2_Mmm_t_vs_standard | 50 kb | 467 | +0.587 [+0.481, +0.677] | +0.536 | 0.567 | 0.017 |
| D2_Mmm_t_vs_standard | 100 kb | 248 | +0.629 [+0.507, +0.723] | +0.564 | 0.617 | 0.017 |
| D2_Mmm_t_vs_standard | 250 kb | 108 | +0.679 [+0.558, +0.773] | +0.641 | 0.611 | 0.017 |
| D2_Mmm_t_vs_standard | 500 kb | 58 | +0.664 [+0.500, +0.789] | +0.628 | 0.707 | 0.017 |
| D3_Mmc_t_vs_standard | 50 kb | 467 | +0.641 [+0.558, +0.710] | +0.591 | 0.675 | 0.017 |
| D3_Mmc_t_vs_standard | 100 kb | 248 | +0.711 [+0.619, +0.778] | +0.646 | 0.673 | 0.017 |
| D3_Mmc_t_vs_standard | 250 kb | 108 | +0.776 [+0.670, +0.855] | +0.726 | 0.713 | 0.017 |
| D3_Mmc_t_vs_standard | 500 kb | 58 | +0.798 [+0.672, +0.897] | +0.773 | 0.724 | 0.017 |

At 100 kb, the within-inversion correlations are D1 r=+0.674, D2 r=+0.620, D3 r=+0.707. Only 15 100-kb bins are wholly outside the inversion union in the matched tree span; outside-only correlations are therefore imprecise and have too few 500-kb blocks for a circular-shift null. At the same resolution, the exact D1 correlation is positive (+0.688), unlike the prior r=-0.097 against q_t_alt. Exact q3−q2 matching removes that weak negative result, but the change also reflects the different full-tip q_t_alt estimand and corrected polarization, so the cause cannot be assigned to quartet definition alone.


At the fundamental 5-kb locus unit, split at each design’s median induced internal branch length, the D_site–D_topology Pearson correlation is short/long: D1 +0.303/+0.404, D2 +0.370/+0.535, and D3 +0.331/+0.565. Concordance is stronger for long-branch quartets, but remains positive in the short-branch half. Exact estimates are in `branch_length_sensitivity.tsv`.

## Matched inversion contrast

| Design | D_site inside | D_site matched complement | delta D_site | D_topology inside | D_topology matched complement | delta D_topology | 500-kb block CI for delta D_site | 500-kb block CI for delta D_topology |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| D1_Mmd_t_vs_standard | +0.353 | +0.463 | -0.110 | +0.075 | +0.142 | -0.067 | [-0.473, +0.344] | [-0.336, +0.073] |
| D2_Mmm_t_vs_standard | +0.197 | +0.347 | -0.150 | -0.175 | +0.068 | -0.242 | [-0.614, +0.252] | [-0.618, -0.066] |
| D3_Mmc_t_vs_standard | +0.242 | +0.216 | +0.025 | +0.066 | +0.154 | -0.088 | [-0.495, +0.218] | [-0.421, +0.085] |

Concordance is interval-specific. Both statistics rise together inside In(17)1 for all designs and generally also in In(17)2; In(17)3 is mixed, and In(17)4 shows lower D_site inside than in its matched complement for all three designs. At the four-inversion union, D1 and D2 are lower inside than matched collinear windows for both statistics; D3 has a small positive D_site difference but lower D_topology. See `matched_inside_outside.tsv` for all individual inversions, flanks, and gaps. The small collinear gap-1 coverage (one resolved 5-kb locus per design) is not interpretable as a stable interval estimate.

The previous D1 r=-0.097 was computed against q_t_alt, a different full-tip topology statistic. The exact q3−q2 contrast is mathematically matched to ABBA−BABA, and D1 is positive across tested resolutions (+0.615 to +0.806). This is evidence that the old negative/near-zero association is not robust to the corrected estimand; it is not proof that sample grouping alone caused the discrepancy.

## Limitations

The carriers in the archived trees are reconstructed pseudo-t-haplotype sequence tips generated from t-specific SNP subtraction, whereas VCF D uses unphased diploid genotypes containing both homologs. Thus the individual labels match exactly but chromosome-level observations are not identical. Population and geographic structure, cross-subspecies history, introgression/admixture, ancestral structure, ILS, and SNP ascertainment can all create correspondence or discordance. Local trees stop at chr17:37,055,000; all nucleotide/tree comparisons use exactly matched 5-kb starts within that range. Full-contig nucleotide summaries are separate. The 8 M. spretus samples are the only designated outgroup in the frozen design/source cohort; no alternative outgroup was selected.
