# Anopheles 2La Tajima’s D analysis

**Status: VALIDATED for the within-group callable-site SFS; arrangement contrasts are descriptive and geographically confounded.** The run used the MalariaGEN Ag3.10 `fontaine-2015-rebuild` all-site genotype stream with AgamP4 reference and the `gamb_colu` site mask. All queried coordinates were checked against non-N bases in the reference, and returned site-mask values were checked against the official mask. The analysis covers the complete 2L arm, 1–49,364,325, and uses the existing 50-kb grid plus a 100-kb sensitivity grid. The 2La interval is 20,524,058–42,165,532 (1-based inclusive).

## Samples and calculation

There are 37 diploid individuals: 26 *A. gambiae* and 11 *A. coluzzii*. Groups include species pooled, standard homozygotes (A0), inverted homozygotes (A1), and locality groups with at least two individuals. Heterokaryotypes contribute only to pooled groups. Species arrangement strata span geography; no regional arrangement contrast is interpreted as a within-locality causal estimate. Twenty-two groups were analyzed. Each group-window has fixed sample size, requiring complete diploid genotypes at each included site. Positions must pass the existing mask; observed multiallelic sites are excluded from both SFS and callable denominator. Callable bases are therefore masked, non-N reference positions with complete group genotypes and biallelic or invariant state. Per-base pi and theta_W use that same group/window denominator. Tajima’s D uses the classical fixed-n S, pi, a_n and standard variance expression; n<=3 or undefined variance gives NA. No classical neutrality-test p-values are interpreted.

The regional summaries use 500-kb spatial block resampling (5,000 replicates) across 50-kb windows. The 50-kb table has 988 windows per group, 432 inside and 555 outside for most groups (the edge-window counts differ by one where appropriate). Spatial intervals describe uncertainty across linked sequence blocks, not independent-window inference.

## Regional results

| Group | D inside | D outside | Difference (95% 500-kb block interval) |
|---|---:|---:|---:|
| *A. gambiae* pooled (n=26) | −1.863 | −2.121 | +0.259 (+0.204, +0.319) |
| *A. gambiae* standard (n=8) | −1.292 | −1.544 | +0.251 (+0.196, +0.308) |
| *A. gambiae* inverted (n=12) | −1.674 | −1.677 | +0.003 (−0.059, +0.068) |
| *A. coluzzii* pooled (n=11) | −1.250 | −1.626 | +0.376 (+0.326, +0.431) |
| *A. coluzzii* standard (n=8) | −1.112 | −1.297 | +0.185 (+0.119, +0.252) |
| *A. coluzzii* inverted (n=3) | −0.655 | −0.649 | −0.006 (−0.077, +0.059) |

Per-callable-base diversity increases inside 2La. Pooled *A. gambiae* pi is 0.01717 inside vs 0.01150 outside (difference +0.00566; block interval +0.00471 to +0.00666); pooled *A. coluzzii* is 0.01655 vs 0.01151 (difference +0.00504; +0.00422 to +0.00591). Arrangement-stratified pi also rises inside: standard/inverted estimates are 0.01484/0.01401 vs 0.01262/0.01166 for *A. gambiae*, and 0.01239/0.01409 vs 0.01118/0.01218 for *A. coluzzii*. These are means of per-window estimates on validated callable denominators, not mutation-rate corrected quantities.

The pooled-minus-unweighted-background D contrast remains negative on both sides but is less negative inside: *A. gambiae* −0.380 vs −0.511 (inside-minus-outside +0.132; block interval +0.111 to +0.157); *A. coluzzii* −0.366 vs −0.654 (+0.287; +0.236 to +0.343). The sample-size-weighted contrast is +0.156 for *A. gambiae* and +0.244 for *A. coluzzii*. These contrasts use distinct species-wide pooled and homozygote sample panels and thus are descriptive, not pure within-population background effects.

## Integration and interpretation

Tajima D joins exactly by existing 50-kb `window_id` to callable pi, dXY/FST, quartet-frequency, branch/tree-diversity, and Patterson’s D tracks. In pooled species groups, D correlates positively with within-background pi (Spearman rho about 0.37–0.42) and between-arrangement dXY (rho 0.78 *A. gambiae*, 0.85 *A. coluzzii), and negatively with q_species (−0.54, −0.70) while positively with q_arrangement (+0.53, +0.70). Each estimate has a 500-kb block bootstrap interval in `anopheles_tajima_d_track_relationships.tsv`. These are descriptive associations among overlapping window statistics, not independent tests; Patterson’s D remains explicitly distinct from Tajima’s D. The existing gene-tree metrics are topology or substitution summaries, not calibrated genealogical-time estimates.

The combination of increased callable pi and shifted pooled D adds SFS shape information beyond pi/dXY/quartet tracks. It is compatible with altered linked genealogies under suppressed recombination but is not unique to MSRC. Geography, population structure, selection, demography, and arrangement ascertainment remain alternative explanations. Arrangement-sampled proportions are not natural population frequencies; the analysis does not use them as such. It also does not imply that Tajima’s D must be positive under MSRC.

Reproduce the all-site scan with `MPLBACKEND=Agg` and the repository’s `.venv-malaria310` using `run_anopheles.py`; if the sufficient-statistics checkpoint exists, `--summary-only` regenerates outputs without remote queries. Then run `calculate_background_contrasts.py` and `analyze_track_relationships.py`. Outputs and figures are in this dataset’s `analysis/tajima_d/` directory.
