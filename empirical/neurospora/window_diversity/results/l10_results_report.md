# Neurospora L10 diversity and differentiation (SNP-panel analysis)

Date: 2026-11-20. This is a preliminary SNP-ascertained analysis, not an estimate per callable nucleotide. The local genotype matrix contains filtered variable sites only; invariant positions and a validated callable-base mask are absent. The requested per-base quantities remain blocked until the published all-sites VCFs or equivalent callable sequence are available.

## Sources, samples, and coordinates

Input calls are from `data/neurospora/raw/neurospora/geno/cras_sit_hist_tet_L10.min20MAC1.geno.gz`, the same L10 panel used for the local tree windows in `empirical/neurospora/results/gene_tree_inventory.tsv`. The 22 CJ columns are haploid base calls. The mating-type audit in `empirical/neurospora/analysis/patterson_d/results/cj_mating_type_resolution.tsv` supports 8 mat-A calls (CJ01, CJ03, CJ05, CJ09, CJ13, CJ15, CJ19, CJ23) and 10 mat-a calls (CJ02, CJ04, CJ06, CJ08, CJ10, CJ14, CJ16, CJ20, CJ22, CJ24); four isolates (CJ07, CJ11, CJ12, CJ18) remain unresolved and were excluded. This is a verified subset, not the source's full 10-vs-12 sample design. No diploid or heterokaryon samples were treated as haploid population chromosomes.

The published Corcoran et al. all-sites L10 VCFs and consensus sequences are listed in Dryad dataset [10.5061/dryad.162mh](https://datadryad.org/dataset/doi:10.5061/dryad.162mh). Dataset metadata describes all-sites scaffold VCF archives and 2509 mat-a-reference consensus coordinates with repeats set to N. The archives were not locally available and direct file-stream retrieval returned HTTP 403 in this environment. Thus the current run uses local *Neurospora* linkage-group coordinates from the matched genotype/tree-window inventory, not the Dryad 2509 assembly. It does not assert a 2509-to-LGI lift.

For 44,986 of 45,303 tree windows, at least one matrix position in the recorded window coordinates had A/C/G/T calls in every one of the 18 included samples. Windows without any such complete panel site are omitted (317 windows). Each retained window is the existing 50-variable-site tree-window span; the estimate includes only complete, observed variable SNPs in that span. It does not include invariant callable positions. For the descriptive LGI region contrast, the approximate suppressed-recombination interval is 1.5–8.0 Mb; boundaries follow the existing Neurospora analysis convention and are marked on the tracks.

## Statistics

At each retained SNP, group frequencies are computed over the 8 mat-A and 10 mat-a haploid calls. Within-group gene diversity uses the unbiased haploid-sample correction, `pi = n/(n-1) * (1 - sum(p_k^2))`. `dXY = 1 - sum(p_A,k * p_a,k)`. The reported Hudson SNP-panel ratio is `sum[dXY - (pi_A + pi_a)/2] / sum[dXY]` across retained variable sites in a window. The window `pi` and `dXY` tracks are arithmetic means across those same complete observed SNPs. They are all **per observed variable SNP**, conditional on the source's global filtering/ascertainment, and must not be interpreted as per-base mutation/diversity rates. Negative or undefined FST values are retained as estimator outputs where applicable; zero denominator gives missing.

Inside/outside effect intervals are percentile 95% intervals from 10,000 independent resamples of 100-kb LGI blocks, with windows averaged within each block before resampling. These are descriptive block-bootstrap intervals; no independent-window p-value is reported. Track correlations are descriptive Pearson and Spearman correlations against existing B(w) (mean of the two matched-pair six-distance estimates) and the published q3−q1 quartet-frequency contrast, joined by window ID. The pre-existing Neurospora A(w) in manuscript terminology is arrangement-informed; it is not estimable for these mating-type groups and is not present as a matching track. The available pairwise branch contrast is also not A(w).

## Numerical summary

| Statistic (per observed variable SNP) | Inside mean | Outside LGI mean | Difference (inside − outside) | 95% 100-kb block-bootstrap interval |
|---|---:|---:|---:|---:|
| Hudson FST | 0.9950 | 0.1897 | +0.8053 | [+0.7184, +0.9875] |
| dXY | 0.3427 | 0.1101 | +0.2326 | [+0.2018, +0.2742] |
| pi mat-A | 0.00136 | 0.05910 | −0.05774 | [−0.06722, −0.03834] |
| pi mat-a | 0.00116 | 0.06010 | −0.05894 | [−0.06922, −0.03906] |

Each contrast uses 6,634 inside and 3,480 outside LGI windows with estimates; the bootstrap resamples 65 inside and 30 outside 100-kb blocks. The large within-region differentiation and very low within-group SNP diversity are consistent with strong allele-frequency separation in the ascertained SNP panel. They may also reflect ascertainment, missing-sample completeness filtering, and the known mating-type/introgression history; these values do not establish absolute nucleotide diversity or the magnitude of per-base differentiation.

Correlations with local genealogy are mixed and descriptive. Within the interval, SNP-panel dXY has Pearson r=0.0169 and Spearman rho=−0.0124 with B(w), and r=0.1563 / rho=0.1580 with q3−q1. Hudson FST has r=0.0494 / rho=0.1160 with B(w), and r=0.0856 / rho=0.0352 with q3−q1. The low/weak correlations mean the topology shift and SNP-panel differentiation do not track one another window-by-window in a simple linear fashion. The equivalent per-window joins and all correlations are retained in `results/l10_tree_diversity_join.tsv` and `results/l10_tree_track_correlations.tsv`.

## Outputs and interpretation

- Reproduction script: `scripts/analyze_l10_snp_panel.py`.
- Per-window statistics and quartet-frequency covariates: `results/l10_snp_panel_window_stats.tsv`.
- Joined local-tree tracks: `results/l10_tree_diversity_join.tsv`.
- Inside/outside effects and spatial intervals: `results/l10_inside_outside_spatial.tsv`.
- Correlations: `results/l10_tree_track_correlations.tsv`.
- LGI tracks with 1.5 and 8.0 Mb boundaries: `figures/l10_snp_panel_tracks.pdf` and `.png`.

The direction of high between-group and low within-group diversity inside the reported region is qualitatively compatible with restricted recombination maintaining divergent mating-type haplotypes (an MSRC prediction). This is only a conditional SNP-panel signal. The all-sites VCF/sample identity/callability validation, per-callable-base FST/dXY/pi, full 10-vs-12 sample comparison, and robust population-genetic inference remain outstanding. Do not use these conditional magnitudes as publication estimates of genomic diversity.
