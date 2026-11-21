# Atlantic cod root-MRCA analysis

This analysis calculates T_MRCA(w) as the age of the MRCA of all Stage-2 Atlantic cod population tips in each published 250-kb SNAPP MCC tree. Any outgroup tips are excluded by selecting the ingroup MRCA; all 426 source trees contain exactly the 12 ingroup tips and no outgroup tips.

## Tree and time validation

Validated: 425/426. One materially negative-branch tree was excluded from age estimation: LG07_023000001_023250000.
The tree strings lack an explicit Newick [&R] marker, but their source is rooted BEAST/SNAPP MCC output generated with a crown/root-age prior (`lognormal(0,3.83,0.093)` in the local XML-generation script). All usable ingroup trees are ultrametric within 1e-7 time units. The branch-length scale is verified against the original Matschiner et al. (2022) Figure 4 caption: it describes these 250-kb sliding-window profiles, gives a 0.5 Myr scale bar, and states that inset branch bars represent 50,000 years. Thus the Newick branch lengths are in Ma (million years), not years; values are plotted directly without rescaling.
The published source caption is available at https://www.nature.com/articles/s41559-022-01661-x. The window trees and frozen supergene masks both use gadMor2 coordinates, following the cod_phylogenomics window-tree workflow; no gadMor_Celtic inversion boundaries are mixed into these tracks.

## Root-MRCA inside/outside comparison

Intervals use the frozen cod_phylogenomics window classifications; boundary-overlap windows are excluded. Confidence intervals are 95% percentile intervals from a spatial moving-block bootstrap (8 adjacent 250-kb windows = 2 Mb); they quantify genomic-window uncertainty, not SNAPP posterior uncertainty. Significance uses exact circular shifts of the root-age track relative to the fixed inside/outside mask over the full 250-kb coordinate grid, retaining gaps and autocorrelation; BH adjustment covers the four linkage groups.

| LG | n in/out | mean inside/out | median inside/out | Δ mean [95% CI] | Δ median [95% CI] | mean ratio [95% CI] | median ratio [95% CI] | circular p (BH) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| LG01 | 66/43 | 0.8171/0.06138 | 0.827/0.05924 | 0.7558 [0.711, 0.7979] | 0.7678 [0.7009, 0.7907] | 13.3 [12.2, 14.6] | 14 [12.7, 15.7] | 0.00885 (0.01389) |
| LG02 | 21/70 | 0.6952/0.06131 | 0.7179/0.06035 | 0.6339 [0.5311, 0.7389] | 0.6576 [0.6051, 0.8156] | 11.3 [9.48, 13.4] | 11.9 [10.7, 14.9] | 0.01042 (0.01389) |
| LG07 | 36/81 | 1.165/0.07054 | 1.14/0.06673 | 1.095 [1.063, 1.129] | 1.073 [1.023, 1.136] | 16.5 [15.6, 17.4] | 17.1 [15.9, 18.5] | 0.008 (0.01389) |
| LG12 | 53/51 | 0.3636/0.07252 | 0.3632/0.06731 | 0.2911 [0.2477, 0.3345] | 0.2959 [0.2669, 0.3444] | 5.01 [4.34, 5.76] | 5.4 [4.91, 6.43] | 0.05505 (0.05505) |

## Relation to pairwise and topology tracks

A(w) is the existing opposite-arrangement minus same-arrangement mean pairwise MRCA contrast. The repository has no Atlantic-cod file literally named B(w); here B(w) denotes the already-produced Stage-7 `mean_time_all`, i.e. the arrangement-blind mean across all 66 pairwise MRCA times in each 12-population tree. Quartet frequencies are the existing Stage-4 arrangement, baseline, and third-topology fractions.

| LG | Δ root age | Δ B(w) | Δ A(w) | Δ arrangement quartet frequency | ΔD | topology circular p |
|---|---:|---:|---:|---:|---:|---:|
| LG01 | 0.7558 | 0.3871 | 0.757 | 0.7411 | 1.233 | 0.009174 |
| LG02 | 0.6339 | 0.3137 | 0.654 | 0.5805 | 0.8923 | 0.01099 |
| LG07 | 1.095 | 0.5872 | 0.4355 | -0.0268 | -0.005981 | 0.5556 |
| LG12 | 0.2911 | 0.1232 | 0.3246 | 0.6411 | 1.043 | 0.03846 |

LG07 specifically has a strong root-age shift despite little change in its dominant quartet/topology contrast: its root-age track and B(w) are evaluated in separate panels from quartet frequencies. This is a time-depth change without a matching shift in which topology dominates; the finer A(w) signal can still change because it measures arrangement-conditioned pairwise MRCA contrasts rather than topology frequencies. LG12 has a positive inside/outside root-age difference and block-bootstrap interval, but its circular-shift test is not significant after correction (BH p=0.055).
A higher root age and higher B(w) indicate deeper overall genealogies and deeper average pairwise coalescence, respectively. A(w) isolates the difference between opposite- and same-arrangement pairs. Their agreement is not a mathematical identity: root age can change through the deepest split even when most pairwise/topology summaries move little. Root age and B(w) are highly correlated within each LG (Pearson r 0.993–1.000); on LG07 both vary strongly despite essentially stable quartet frequencies.

## Limitations

The 426 values are published MCC point estimates, not independent posterior trees or independent genomic replicates. One invalid tree was excluded; no branch lengths were repaired. The root-age prior is shared across windows, root-time units are undocumented in the local source, and the block confidence intervals and cyclic-shift tests quantify spatial rather than posterior uncertainty. Cyclic shifts rely on a circular/stationary placement null. The window-level confidence intervals are sensitive to block length; see `root_mrca_block_sensitivity.tsv`.

## Reproduction

From the repository root run `python3 empirical/atlantic_cod/scripts/13_root_mrca_time_analysis.py`. Outputs are written under `empirical/atlantic_cod/results/root_mrca/` and the two PDF/PNG figure pairs under `empirical/atlantic_cod/figures/`.
