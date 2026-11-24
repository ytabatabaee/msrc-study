# Atlantic cod window treeness

Treeness follows Tabatabaee et al. (2026), Systematic Biology: **Treeness(T) = sum of internal branch lengths / sum of all branch lengths**. For each rooted MCC tree, internal means non-terminal edges whose child node is not the root; denominator sums every non-root edge (internal plus terminal). Zero-length edges contribute zero and are retained. Missing, nonfinite, negative, or zero-total lengths invalidate a tree; no values are imputed. A root stem is an edge above the root and is excluded; the input Newick has no such edge, and root.branch_length is not counted.

## Input and validation

Reused `data/atlantic_cod/processed/cod_window_trees.tsv`: 426 250-kb published MCC trees over LG01/LG02/LG07/LG12. Valid treeness values: 425/426. Taxa: 12 Atlantic cod population tips, consistent across all source trees; outgroups: zero. Rootedness comes from the BEAST/SNAPP time-tree workflow and crown/root-age prior; Newick strings lack an explicit [&R] marker. Branch lengths are time units (Ma), verified against the source Figure 4 scale. No posterior samples for each window are locally available, so posterior median and credible intervals are unavailable.

Zero-length edges across valid trees: 0. Invalid windows: LG07:LG07_023000001_023250000 (negative_branch_length:-0.003195445245). Ultrametric deviation is included per tree in the TSV. Outgroup sensitivity is not estimable from these trees because none include outgroups; thus outgroups cannot inflate the reported treeness.

## Inside/outside comparison

Inside means wholly inside the frozen suppressed-recombination interval; outside means no overlap; boundary windows remain in tracks but are excluded from contrasts. Contrasts report mean/median differences and a percentile interval from 5,000 moving-block resamples (8 adjacent windows; 2 Mb), plus exact circular shifts of the treeness series relative to the region mask with BH correction across four LGs. These methods retain local serial dependence; the shift test assumes a stationary circular placement null.

| LG | n inside | n outside | mean inside | mean outside | mean difference [95% block CI] | circular p | BH p |
|---|---:|---:|---:|---:|---:|---:|---:|
| LG01 | 66 | 43 | 0.877699157 | 0.428918338 | 0.44878082 [0.429992636, 0.4682053] | 0.018018018 | 0.0286738351 |
| LG02 | 21 | 70 | 0.853956778 | 0.440348259 | 0.413608519 [0.38644629, 0.440751182] | 0.0215053763 | 0.0286738351 |
| LG07 | 36 | 81 | 0.926831339 | 0.412113214 | 0.514718125 [0.497723187, 0.53134604] | 0.0166666667 | 0.0286738351 |
| LG12 | 53 | 51 | 0.777680895 | 0.471254144 | 0.306426751 [0.24518107, 0.365167007] | 0.113207547 | 0.113207547 |

Full details: `inside_outside_treeness.tsv`.

## Association with branch/time statistics

Per-LG Pearson and Spearman correlations compare treeness with the existing mean pairwise MRCA and ingroup root-MRCA times, plus the total branch-length sum of each same MCC tree; 95% uncertainty intervals use moving-block resampling (8 windows). The correlation with total length is mathematically coupled to treeness because that sum is its denominator. These associations are descriptive and do not establish causal mechanisms or phylogenetic bias. See `treeness_correlations.tsv`.

## Reproduction

Run `python3 empirical/atlantic_cod/analysis/treeness/scripts/run_treeness.py` from the repository root. Outputs: TSVs and figures under `empirical/atlantic_cod/analysis/treeness/`. Tests: `python3 -m unittest discover -s empirical/atlantic_cod/analysis/treeness/tests`.
