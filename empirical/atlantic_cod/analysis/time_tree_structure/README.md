# Atlantic cod SNAPP time-tree structure

This analysis uses only the existing 426 coordinate-indexed 250-kb SNAPP
maximum-clade-credibility (MCC) population time trees. It does not rerun SNAPP
or alter existing window analyses. Local trees are rooted time trees with a
consistent 12-population set; branch times are in Ma. Per-window posterior
samples were not identified, so every window contributes one MCC point tree.
The existing inventory lists supergene-wide posterior samples but not
window-matched samples.

## Method

For each valid ultrametric tree, pairwise MRCA age is measured backward from
the present as root age minus root-to-MRCA path length. The root stem is
excluded. Trees with missing or negative branch lengths, zero root age,
inconsistent taxa, or non-ultrametric root-to-tip distances are excluded, not
repaired. Zero branches are preserved. The validation table includes the one
published LG07 tree with a negative branch; it is not assigned a PCA vector or
triplet results. No outgroups occur in the 12-tip local trees, so outgroup-only
PCA features are not applicable.

Absolute PCA uses the 66 unscaled MRCA ages. Relative PCA divides every age by
that tree's root age. Feature-wise z-scored PCA is retained only as a scaling
sensitivity analysis. PCA effects by LG use cyclic shifts of the frozen
inversion mask over ordered available windows. Hierarchical clustering uses
root-age-normalized vectors with average linkage and Euclidean and Manhattan
(SciPy `cityblock`) distances. The selected k is the highest silhouette among
k=2 through 7. Stability uses 80% contiguous-block subsamples (8 windows per
block) and adjusted Rand agreement with full-data labels. Cluster enrichment
uses independent cyclic mask shifts within each LG, equally weighting the four
linkage groups.

All 220 triples of the 12 populations are analyzed. For each triplet, r1, r2,
and r3 respectively mean (taxon1+taxon2)|taxon3,
(taxon1+taxon3)|taxon2, and (taxon2+taxon3)|taxon1, with taxa ordered
lexicographically. The supported topology is the unique pair with the youngest
MRCA. Its temporal separation is the oldest of the three pairwise MRCA ages
minus the sister-pair MRCA age. The `S_r1`–`S_r3` columns give that contrast for
each candidate resolution; the `Q_r1`–`Q_r3` values are the corresponding
signed contrasts. Values are also root-age normalized. A triplet is unresolved
when its age spread or nearest-pair difference is no greater than
max(1e-10 Ma, 1e-6 × root age); no polytomy is resolved artificially.

The published genome-wide 12-population tree provides a topological reference
only; its induced rooted triplet is compared with each local MCC topology.
Window-tree topology frequencies are frequencies among MCC point trees, not
posterior support for the genome-wide tree. Per-LG circular-shift tests preserve
the observed ordered statistic track and rotate the interval mask. BH-adjusted
values are provided by metric family. The four LGs are the genomic units; the
426 windows are not independent replicates.

## Main findings

- 425 of 426 local trees are valid, consistently rooted, and ultrametric. The
  excluded LG07 tree has a negative branch. No zero-length non-root branches
  were found.
- Absolute PCA PC1 explains 71.1% of variation and correlates strongly with
  root age (r=0.991), so it primarily captures absolute divergence-time scale.
  The normalized PCA retains structure (PC1 35.3%, PC2 28.0%), so the patterns
  are not only overall tree age.
- Root ages are higher inside LG01, LG02, and LG07 regions (inside-minus-outside
  means +0.749, +0.634, +1.091 Ma; BH q=0.0145 each). LG12 is +0.292 Ma but its
  cyclic-shift p=0.0857.
- Normalized-tree clustering is metric-dependent (best silhouette k=5 for
  Euclidean, 0.581; k=3 for Manhattan, 0.601), with high contiguous-subsample
  ARI. Several clusters are enriched/depleted inside intervals under the
  within-LG spatial null; exact counts and effects are in the results tables.
- The modal rooted triplet topology differs inside versus outside for 120/220,
  125/220, 106/220, and 130/220 triplets on LG01, LG02, LG07, and LG12,
  respectively. Mean agreement with the published tree is lower inside on all
  four LGs. See the report for per-LG percentages, separation shifts, and the
  most informative triplets.

These are local-tree contrasts, not a new joint population-tree inference and
not proof of genome-wide dating bias. Posterior uncertainty at the window level
cannot be assessed from the currently available data.

## Outputs

- `results/tree_availability_validation.tsv`: all 426 source rows and QC.
- `results/window_mrca_age_vectors_absolute.tsv` and
  `results/window_mrca_age_vectors_relative.tsv`: one row per valid window.
- `results/pca_scores.tsv`, `pca_loadings.tsv`, `pca_top_population_pairs.tsv`,
  and `pca_inside_outside_spatial_tests.tsv`.
- `results/cluster_assignments.tsv`, `cluster_silhouette.tsv`,
  `cluster_stability.tsv`, `cluster_composition.tsv`, and
  `cluster_spatial_enrichment.tsv`.
- `results/rooted_triplets_per_window.tsv` (425 × 220 rows) and
  `rooted_triplet_inside_outside.tsv` (220 × 4 LG summaries), including fixed
  r1/r2/r3 Q and S columns, topology frequencies, unresolved rates, and
  published-tree agreement.
- `results/rooted_triplet_rankings.tsv` and
  `selected_triplets_for_figures.tsv` document non-visual selection criteria.
- `figures/` contains Figure A (PCA, genomic PC tracks, cluster tracks,
  composition and observed representative trees) and Figure B (triplet
  frequencies and signed contrasts), each as PDF and PNG.
- `results/interpretation_report.md` provides quantitative findings and
  interpretation limits.

## Reproduction and tests

From the repository root:

```sh
python3 empirical/atlantic_cod/analysis/time_tree_structure/scripts/run_time_tree_structure.py
python3 -m unittest discover -s empirical/atlantic_cod/analysis/time_tree_structure/tests -v
```

Eight tests cover identical and different rooted trees, all three triplet
resolutions, uniform rescaling, zero-branch polytomies, non-ultrametric trees,
missing/inconsistent taxa, invalid root age, and root-stem exclusion.
