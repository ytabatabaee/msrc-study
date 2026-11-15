# Neurospora Stage 0/1 empirical analysis

This analysis uses the local Newick windows, coordinate table, and TWISST weights in Martin & Van Belleghem's deposited `neurospora.tar.gz`. The raw archive and extracted files live below `data/neurospora/raw/` and are excluded from Git.

## Reproduce

From the repository root, with Python 3, Biopython, NumPy, and Matplotlib installed:

```sh
python3 empirical/neurospora/scripts/run_stage1.py
python3 -m unittest empirical.neurospora.tests.test_stage1 -v
```

The script checks alignment of coordinate, weight, and Newick rows, normalizes the deposited TWISST subtree counts, records every window, and calculates four-tip branch metrics from the deposited trees. Topology weights in the data are counts out of the possible 120 L10 A-by-a pairings, not frequencies; the parser divides each complete row by its sum. `NA` topology rows are retained in the inventory and omitted only from topology summaries.

The primary independently published L10 suppressed-recombination interval is represented as 1.5–8.0 Mb on LGI (the approximate boundary shown for this lineage in Corcoran et al. and summarized in Martin & Van Belleghem). Flanks are 250 kb outside each end. Each source row is one 50-variable-site window; genomic windows are linked observations, so block summaries, block bootstrap, equal-count matching, and one-per-100-kb thinning accompany ordinary window-weighted summaries.

`CJ03` (mat A) and `CJ04` (mat a) are the focal matched L10 homokaryons for branch contrasts; `CJ15`/`CJ16` are a sensitivity pair. The quartet topology weights are the published TWISST estimates over the full L10 mating-type groups plus *N. crassa* and *N. hispaniola*. The topology labels are read from the deposited `.topos` file: q1 is the putative species topology; q3 groups mat A with *N. hispaniola*; q2 is the remaining resolution.

## Output map

- `results/stage0_data_audit.md`: archive and local-tree audit.
- `results/stage1_analysis_report.md`: quantitative findings and manuscript recommendation.
- `results/source_manifest.tsv`: file sizes and checksums.
- `results/gene_tree_inventory.tsv`: per-window coordinates and published topology weights.
- `results/tree_quality_inventory.tsv` and `tree_taxon_sets.tsv`: per-tree taxon-set ID, branch completeness, polytomies, and the deduplicated tip labels.
- `results/region_comparison.tsv`, `branch_length_inventory.tsv`, `branch_length_summary.tsv`, `lgi_100kb_blocks.tsv`, `linked_window_thinning.tsv`, and `sensitivity_summary.tsv`: derived analysis tables.
- `figures/`: four figures, each in PNG, PDF, and SVG formats.
