# Stage 13 extension: Phlag cross-system benchmark

## Result

The original Phlag implementation was preserved and audited at commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6` (v0.2.0). The bundled example has a finite reference output in `raw_results/upstream_bundled_example_stage12b.txt`; this is a reduced-iteration software check rather than a fresh full-default benchmark.

Atlantic cod reuses the unchanged Stage 11 default runs. Those runs are numerically degenerate on SNAPP MCC population-tree inputs: posterior anomaly probabilities are nonfinite and state paths are not interpretable. They are therefore classified as `degenerate_numeric`, with no Phlag accuracy score. This is not evidence that Phlag found no anomaly.

The fire-ant, Anopheles, and house-mouse inputs are classified as `incompatible_input` because their local trees contain individual, phased haplotype, or pseudo-t-haplotype tips while the available reference trees use group/species/population taxa. The original implementation does not provide a scientifically neutral multi-sample mapping layer. No arbitrary representative selection or silent collapsing was performed.

Consequently, this benchmark does not support a cross-system accuracy claim for Phlag. Existing unified MSRC detector scores are carried into `accuracy_comparison.tsv` for context, while invalid/incompatible Phlag runs are explicitly marked `NA`. A future valid benchmark requires gene-tree and species-tree taxon sets designed for Phlag, or a method extension validated independently of the target intervals.
