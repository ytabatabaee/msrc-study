# Diagnostic report

## Classification

**Verified:** extreme terminal branches, the 14.7375-Mb discontinuity, and the effect on B(w) are present in and reproducible from the frozen source Newick trees.

**Supported but not conclusive:** the downstream increase is driven by the focal sample's terminal branch and its pairwise distances; the arrangement-informed numerator changes much less than B(w).

**Verified filtering-specific anomaly:** in representative windows, the extreme focal branch is present in the coverage-plus-allele-ratio-filtered ML archive but not in the PASS-SNP or coverage-only ML archives. This rules out a generic property of the sample label or chromosome-wide biology as the sole explanation.

**Unresolved mechanism:** the available archive does not include the per-window alignments or complete IQ-TREE optimization provenance needed to distinguish an allele-ratio filtering consequence from a sparse/ambiguous alignment or ML optimization instability. The leave-one-out track is therefore a sensitivity analysis, not a correction.

`source_filter_comparison.tsv` compares the focal branch in representative windows across the three published ML filtering products where the nested archives were readable.
