# House-mouse Stage 8 report

The approximately 40-substitution-unit terminal branches are present in the source ML Newick trees and are not an arithmetic artifact in the Stage 4–7 calculations. The focal branch's downstream distribution and threshold counts are recorded in `branch_length_distribution.tsv` and `extreme_windows.tsv`; `source_filter_comparison.tsv` tests representative windows across the three published SNP-filtering products.

Among 1,298 upstream windows, focal-branch median, 90th, 95th, 99th percentile, and maximum values are 0.0002, 0.0012013, 0.001804, 0.003211, and 0.006939. Among 2,748 downstream windows, the corresponding values are 5.470856, 10.9336186, 26.7557431, 34.3524561, and 52.733546. The transition is therefore abrupt and broad downstream rather than a few isolated values at a fixed numerical cap.

The maximum terminal branch among all other tips is 0.016626 and the maximum internal branch is 0.028975 in the frozen set. The extreme pattern is therefore concentrated in the focal terminal branch rather than reflecting a general inflation of every tree branch.

The observed change point is 14.7375 Mb. All-sample B(w) changes from 0.00912753 upstream to 0.322949 downstream, while leave-one-out B(w) changes from 0.00916046 to 0.00785359. A_distance is recalculated consistently for both sample sets, and the ALL_TIPS topology track is recalculated separately.

The three-filter source comparison shows that the extreme focal branch is specific to the coverage-plus-allele-ratio-filtered ML product in representative windows; PASS-SNP and coverage-only trees do not show the same branch. This is a verified filtering-product-specific anomaly, while its exact mechanism remains unresolved without alignments and complete IQ-TREE settings. The recommended manuscript treatment is to retain the original all-sample B(w) as a reproducible audit track with an explicit filtering/provenance caveat, show leave-one-out and robust tracks as sensitivity analyses, and avoid presenting the extreme downstream B(w) as a calibrated time or ordinary within-species divergence estimate.
