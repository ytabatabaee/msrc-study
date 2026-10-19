# Methods

We used the original Phlag 0.2.0 implementation at commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6` with its default command-line parameters. The rooted 12-population baseline topology was copied into a benchmark input and given deterministic internal labels based on descendant taxon sets; its topology, taxa, and branch lengths were preserved. Local SNAPP MCC trees were ordered by genomic midpoint and supplied one Newick tree per line, independently for LG01, LG02, LG07, and LG12.

Phlag was run once per eligible focal branch per linkage group. QQS files were computed by Phlag once per linkage group and then read by the subsequent branch runs. The upstream default state sequence (`1`) was the prediction; no inversion coordinate or truth label was passed to Phlag. A window was a primary positive truth if its base-pair interval overlapped the frozen inversion interval.

The union predictor flags a window if any eligible focal branch flags it. Branch-specific precision, recall, F1, and interval IoU are reported independently; the maximum branch result is explicitly an oracle/post hoc benchmark. Topology-only, time-only, and joint comparators use fixed sign rules based on the frozen Stage-4A and Stage-7 window tables.

Because these inputs are SNAPP population-tree summaries rather than independent gene trees, Phlag's inferred HMM states and MSC calibration are methodologically provisional. This benchmark evaluates execution and predictive behavior under the original implementation; it does not claim that the Phlag state has its original gene-tree interpretation for these data.
