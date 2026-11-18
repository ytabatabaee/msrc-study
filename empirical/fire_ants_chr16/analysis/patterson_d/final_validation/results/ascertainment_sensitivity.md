# BUSCO-gene ascertainment and influence

The nucleotide observations are SNPs in four selected BUSCO gene bodies per
author window. The 852 intervals span 1.75 Mb of target gene bodies in the
supergene and about 1.02 Mb in chr16A+B; these spans are not callable-base
denominators. Variant density reflects both the biology and the upstream
ascertainment/filter pipeline.

The frozen window-level D table allows delete-one-window calculations. Removing
any one of the 52 supergene windows leaves pooled supergene D between +0.8124
and +0.8311 for Q1 and +0.8370 to +0.8512 for Q2. Those deletions are
descriptive leverage checks on linked windows, not independent replication.

Excluding windows with variant-site call rate below 0.80 leaves supergene D
+0.8266 (Q1; 3,412 informative sites) and +0.8439 (Q2; 3,428 sites). The
corresponding outside D remains near zero for chr16A and is based on limited
chr16B data. More stringent callability cutoffs leave few or no outside
informative sites, so those estimates are not evidence of broad robustness.
The table `window_callability_exclusion.tsv` reports retained windows and counts
at 0.80, 0.85, and 0.90 thresholds.
The original scanner did not
retain per-gene ABBA/BABA summaries, so gene-level leave-one-out and balanced
gene resampling cannot be reconstructed from saved tables alone. No gene is
asserted to dominate or not dominate. Matched-information comparisons and
per-gene callability effects are likewise not established by this package.

There are 37 overlapping gene interval pairs in the coordinate bridge, ten
crossing author-window boundaries. Within-window overlaps are collapsed to one
window hit before pooling. Cross-window overlaps mean some sites can contribute
to two distinct window summaries; such windows are correlated, not duplicate
independent observations. This convention matches the existing analysis unit.
