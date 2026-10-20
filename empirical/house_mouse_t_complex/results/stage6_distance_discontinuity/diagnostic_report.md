# Stage 6 distance-discontinuity diagnostic

## Main result

The largest native 5-kb mean-shift statistic for B(w), scanned without using inversion coordinates, occurs at approximately **14.7375 Mb**. The top candidate score is 56.32; the downstream/upstream B mean fold change is 35.38.

The arrangement-blind B(w) recomputation agrees with the frozen Stage 4 values to a maximum absolute difference of 0. The downstream maximum-terminal branch is consistently `Mmd_FRA4_18B.variant`; the B(w) sensitivity excluding that tip is recorded in the diagnostics. A_relative is retained as A_distance/B and the source Stage 4 value is used for the normalization audit.

The transition is accompanied by the local-tree branch diagnostics in `window_tree_diagnostics.tsv`; zero/negative/nonfinite branch counts and retained-tip counts are reported per window. Alignment length, SNP count, missing-data fraction, and IQ-TREE command-line settings were not recoverable from the committed processed files, so a preprocessing explanation cannot be verified from this repository state.

The upstream/downstream Pearson correlations are: A_distance versus B = 0.7094 upstream and 0.0483 downstream; A_relative versus B = 0.2664 upstream and -0.2962 downstream; A_distance versus total tree length = 0.5248 upstream and 0.0416 downstream. The B sensitivity excluding `Mmd_FRA4_18B.variant` changes the downstream mean from 0.322949 to 0.007854, while the upstream means are 0.009128 and 0.009160.

## Interpretation status

- **Verified:** the discontinuity is present in the frozen ML patristic-distance track; it is reproduced directly from the 4,046 source Newick trees; it is not caused by an arithmetic mismatch in B(w) or A_relative.
- **Supported but not conclusive:** the transition's relationship to local branch-length components and within/between-subspecies distance components is quantified in the diagnostic tables and figure.
- **Unresolved:** whether the shift reflects a biological divergence pattern, SNP ascertainment, alignment/filtering change, or another tree-estimation difference. The available metadata do not establish uniform upstream IQ-TREE settings or alignment quality across windows.

These are substitution-unit distances from local ML trees. They are not calibrated coalescence times. The original A_distance, A_relative, and B(w) tracks remain scientifically interpretable as their defined local-tree summaries, with this unit/provenance limitation. The existing manuscript figures were not changed.
