# Arrangement-blind track validation

- Unified window rows: 4,046; source Stage 2 ALL_TIPS rows: 4,046; source Stage 4 rows: 4,046.
- Coordinate sets and ordering match: True.
- `A_distance`, `A_relative`, and `B=mean_pair_distance_no_spretus` match the frozen Stage 4 values to the 1e-15 check threshold.
- The arrangement-blind track uses all non-`M. spretus` pairwise distances and does not use t-versus-standard labels.
- Native 5-kb and separate 250-kb figures were generated. The prior Stage 5 figures and `binned_250kb_tracks.tsv` were preserved.
- 250-kb aggregation uses pooled ALL_TIPS quartet counts and arithmetic means of valid distance observations; bin coverage is recorded in the new arrangement-blind table.
