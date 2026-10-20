# Stage 13 unified spatial boundary detection

This analysis applies one robustly standardized, penalized piecewise-constant detector to frozen topology and continuous tracks from Atlantic cod, fire ants, Anopheles, and house mouse. The fixed penalty is `2 log(n)` per added free-mean segment after MAD scaling. Exact dynamic programming optimizes each fixed segment count; the selected number of intervals is chosen from zero through four. Observed rows retain their coordinates and are not interpolated across gaps.

The circular-shift tests assess alignment with the frozen reference intervals. The block-permutation change-point tests assess existence of a scanned change point and are reported separately; they are heuristic under strong spatial autocorrelation. Branch-length quantities are kept in their source units. Atlantic cod uses SNAPP population-tree time summaries; the other distance tracks use substitution-distance units.
