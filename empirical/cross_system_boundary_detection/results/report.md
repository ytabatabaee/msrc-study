# Stage 13 unified spatial boundary detection

This analysis applies one robustly standardized, penalized piecewise-constant detector to frozen topology and continuous tracks from Atlantic cod, fire ants, Anopheles, and house mouse. The fixed penalty is `2 log(n)` per added free-mean segment after MAD scaling. Exact dynamic programming optimizes each fixed segment count; the selected number of intervals is chosen from zero through four. Observed rows retain their coordinates and are not interpolated across gaps.

The fixed detector gave mean union IoU/F1 of 0.627/0.678 for cod topology versus 0.784/0.846 for continuous features; 0.530/0.693 versus 0.000/0.000 for fire ants; 0.995/0.997 versus 0.997/0.999 for Anopheles; and 0.143/0.269 versus 0.117/0.200 for mouse. The joint result was essentially identical to continuous. These averages are descriptive because systems have different numbers of groups and mouse has four reference intervals.

The circular-shift tests assess alignment with frozen reference intervals. The separate change-point null is a quick cumulative-sum permutation diagnostic for existence of a change point and is heuristic under strong spatial autocorrelation. Atlantic cod uses SNAPP population-tree time summaries; the other tracks use source-specific substitution-distance units.
