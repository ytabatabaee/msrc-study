# Stage 13 unified spatial boundary detection

This analysis applies one robustly standardized, penalized piecewise-constant detector to frozen topology and continuous tracks from Atlantic cod, fire ants, Anopheles, and house mouse. The fixed penalty is `2 log(n)` per added free-mean segment after MAD scaling. Exact dynamic programming optimizes each fixed segment count; the selected number of intervals is chosen from zero through four. Observed rows retain their coordinates and are not interpolated across gaps.

Continuous information improved the common-detector mean IoU over topology for cod and Anopheles, while topology was better for fire ants and mouse. The joint result did not improve over the better modality. Arrangement-blind tracks worked well in cod and Anopheles but were weaker or differently shaped in fire ants and mouse. These results support complementarity as a dataset-specific possibility, not a general claim of joint-model superiority.
