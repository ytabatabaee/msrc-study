# Recommended J48 mixture follow-up

The current grid has one independent with-replacement realization per epsilon and is sufficient to document sensitivity, but not to estimate transition uncertainty. Use a minimal two-part replication:

1. Run 5 independent seeds at epsilon = 0.02, 0.05, 0.075, 0.10, and 0.15. These values bracket the whole-tree transition and resolve whether the early individual split changes are stable.
2. Add 3 independent seeds at epsilon = 0.025, 0.05, 0.075, 0.10, and 0.125 if the first set shows disagreement in topology or focal split presence. Report the probability of each topology and of each focal split, rather than a single fitted threshold.

Use the same 62,945-record target and frozen taxon set. Treat the independent mixtures as independent resamples; they are not a nested sequence of accumulated loci. A fixed-background progressive-addition experiment and random genomic-region controls are useful secondary designs after this replication, but are not needed to interpret the current coarse grid.

Distinguish three outcomes: first appearance of an individual focal split, recovery of the complete J48_T0 topology, and the across-replicate probability of a topology change. Do not call any of these a recombination-suppressed genome fraction without an independent genomic denominator.
