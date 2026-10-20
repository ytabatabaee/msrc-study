# Multi-resolution methods

The native 5-kb rows are the frozen ALL_TIPS quartet and Stage 4 distance observations. For each 25-, 50-, 100-, and 250-kb bin, all source windows whose starts fall inside the common bin boundaries anchored at 5 Mb were included. Quartet support was pooled by reconstructing category counts from each window's support fraction times its total quartet count (`n_resolved_quartets + n_unresolved_quartets`), then dividing pooled counts by the pooled denominator. This preserves the source unresolved-quartet denominator. Distance features were arithmetic means of valid local-tree observations.

Each bin records observed window count, expected 5-kb window count, observed base-pair union, and coverage fractions. Missing sequence is retained as missing; no values are interpolated. Inversion boundaries are plotted at their frozen coordinates and are not snapped to bin edges. This is a figure-readability sensitivity analysis only; boundary detection was not rerun.
