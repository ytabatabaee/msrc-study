# Theory comparison and scope

The exact finite-M quartet operator in
`simulations/quartet/moran_theory/` predicts `P(Q,Z)` for a shared
four-lineage hierarchical population and the full-species-tree recursion
extends this to the supported four-taxon tree.  It does not provide an exact
30-taxon distribution for the Mai pilot.  The present tables therefore use
induced four-taxon panels from the realized 30-taxon gene trees and are
simulation diagnostics, not a direct theorem evaluation.

The diagnostic decomposition uses the identity

```text
P(Q) = sum_z P(Z=z) P(Q | Z=z).
```

It separates natural terminal-pattern frequency from the topology composition
within each pattern.  The panel includes the four pattern classes requested
for a reference split, but the random quartet panel is not a replacement for
an all-quartet 30-taxon theory.

The exact finite-state forward Wright--Fisher persistence distribution was not
recomputed here.  The retained frequency histories provide direct empirical
estimates of polymorphism at branch boundaries and through recorded
generations.  A future calculation can compare those estimates with a
branch-specific finite-state transition matrix without changing the backward
simulator.

The expanded logarithmic switching-rate grid is a planned sensitivity design;
its rows are explicitly marked `planned_not_executed`.  No claim about a
panmictic limit or a parameter regime with topology reversal is made from
that grid.
