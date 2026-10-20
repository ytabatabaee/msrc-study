# Branch-length validation

The deposited chromosome-4 Newick records contain finite branch lengths in the original ML trees. This stage recalculates branch-length and patristic-distance summaries from the named trees. Values are substitution-unit distances from the deposited ML analysis, not calibrated divergence times. The fixed mapped-taxon track requires all mapped taxa to be present; incomplete trees are retained in the diagnostic table but are not treated as complete fixed-set observations. No branch rescaling is performed.
