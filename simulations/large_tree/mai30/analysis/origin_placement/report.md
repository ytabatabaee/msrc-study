# Rearrangement-origin placement report

## What was completed

- Four deterministic origin scenarios were defined for all four Mai pilot
  trees.
- All 27,405 quartet subsets were enumerated for structural exposure.
- Forward Wright--Fisher persistence and terminal arrangement sorting were
  simulated with the authoritative `msrc-sim` implementation.
- The requested rho grid was recorded, including 0, 1e-4, 1e-3, 2e-3,
  1e-2, 2e-2, 1e-1, and 1.

## What was not completed

The 30-taxon backward genealogy discovery run did not complete locally even
after reducing to one history and one locus. The long-root and large-
population cases make direct genealogy simulation computationally expensive.
No quartet reversal, statistically supported negative margin, or ASTRAL bias
is therefore reported. The corresponding TSV files contain explicit
`not_executed` status rows.

## Scientific interpretation

Deeper and root origins clearly increase structural exposure. They do not, by
themselves, establish persistent polymorphism or quartet reversal: the
forward process can fix or lose the arrangement during the longer exposure.
The original internal origins should remain the Mai30 production setting until
a cluster-scale, checkpointed origin-placement experiment supplies validated
genealogy-level results. The doubled-root scenario is especially unsuitable
for local execution and should first be benchmarked on Bridges2.

This is exact only for the hierarchical simulator components actually run:
forward Wright--Fisher frequency histories and the specified origin rules.
It is not a full 30-taxon exact MSRC theorem.
