# Bridges2 stage contract

The driver dispatches to the versioned Python adapters in this directory. No manual `MAI_*_COMMAND` variables are supported.

| Stage | Entry point | Required validated output |
|---|---|---|
| forward / benchmark_forward | `run_forward_histories.py` | One gzip JSON checkpoint per independent history, with simulator commit and configuration hashes. |
| persistence | `run_forward_persistence.py` | Separate boundary persistence, joint ancestral persistence, and sampled terminal-pattern tables from completed checkpoints. |
| panel | `select_quartet_panel.py` | Deterministic quartet panel with canonical induced unrooted species splits. |
| backward / benchmark_backward | `run_backward_genealogies.py` | One parseable Newick genealogy per completed forward history, with history IDs and rho. |
| analysis | `analyze_quartets.py` | Normalized three-topology probabilities, both pairwise margins, and history-level observations. |
| validation | `validate_quartet_reversal.py` | Multiple-testing-adjusted independent-history gate requiring a discordant topology to exceed both alternatives. |
| astral | `run_optional_astral.py` | Validated all-tree input and nonempty Newick species tree from the site ASTRAL executable. |

Each adapter returns nonzero on missing prerequisites, invalid files, count mismatches, duplicate IDs, taxon mismatches, or hash mismatches. The driver then performs a second output-manifest and checksum validation before writing the completion marker.
