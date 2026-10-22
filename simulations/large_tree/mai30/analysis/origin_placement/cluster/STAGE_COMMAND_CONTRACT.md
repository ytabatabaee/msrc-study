# Bridges2 stage-command contract

The repository does not contain the optimized `msrc-sim` checkout, so the
cluster wrapper does not invent a simulator API. Before submission, provide
these environment variables using adapters tested against commit
`bfba58d627f6f467af649cdb0b37f969349e0aab`:

| Variable | Required output contract |
|---|---|
| `MAI_BENCHMARK_FORWARD_COMMAND` | One forward history, atomically checkpointed, with CPU time, walltime, peak RSS, output bytes, history count, and simulator commit recorded. |
| `MAI_BENCHMARK_BACKWARD_COMMAND` | One genealogy from each completed benchmark history, with the same resource fields and the rho value recorded. |
| `MAI_FORWARD_COMMAND` | The assigned manifest row's forward histories, one checkpoint per history, safely resumable. |
| `MAI_PANEL_COMMAND` | The deterministic geometry-only quartet panel, with its source tree/config checksums. |
| `MAI_BACKWARD_COMMAND` | The assigned row's six-rate or single-rate backward genealogy collection, reusing completed forward checkpoints and never regenerating them. |
| `MAI_ANALYSIS_COMMAND` | Existing-output analysis only: persistence, pattern frequencies, conditional and marginal quartet probabilities, margins, and uncertainty intervals. |
| `MAI_VALIDATION_COMMAND` | Independent validation gate using history-level replication and the prespecified multiple-comparison criterion. |
| `MAI_ASTRAL_COMMAND` | Optional full-tree ASTRAL analysis only after the validation gate; records ASTRAL version and command. |

Every adapter receives `MAI_STAGE`, `MAI_TASK_ID`, `MAI_TASK_JSON`,
`MAI_OUTPUT_DIR`, and `MAI_BASE_SEED=20261022`. It must honor the task's seed
policy, use atomic writes, and return nonzero on incomplete or invalid output.
The wrapper then writes the authoritative `.done.json` marker. A failed or
timed-out task must leave checkpoints and logs usable by a later resubmission.
