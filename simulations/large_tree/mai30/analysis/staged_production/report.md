# Mai30 staged production report

## Status

The four-replicate scaled Wright--Fisher pilot remains complete and unchanged.
The branch-length discrepancy in the separate suppression analysis was
corrected: raw SimPhy generation lengths had been used where scaled pilot
lengths were required. The corrected audit uses `scaled_length/(2*N_scaled)`
and matches the existing pilot ASTRAL convention.

All 300 original no-outgroup Mai collections were audited from their retained
500-tree files. A deterministic selection retained ten species-tree
replicates per condition. This produced a 90-row staged MSRC manifest, but no
staged collection or ASTRAL job was executed here.

## Answers at the current stopping point

1. The branch discrepancy was an analysis-layer factor-1,000 raw-generation
   versus scaled-time error, not evidence of a new demographic effect.
2. `kappa=1` is an unsuppressed MSRC treatment with finite structured
   backgrounds, not ordinary MSC.
3. The 0.15 condition is slow because long Wright--Fisher/structured
   genealogical histories require many backward event iterations; the
   incomplete local run does not justify a precise cost extrapolation.
4. Feasibility of the full 30-replicate staged experiment is not established
   until one-history Bridges2 benchmarks are run for all three conditions.
5. No staged collections were completed in this task; the prior 12 pilot
   collections remain the only completed MSRC collections.
6. Quartet distortion must be evaluated from the staged collections before
   deciding whether it precedes ASTRAL topology changes.
7. The four-replicate pilot does not support a systematic suppression claim.
8. The 300-replicate campaign is not ready to launch. The staged 90-task
   workflow is ready for benchmark-gated submission.

The staged analysis is therefore a validated preparation and correction
stage, not a completed production result.
