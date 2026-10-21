# Stage 13d Phlag benchmark report

## Completed

- House mouse: 24/24 matched-taxon runs completed with finite, length-matched
  posteriors and nonconstant state paths: six representative sets and four
  focal edges.
- Valid aggregate performance across representatives is summarized in
  `house_mouse_replicate_summary.tsv`; union-state-1 mean IoU/F1 are
  0.24749412238111126 / 0.38718836762449044,
  and max-posterior mean IoU/F1 are
  0.25398101078063656 / 0.39610415834298857.
- Inv5 was evaluated only through 37,054,999 bp and is explicitly right-censored.
- The repaired Neoaves Stage 11B pilot remains preserved as a separate valid
  5,619-locus result; this Stage 13d run does not relabel it as a multi-replicate
  result.

## Numerical failures

The Atlantic cod, Anopheles, and fire-ant Phlag runs remain unavailable for
accuracy. Their raw outputs contain a focal edge at the Phlag floor (`1e-6`)
and nonfinite state-divergence/posterior values. These are numerical or input
calibration failures, not zero detections. The diagnostics identify the first
observable failure and preserve the raw outputs.

## Fair interpretation

The house-mouse comparison is valid for the matched representative inputs and
uses the original Phlag algorithm. It does not support a cross-dataset claim
that MSRC outperforms Phlag, because the other three datasets do not produce
valid Phlag posteriors and the Neoaves result is a separate one-pilot
experiment. The house-mouse results belong in a supplementary benchmark until
additional independent biological calibration is available.

## Reproduction commands

```bash
/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria/bin/python empirical/cross_system_boundary_detection/stage13d/scripts/complete_phlag_benchmark.py --replicates 6 --workers 4 --run
/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria/bin/python empirical/cross_system_boundary_detection/stage13d/scripts/complete_phlag_benchmark.py --score-only
/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria/bin/python empirical/cross_system_boundary_detection/stage13d/scripts/finalize_phlag_benchmark.py
```
