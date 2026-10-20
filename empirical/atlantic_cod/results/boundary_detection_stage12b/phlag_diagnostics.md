# Phlag diagnostic

The Stage 11 default runs were not rerun or thresholded. This audit reads their saved upstream outputs.

Raw runs inspected: 36. Every run has an all-zero Viterbi/state path, no finite posterior anomaly probabilities, `State divergence: nan`, and a final focal-edge estimate at `1e-06` for the focal edge reported by Phlag.

The output is therefore a numerical/calibration degeneracy on these inputs rather than evidence that a calibrated Phlag null state was preferred. Phlag's `compute_output()` writes smoothed posterior probabilities, and the saved values are all `nan`; its default state path is consequently not a usable anomaly probability benchmark here.

The likely mechanism is input-model mismatch: Phlag's MSC emissions and branch-length update are calibrated for gene-tree quartet concordance, while the benchmark supplies SNAPP MCC population-tree summaries. The local QQS values are highly discrete and the focal branch update collapses to the minimum branch length. No threshold was lowered and no original Stage 11 output was altered.

The upstream implementation does not emit per-window emission likelihoods or transition matrices in its result file. Those internals would require instrumentation of the external package, so this diagnostic reports the available posterior/state/path, edge-length, and state-divergence evidence without modifying Phlag.

## Phlag upstream toy gene-tree check

The bundled Phlag Neoaves toy example was run separately from the Atlantic-cod
benchmark using the original gene-tree inputs and precomputed QQS file. The
reproduction command used one outer iteration, five inner increment steps, and
10 simulation replicates to keep the diagnostic bounded:

`/private/tmp/phlag-venv/bin/phlag -s /private/tmp/phlag-upstream/test/neoaves.nwk -g /private/tmp/phlag-upstream/test/emission.gtrees -e N159 -L 1 -l 5 --n-replicates 10 --read-qqs-path /private/tmp/phlag-upstream/test/qqs.tsv -o /private/tmp/phlag-stage12b-toy.txt`

The run completed successfully on the bundled gene-tree example. This
distinguishes a working installation from the cod-specific numerical behavior;
the cod result remains classified as a degenerate application to SNAPP MCC
population-tree summaries rather than a reliable negative anomaly result.

## LG01
9 branch runs; state1 total=0; finite posterior values=0.

## LG02
9 branch runs; state1 total=0; finite posterior values=0.

## LG07
9 branch runs; state1 total=0; finite posterior values=0.

## LG12
9 branch runs; state1 total=0; finite posterior values=0.
