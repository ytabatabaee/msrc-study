# Stage 11C Neoaves Phlag representative-sampling robustness

This benchmark extends the validated Stage 11B result without changing Phlag, the matched detector, or the segmentation rule. It evaluated 10 deterministic representative sets, including the original frozen replicate zero. The original seven Stage 11B intervals were preserved as the replicate-zero reference.

## Inputs and validity

The selected taxa, eligible-locus counts, and repository-relative input paths are in `representative_sampling_manifest.tsv`. Phlag execution status, exit code, state/posterior lengths, finite-posterior validation, runtime, and raw-output paths are in `run_manifest.tsv`. All 10 planned Phlag runs are valid. Invalid runs are excluded from accuracy summaries rather than counted as negative predictions.

## Aggregate comparison

| method | mean IoU | median IoU | SD IoU | range IoU | mean F1 | median F1 | SD F1 | range F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Phlag | 0.792109 | 0.792416 | 0.003282 | 0.783188-0.794614 | 0.883993 | 0.884188 | 0.002052 | 0.878413-0.885554 |
| matched four-taxon detector | 0.638987 | 0.638987 | 0.000000 | 0.638987-0.638987 | 0.779734 | 0.779734 | 0.000000 | 0.779734-0.779734 |

## Stability and interpretation

All seven original Phlag intervals are called in all 10 valid replicates. Six retain identical boundaries; one replicate splits the long seventh interval, and two nearby interval starts move by one or a few windows. The broad regions are therefore stable, with small fragmentation and boundary sensitivity. The matched detector returns the same single interval for all 10 representative sets and remains less accurate under the same matched inputs. Pairwise agreement uses the common eligible-locus support for each replicate pair; region stability and anomaly-frequency denominators are reported explicitly.

The benchmark measures sensitivity to representative selection, not ten independent biological datasets. Structural candidate regions are synteny-derived annotations, not experimentally verified inversion breakpoints, and were excluded from representative selection and detector fitting. The benchmark is complete for the available valid representative sets; the next methodological priority is the separately specified unified multi-interval detector.
