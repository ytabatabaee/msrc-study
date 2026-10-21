# Neoaves Phlag benchmark reconciliation

## Audit result

Commit `3b71262e7be30eed2b503e6bd3d28e9bf1edbdba` replaced the Stage 11B
summary metadata, comparison tables, text, and figure with a Stage 11-style
detector audit. The replacement detector rows report 6,096 original
multi-taxon observations, while official Phlag is defined on the 5,619-entry
matched four-taxon manifest. The replacement matched manifest and Phlag run
metadata therefore mixed valid Phlag counts with an unmatched detector
analysis and absolute paths.

The retained raw Stage 11B output, matched input, coordinate manifest, and
window predictions were inspected before restoration. The Stage 11B versions
from commit `2e59ef6` were restored for benchmark metadata and summaries.
Derived comparison tables and the figure were regenerated from source
predictions by `scripts/reconcile_phlag_benchmark.py`.

The file-level classification is:

| files | classification | disposition |
|---|---|---|
| `manifest.json`, `matched_locus_manifest.tsv`, `phlag_run_manifest.tsv`, `report.md`, `methods_text.md`, `results_text.md` | outdated or inconsistent replacement | restored to Stage 11B and path-normalized |
| `interval_comparison.tsv`, `boundary_accuracy.tsv` | inconsistent mixed-domain summaries | regenerated from all four labelled prediction sources |
| `our_detector_predictions.tsv` | scientifically valid detector output, but multi-taxon | historical Stage 11 copy retained; latest copy moved to the clearly labelled multi-taxon directory |
| `feature_scaling.tsv`, `spatial_sensitivity.tsv`, `phlag_reproduction.tsv`, `published_phlag_predictions.tsv` | harmless formatting or unchanged audit metadata | retained without scientific reinterpretation |
| `runtime_memory.tsv` | mixed run metadata | restored for the validated Phlag/matched benchmark; latest detector metadata retained in the multi-taxon directory |
| comparison PDF/PNG | inconsistent figure | regenerated with disconnected intervals and separate input-domain tracks |

## Retained analyses

1. **Official Phlag, matched four-taxon input.** Upstream Phlag v0.2.0 at
   commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6`, 5,619 input trees,
   5,619 finite state calls and posterior probabilities, 1,200 state-1
   windows, and seven intervals using the 250-kb physical-gap rule.
2. **Our matched four-taxon topology detector.** The categorical detector in
   `matched_detector_predictions.tsv`, using the same 5,619 local trees.
3. **Our original multi-taxon topology detector.** The preserved
   `stage11_original_multitaxon_detector/our_detector_predictions.tsv` row
   reports 25.02--56.82 Mb from 6,096 locus-level observations.
4. **Our original multi-taxon topology-plus-B detector.** The same preserved
   output reports 44.13--56.82 Mb from 6,096 observations and adds
   `B_observed`.

The two multi-taxon analyses are not substitutes for matched-input analyses.
Their copied 5,619-entry manifest and Phlag metadata are retained only as
provenance because they do not provide a verified 6,096-observation
coordinate manifest. They are included in the figure and tables as
unmatched-input exploratory comparisons.

## Validated official intervals

| interval | start (Mb) | end (Mb) | anomalous windows |
|---:|---:|---:|---:|
| 1 | 25.02 | 32.30 | 503 |
| 2 | 32.67 | 32.68 | 1 |
| 3 | 33.46 | 34.48 | 40 |
| 4 | 44.13 | 46.42 | 115 |
| 5 | 46.77 | 47.54 | 47 |
| 6 | 47.80 | 48.58 | 51 |
| 7 | 49.05 | 56.82 | 443 |

The official union comparison is IoU 0.792416, precision 0.905185,
recall 0.864142, and F1 0.884188 against the canonical structural-region
union. These synteny-based regions are candidate annotations used only for
post-fit evaluation; they are not experimentally verified inversion
breakpoints.

## Provenance and remaining issue

Repository-relative paths are used in the retained run manifest. The raw
Phlag outputs and prior invalid Stage 11 evidence remain in the benchmark
directory. The only unresolved issue is the provenance of the 6,096-locus
multi-taxon detector coordinate mapping: the committed detector summary has
no independently verified 6,096-row coordinate manifest. It is preserved,
labelled, and excluded from claims about a matched Phlag input or from fitting
against structural annotations.
