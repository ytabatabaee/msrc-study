# Atlantic cod Phlag benchmark

The benchmark ran the upstream Phlag implementation independently for LG01, LG02, LG07, and LG12, with one HMM reset per linkage group and one run per eligible focal internal branch. The source commit was `76597d0e4ebd0e507fdafebf7f5163ef634a36d6` and the package version was 0.2.0.

The four linkage groups completed successfully. The labeled basal three-taxon node N009 was excluded because Phlag's own `collapse_basal_bifurcation()` removes it before focal-edge validation. Ten internal nodes were identified in the benchmark tree; nine were Phlag-eligible.

Phlag analyzed SNAPP MCC population trees, which are population-tree summaries rather than independent gene trees. Phlag's MSC-based QQS calibration and anomaly states therefore require cautious interpretation and are not equivalent to gene-tree anomaly calls.

The primary evaluation labels a window positive when its genomic interval overlaps the frozen inversion interval. Primary comparator rules were fixed before scoring: topology D > 0, time A > 0, and their conjunction. Branch-specific best values are oracle/post hoc summaries; the union is positive when any eligible branch flags a window.

No Phlag focal branch flagged any window under the upstream default state call. Consequently the Phlag union and all branch-specific Phlag metrics have F1=0 and IoU=0 in every linkage group. The full metrics are in `detector_metrics.tsv`, `phlag_branch_metrics.tsv`, `metrics_iou.tsv`, and `metrics_f1.tsv`.

## Reproduction

The isolated environment is `/private/tmp/phlag-venv`; the executable was installed from the upstream repository. Run:

`/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria/bin/python empirical/atlantic_cod/scripts/11_phlag_benchmark.py`
