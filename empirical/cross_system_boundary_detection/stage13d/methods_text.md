# Stage 13d methods

The primary house-mouse benchmark uses the unmodified Phlag v0.2.0 executable
at commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6`. Six prepared
representative-taxon sets were evaluated on four eligible focal edges (N002,
N003, N004, N005), with the same representatives in all 4,046 windows of a
replicate. State and posterior sequences were accepted only when both had the
expected length and finite probabilities in [0, 1].

The five inversion references are scored on the observed chr17 domain
5,000,000--37,054,999 bp. Inv5 is right-censored. Accuracy uses merged
interval-set union/intersection metrics. Two fixed aggregations are reported:
the union of state-1 calls across edges and the maximum posterior across edges
with a prespecified 0.5 threshold. These rules were set before scoring.

The macOS `/usr/bin/time -l` wrapper returned nonzero because its optional
`sysctl` memory probe was denied. Complete Phlag output was nevertheless
present and finite for all 24 runs; peak memory is therefore unavailable rather
than imputed.
