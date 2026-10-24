# ASTRAL validation gate

No full 30-taxon ASTRAL analysis was run for the validated 100-history
extension. The prespecified gate requires an independently supported
discordant dominant quartet after accounting for the many quartet, switching,
and origin comparisons.

The current summaries contain one nominal history-level negative-dominant row,
but zero rows pass the conservative multiple-comparison gate (`z = 4.5`) in
`validated_100/quartet_reversal_validation.tsv`. That nominal row is therefore
not evidence of a population-level quartet reversal and does not authorize an
ASTRAL run.

When the gate is met in an independent validation batch, the analysis must
include the matched scenario A, an ordinary MSC control, identical full-tree
sampling, and the corrected coalescent-unit branch-length conversion. Until
then, `astral_results.tsv` remains the prior not-executed record.
