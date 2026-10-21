# Mai30 production readiness

## Completed pilot

The four-replicate scaled-Wright–Fisher pilot is complete. There are 12
collections, three suppression treatments for each of replicates 007, 042,
074, and 005, with 25 independent histories and 20 loci per history. Every
collection has 500 validated 30-taxon gene trees, 14,500 coalescent events,
and a history-to-locus mapping. The six collections for 007 and 042 were
reused; missing histories for 005 and 074 were filled by the existing
resumable runner.

ASTRAL-IV v1.25.4.8 was rerun for all 12 collections with
`-t 2 -u 1 --length CULength`. The result is in
`analysis/results/mai30_scaled_wf_astral.tsv`. Gene-tree and distribution
summaries are in `mai30_msrc_vs_species.tsv` and
`mai30_msrc_vs_original_distribution.tsv`. These rows are the pilot result;
they do not support population-level inference from four species-tree
replicates.

## Resource observations

The completed collections had summed history runtimes of approximately 96--
2,130 seconds per 500-tree collection and recorded peak RSS of 71--174 MB.
The spread is driven by the recovered Mai demographic histories and finite
Wright–Fisher population sizes, not only by the suppression fraction. The
original-scale generation-by-generation process remains infeasible locally.

A one-history, 20-locus probe for the third no-outgroup condition,
`outgroup.0.species.0.15.genes.0.15`, replicate 001, no suppression, was
stopped after 562.8 seconds without producing a validated history. It remains
an incomplete local benchmark. This probe did not modify the pilot outputs.
Measured pilot timings and this lower-bound probe are recorded in
`mai30_production_resource_estimate.tsv`.

## Planned production

`datasets/processed/mai30_production_manifest.tsv` contains exactly 900
collection tasks: three no-outgroup conditions, 100 replicates per condition,
and three treatments. It uses the validated mapping
`N_e = N_H/2`, `c=0.001`, initial frequency 0.5, 25 histories, and 20 loci
per history. Paths are repository-relative and task seeds are deterministic.

`cluster/run_mai30_production.slurm` runs one collection per array task with a
12-task concurrency limit. The collection runner writes per-history logs,
reuses only validation-marked histories, validates Newick/taxon/coalescence
contracts before aggregation, and writes `collection_complete.json` only after
the combined 500-tree file is checked. `cluster/run_mai30_production_astral.slurm`
is a separate validation-gated ASTRAL array.

Neither production array has been submitted. Before submission, benchmark a
representative condition-0.15 task on Bridges2 and confirm walltime, memory,
and storage allocations. The pilot demonstrates that local completion of all
900 tasks should not be assumed.
