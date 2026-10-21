# Neoaves Stage 11B repair report

The original Stage 11 Phlag comparison was invalid for genomic evaluation: its committed `matched_locus_manifest.tsv` contained only a header and `matched_inputs/Columbea.gtrees` was empty, while the retained raw output contained 6,096 state and posterior positions. The resulting zero Phlag interval was therefore a failed coordinate mapping, not a biological negative result. The invalid files are archived under `provenance_invalid_stage11/`.

Stage 11B regenerates the four-taxon local-tree input and manifest together from the original chromosome-4 Newick file. The validated input contains 5,619 trees and the same number of manifest rows, states, and finite posterior probabilities; 477 of 6,096 source Columbea loci lacked at least one selected representative. The manifest preserves the exact written-tree order, source identifiers, coordinates, selected taxa, and hashes.

Unmodified upstream Phlag v0.2.0 (commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6`) flagged 1,200 windows. At the prespecified 250-kb maximum physical-gap rule, these form seven predicted intervals: 25.02–32.30, 32.67–32.68, 33.46–34.48, 44.13–46.42, 46.77–47.54, 47.80–48.58, and 49.05–56.82 Mb. The unmerged run and all gap thresholds are retained in the Stage 11B TSV outputs.

The matched detector is explicitly a categorical four-taxon adaptation: one local tree contributes one resolved topology indicator, whereas the original 6,096-locus detector uses multi-taxon fractional QQS summaries. It selects 25.02–56.82 Mb under the fixed scan configuration. The original multi-taxon detector remains an unmatched exploratory result at 44.13–56.82 Mb. Structural-region metrics are descriptive comparisons to overlapping synteny-derived annotations, not independent proof of inversion breakpoints, and no method-superiority claim is made.

Phlag runtime and peak child memory are measured in `runtime_memory.tsv`; the corrected run used 5,619 loci and returned successfully. The official Phlag output is parsed only when the state and posterior sequences are each exactly the manifest length.
