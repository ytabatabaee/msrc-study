# Real nucleotide Patterson's D analysis: Neoaves chromosome 4

## Inputs and validation

The original Stiller 63,430-intergenic-alignment archive was downloaded from the current ERDA path:
`https://sid.erda.dk/share_redirect/ENhZODU9YE/01_alignments_and_gene_trees/intergenic_regions/63430.alns.tar.gz`

SHA-256: `87a2a643083bb567022683d99b4f5446cf9d52b1da0b0d694bd52aca55f37e43`; compressed size 3,645,947,507 bytes. The Stage-4D manifest contains 63,430 loci and 6,247 chr4 loci. All 6,247 chr4 IDs matched alignment members exactly after normalizing only the archive's `.1k.start` and manifest `_1k_start` naming and removing the tree-only suffix. No positional matching was used.

The three frozen PNAS intervals are the published GalGal4 window-start inclusive intervals 25,030,000--32,670,000, 33,510,000--34,470,000, and 44,130,000--56,810,000. The local locus coordinates are retained from the Stiller manifest; interval membership follows the existing repository's verified published predicate and should not be read as a nucleotide-resolution breakpoint lift.

## Frozen designs

Primary designs use one representative per published clade and *Gallus gallus* as outgroup: Columbea = `Columba_livia, Phoenicopterus_ruber, Ciconia_maguari`; N61 = `Tauraco_erythrolophus, Columba_livia, Ciconia_maguari`; N62 = `Columba_livia, Pterocles_gutturalis, Tauraco_erythrolophus`. `Meleagris_gallopavo` was frozen as a galliform sensitivity outgroup. These choices were written before reading nucleotide D values.

## Results

D was calculated directly from nucleotide site patterns with outgroup state A and derived state B: ABBA = A/B/B/A and BABA = B/A/B/A. Gaps, IUPAC ambiguity, and missing calls were excluded. Sites with more than two states or no ABBA/BABA pattern were uninformative. Regional D pools ABBA and BABA counts before taking the ratio. Uncertainty is a genomic block jackknife; the main table uses 1-Mb blocks and sensitivity tables retain 500-kb and 2-Mb results.

* **Columbea**: inside D = 0.2164 (SE 0.0107; 1281 loci, 18152 sites), outside D = 0.2168 (SE 0.0052; 4446 loci, 67582 sites), inside-minus-outside = -0.0004 (SE 0.0122; descriptive z = -0.03).
* **N61**: inside D = -0.2582 (SE 0.0098; 1221 loci, 24683 sites), outside D = -0.0806 (SE 0.0040; 4427 loci, 73729 sites), inside-minus-outside = -0.1776 (SE 0.0106; descriptive z = -16.76).
* **N62**: inside D = 0.0394 (SE 0.0111; 1139 loci, 18871 sites), outside D = 0.0385 (SE 0.0041; 4303 loci, 74218 sites), inside-minus-outside = 0.0009 (SE 0.0120; descriptive z = 0.07).

The interval-specific results are in `manuscript_results_table.tsv`. Columbea is nearly unchanged after pooling all PNAS intervals; N61 shows a negative inside-minus-outside contrast; N62 is essentially unchanged. The three intervals are not equivalent in size or site count, especially PNAS_2, so the interval-specific estimates have different precision.

## Sensitivity and interpretation

The predeclared alternate galliform outgroup and transversion-only estimates are in `patterson_d_summary.tsv`. They are sensitivity analyses, not alternate post hoc sign choices. The signal is a nucleotide discordance statistic and is not by itself evidence for introgression, recombination suppression, or MSRC. It can also reflect incomplete lineage sorting, taxon sampling, ancestral polymorphism, alignment/quality effects, mutation-rate variation, or outgroup polarization error.

The results provide a real, reproducible Patterson-D baseline for comparison with the existing local-tree tracks. They do not independently establish that the PNAS structural intervals caused the observed D differences, and no claim of causal MSRC support is made from this analysis alone.

## Outputs

* `results/input_audit.json`, `results/locus_alignment_mapping.tsv`, `results/patterson_d_locus_results.tsv`
* `results/patterson_d_summary.tsv`, `results/patterson_d_1Mb_windows.tsv`, `results/manuscript_results_table.tsv`
* `figures/patterson_d_spatial_chr4.{svg,png,pdf}`
* `figures/patterson_d_inside_outside.{svg,png,pdf}`
