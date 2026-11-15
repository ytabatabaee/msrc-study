# Neurospora LGI chromosome-context figure

`figure0_chromosome_context` compares published L10 intervals, lineage-specific inversions, and the local topology pattern on a single LGI axis. The plotted axis uses the *N. tetrasperma* 2509a reference orientation and is rounded to 9.5 Mb. The Stage 1 window data reach 9,464,227 bp; that last observed coordinate is not asserted to be the biological chromosome's exact length.

## Intervals and sources

- **L10 recombination suppression, 1.38–7.58 Mb (6.20 Mb):** Corcoran et al. (2016), Table 1, gives these bounds for L10 and defines its coordinates on the *N. tetrasperma* 2509a reference. Sun et al. (2017) independently reports suppression in all analyzed lineages and an overall extent of approximately 8 Mb / 86%, while emphasizing that interval size differs among lineages. The 1.38–7.58 Mb interval shown here is specifically the L10 interval; it is not a universal boundary.
- **L10 introgression-associated region, 1.38–7.58 Mb:** Corcoran et al. (2016), Table 1 and introgression analyses, links L10 mat-A to *N. hispaniola* and reports the interval. Martin & Van Belleghem (2017) describe the associated topology shift across approximately the full ~7 Mb nonrecombining LGI region, without an independent pair of nucleotide breakpoints. These published boundaries overlap the L10 suppression interval in Corcoran's table.
- **Rearrangements:** Sun et al. (2017), Supplementary Table 4, provides exact breakpoint coordinates aligned to the 2509 mat-a reference. The figure shows the LGI inversions from L4, L5/L6, L7, L8, L9, and L10 as separate lineage rows. L10's 3,632,717-bp inversion (4.332431–7.965148 Mb) was found only in strain CJ01 mat-A and is not fixed across L10. The L5 interval (0.732329–7.462114 Mb) corresponds to the reported pattern in L6 reference genomes. These are lineage- and mating-type-specific structural differences, not one shared arrangement across the species. The L10 inversion is 3.63 Mb, smaller than the 6.20-Mb L10 suppression interval, though it extends slightly beyond the reported right boundary.
- **Topology pattern inferred here:** Stage 1 100-kb block means show q3-dominant runs at 1.2–5.7 Mb and 5.8–8.2 Mb, with one intervening 100-kb non-q3-dominant block. This is derived from the same local-tree data being analyzed and is not independent biological validation or an independently published interval.
- **Suppression extent from Sun et al.:** Sun reports approximately 8 Mb / 86% across the chromosome but does not state one common coordinate interval in the main text; boundaries vary by lineage. This length-only datum is retained in `region_annotation_summary.tsv` and not drawn as an invented shared interval.

## Interpretation and use

The figure shows a broad L10 suppression interval and introgression-associated interval, alongside substantially smaller, lineage-specific rearrangements; the topology runs broadly overlap those published regions. This supports the spatial comparison, not a causal claim that introgression caused suppression or that rearrangements are irrelevant. Sun et al. argues that large-scale suppression predates the detected lineage-specific rearrangements and occurs in collinear regions.

This is a standalone supplementary/context figure. Existing spatial-topology and branch-length figures were not changed because adding a multitrack chromosome map to those quantitative panels would reduce their clarity.

## Sources

- Corcoran et al. (2016), *Genome Research*: https://pmc.ncbi.nlm.nih.gov/articles/PMC4817772/
- Martin & Van Belleghem (2017), *Genetics*: https://pmc.ncbi.nlm.nih.gov/articles/PMC5419486/
- Sun et al. (2017), *Nature Communications*, article: https://pmc.ncbi.nlm.nih.gov/articles/PMC5658415/
- Sun et al. (2017), Supplementary Information (Supplementary Table 4): https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-017-01317-6/MediaObjects/41467_2017_1317_MOESM1_ESM.pdf
