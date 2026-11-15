# Stage 0 data audit — Channel Island monkeyflowers

**Status: blocked before local-data audit.** This document records source discovery only. No local trees were available to parse, and no downstream statistic or figure has been produced.

## Sources checked

- Short & Streisfeld (2025), “Disentangling Complex Histories of Hybridisation: The Genomic Consequences of Ancient and Recent Introgression in Channel Island Monkeyflowers,” *Molecular Ecology* 34:e17778, DOI [10.1111/mec.17778](https://doi.org/10.1111/mec.17778).
- Dryad dataset DOI [10.5061/dryad.cfxpnvxjw](https://doi.org/10.5061/dryad.cfxpnvxjw), published 2025-06-25. The Dryad API reports 49,728,125,513 bytes total storage.
- Source manifest: `source_manifest.tsv`.

## Data inventory from the published Dryad index

The Dryad README/index names the 100-kb window tree, metadata and TWISST files. The metadata TSV (92,841 bytes; file ID 4135122) and weight CSV (73,771 bytes; file ID 4135121) are top-level deposited files. The tree filename is described as a windowed output in the README, but its exact archive membership has not yet been verified. Dryad lists `recent_and_ancient.zip` (37,807,629 bytes), among other archives; it is not assumed to be the empirical tree archive without inspecting its contents.

Attempts to retrieve files through `https://datadryad.org/api/v2/files/{id}/download` returned HTTP 401 with “must have current bearer token.” The public `https://datadryad.org/downloads/file_stream/{id}` endpoint returned HTTP 403 in this execution environment. No files were saved locally. This is an external delivery/access limitation, not evidence that the dataset is private or unavailable to the user.

## Biological context established from the paper

The authors report evidence consistent with both recent hybridization between Island longiflorus and parviflorus and older gene flow involving ancestors of parviflorus and clades C/D, while noting alternative explanations including shared ancestral variation and ghost introgression. Simulations show that recent gene flow alone can resemble a recent-plus-ancient history in some statistics. The study reports 11 recent hybrid individuals excluded from subsequent analyses and 63 individuals retained for a species-tree analysis. These facts make a biologically labeled quartet design possible in principle, but do not justify interpreting every alternate local topology as introgression.

The paper identifies an approximately 300-kb chromosome 5 region around 24.7–24.9 Mb as a candidate for adaptive introgression from Island longiflorus into parviflorus. Its boundaries and coordinate convention must be confirmed from the paper's supplement and exact deposited tables before a region comparison.

The paper uses recombination estimates previously calculated in 500-kb windows by Stankowski et al. (2019), and relates these to 100-kb introgression statistics. This supports a recombination covariate in principle; it does not establish a recombination-suppressed interval or a local recombination rate specifically for the chromosome 5 candidate.

## Audit fields not yet determinable

Without the actual Newick and aligned metadata files, counts of windows/trees/chromosomes, coordinate system, taxon labels, branch-length units, missing taxa, topology resolution, duplicate/overlapping windows, and exact TWISST matches cannot be verified. Substitutions-per-site units cannot be assumed for the named files until the tree-generation documentation and branch lengths are inspected.

No quartet frequencies, spatial runs, autocorrelation, patristic distances, chromosome 5 contrasts, genome-wide sensitivity analysis, or Neurospora comparison were computed. No numerical output table or figure has been fabricated.

## Gate and recommendation

**NO-GO for inclusion based on this run; analysis is incomplete, not a biological rejection.** Resume only after obtaining the tree file and matching metadata/weights with verifiable provenance. Then audit the data first; proceed to topology and branch-length statistics only if IDs and coordinates align and the fixed taxon sets support the planned quartet interpretation.
