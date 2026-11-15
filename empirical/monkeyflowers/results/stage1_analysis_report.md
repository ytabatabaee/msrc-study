# Stage 1 analysis report — Channel Island monkeyflowers

## Decision

**NO-GO for inclusion in the MSRC paper from this run.** This is an operational NO-GO because the actual local-tree input could not be retrieved; it is not evidence that the biological dataset lacks a spatial genealogy signal. Do not cite this audit as a negative result.

## Answers to the requested questions

1. **Are the downloaded trees suitable for MSRC quartet statistics?** No trees were downloaded. Suitability is undetermined until taxon sets, topology resolution, coordinate alignment, branch-length definitions, and sampling structure are audited.
2. **Is there a reproducible spatially coherent genealogy shift?** Undetermined; no local-tree analysis was possible.
3. **Are topology and branch-length changes spatially associated?** Undetermined; no branch lengths were inspected.
4. **Are alternative-topology regions short/fragmented or sustained?** Undetermined; no runs were measured.
5. **Does excluding introgression-associated regions change genome-wide results?** Undetermined; no genomic-window frequencies were computed.
6. **How does the spatial pattern compare with Neurospora?** No comparison can be made. The existing Neurospora benchmark describes an approximately 7-Mb alternative-topology signal in a recombination-suppressed region; this monkeyflower run contributes no matched physical-distance measurements.
7. **Is the dataset suitable for inclusion?** Not on current evidence. Reassess after the source access blocker is resolved and all data-audit gates pass.

## Work completed

- Verified the article and Dryad dataset identity, the source file names and the Dryad-declared file sizes via the public landing page/API metadata.
- Confirmed that the dataset is 49,728,125,513 bytes overall and avoided raw sequencing files.
- Attempted the public Dryad API file download (HTTP 401: current bearer token required) and public file-stream endpoint (HTTP 403).
- Recorded provenance and the unresolved access issue in `source_manifest.tsv` and `stage0_data_audit.md`.

## Work not completed

No window-level audit, topology classification, branch-length analysis, spatial fragmentation statistics, chromosome 5 comparison, genome-wide sensitivity analysis, figures, or focused tests were completed because these require the actual Newick and matching data files. No tables were created for unobserved results. Continue once the exact empirical tree file, data TSV, weight CSV, and README are available locally and checksummed.
