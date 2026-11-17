# Sunflower tree availability audit

Audit date: 2026-11-17. No sunflower files were present in the working tree before this audit. The three local Study B inputs are isolated under the git-ignored `data/sunflowers/raw/` directory.

## Findings by data class

1. **Actual deposited local Newick trees — YES.** The public repository `owensgl/helianthus_hybrid_species_2021` has `results/trees/ind_trees.tar.gz` (19,479,575 bytes), which lists 3,242 `.treefile` members. It also has `results/trees/gowens22.mappable.dp4.missing80.10.2kbphysical.treenames.trees.txt` (2,642,093 bytes), containing 3,242 tab-separated coordinate key/Newick records. All records parsed in the local inventory have 22 tips, numeric branch lengths, a terminal semicolon, balanced parentheses, and chromosome/start/end coordinates in the key. Every one matches a `.treefile` member in the archive. The trees are not inferred from a statement in the paper: the tree files themselves are deposited and downloaded.
2. **Topology summaries — YES.** The public analysis repository also contains ASTRAL/concordance tree outputs and Twisst topology-weight files; these are summaries, not substitutes for the coordinate-indexed local Newicks.
3. **Alignments — NOT verified in the public analysis repository.** Its tracked files include scripts/results and tree outputs; no 10-kb window FASTA alignments were found in the repository tree. The paper supplement download returned a Google CAPTCHA HTML response instead of the listed 2-MB ZIP, so it was not treated as an alignment source.
4. **Raw reads/VCFs — not required for this audit.** The papers point to SRA BioProjects and filtered SNP resources. These were not downloaded. Raw reads are not a necessary next step because the local Newicks are available.
5. **Published figures/local phylogenies from Study A.** `owensgl/haploblocks` includes `Figure_data/ED_Fig_3f_pet05.01.nexus` (9,882 bytes) and `Figure_data/ED_Fig_3f_pet11.01.nexus` (9,614 bytes); both contain NEXUS trees with numeric branch lengths. They are block-focused, H. petiolaris sample phylogenies, not the 2023 genome-wide window trees. Study A also presents a Bayesian haploblock phylogeny in its figure, but the two NEXUS files above are the only directly verified Newick-like source tree files in that repository audit. The paper prose reports map-detected recombination suppression at ten haploblocks, while the Extended Data Table 1 transcription contains 11 rows marked `g` (its legend defines `g` as map suppression); this inconsistency needs source-level reconciliation.

## Coordinate and branch-length verification

`10kb_window_sizes.txt` (506,107 bytes) contains 11,797 names of the form `Ha412HOChrNN-start-end.fasta` and called-base counts. The physical-tree file has 3,242 coordinate records. Its first entries include `Ha412HOChr01-101772574-101779479` and corresponding numeric Newick branch lengths. The tree archive contains exactly 3,242 `.treefile` entries. The local script records each tree's bp interval, tip count, branch-length count, and branch-length presence in `tree_coordinate_inventory.tsv`.

The published 2023 paper reports 3,267 retained windows, while the verified public tree table/archive contain 3,242 trees (25 fewer). The coordinate/called-base manifest has 3,273 entries at ≥2,000 called bases, 31 of which have no corresponding tree record. These counts do not reconcile exactly; use the deposited per-tree list, not the paper count, for analysis and resolve provenance of the 31 absent trees before sensitivity claims.

The 2023 paper specifies maximum-likelihood IQ-TREE 2.0.6 gene trees. The Newick branch lengths are valid numeric estimates in substitution units (substitutions/site for IQ-TREE). The inspected Newicks have 41 branches for 22 tips, consistent with fully bifurcating rooted trees; this is structural format verification, not an assessment of branch-length estimation error.

## Exactness caveat

Todesco et al. report the 37 HA412-HOv2 haploblock regions and sizes in Extended Data Table 1, with region endpoints rounded to whole Mbp; several blocks are represented by multiple separated segments. Those table entries support screening-level overlap counts, but they are not base-accurate suppression boundaries. The later Dryad `haploblocks.petfal.txt` (362 bytes, described as H. petiolaris haploblock coordinates) is individually listed by Dryad; its file-download endpoint returned HTTP 401 during this audit. No replacement interval source was downloaded. Thus the compatibility table's assembly match is established, while exact endpoint mapping for all focal blocks is not.

## Verified source URLs and byte sizes

- Study A paper and data statement: [Nature DOI](https://doi.org/10.1038/s41586-020-2467-6); [source scripts/data repository](https://github.com/owensgl/haploblocks).
- Study A exact-size selected haploblock trees: [pet05.01 NEXUS](https://github.com/owensgl/haploblocks/blob/master/Figure_data/ED_Fig_3f_pet05.01.nexus), 9,882 bytes; [pet11.01 NEXUS](https://github.com/owensgl/haploblocks/blob/master/Figure_data/ED_Fig_3f_pet11.01.nexus), 9,614 bytes.
- Study B article/supplement listing: [PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9907532/); supplement is listed as a 2-MB ZIP, but the attempted download returned CAPTCHA HTML and was removed.
- Study B results repository: [Owens et al. GitHub](https://github.com/owensgl/helianthus_hybrid_species_2021).
- [Compressed individual IQ-TREE outputs](https://raw.githubusercontent.com/owensgl/helianthus_hybrid_species_2021/main/results/trees/ind_trees.tar.gz), 19,479,575 bytes.
- [Coordinate-indexed Newick table](https://raw.githubusercontent.com/owensgl/helianthus_hybrid_species_2021/main/results/trees/gowens22.mappable.dp4.missing80.10.2kbphysical.treenames.trees.txt), 2,642,093 bytes.
- [10-kb window/called-base manifest](https://raw.githubusercontent.com/owensgl/helianthus_hybrid_species_2021/main/results/trees/10kb_window_sizes.txt), 506,107 bytes.
- [H. p. fallax interval file](https://doi.org/10.5061/dryad.bcc2fqznn), listed as 362 bytes; direct file download returned HTTP 401.

The `source_manifest.tsv` records the source listing and local file paths. Raw data are excluded from Git.
