# Primula mapping reconciliation

## Finding

The mismatch is real, but the available deposit does **not** establish that 300 trees were missing. It contains 3,923 chromosome-FASTA partition entries and 3,623 valid unlabeled Newick records. The partition files describe alignment blocks in the chromosome FASTAs: each row names a `phylo2_noreps_filtP2_pve_haplotypeT_###:<start>-<end>.min4.fasta` segment (a 50-kb genomic interval) and gives that segment's inclusive column range in the chromosome-concatenated alignment. They are coordinate-to-alignment metadata, not tree IDs or a direct tree-to-coordinate manifest.

The 3,923 rows have unique chromosome/start/end keys. Within each chromosome, the alignment-site ranges tile the alignment without gaps when numerically sorted by site start. Physical rows are commonly spaced 100 kb apart, consistent with selecting every other 50-kb genomic fragment; the partition text itself is lexically ordered and is not physical-coordinate order. The Newick file consists of bare Newick records without locus names, chromosome labels, coordinates, or record annotations. It has no embedded chromosome grouping that can be validated.

Stubbs et al. report selecting every other 50-kb GF and estimating 3,623 alignments/trees with IQ-TREE. The Dryad README says the GF trees came from every other GF and that chromosome FASTAs/partition files were used to make the trees. Neither source supplies the original file list, per-GF tree filenames, filtering/failed-inference list, IQ-TREE logs, concatenation order, or a tree-record-to-GF manifest. Thus, the archive-level 300 count difference cannot be allocated among chromosomes or regions, and its cause cannot be determined exactly from the deposit.

## Hypothesis checks

| Hypothesis | Evidence and contradictions | Result |
|---|---|---|
| A. 3,923 candidate windows, with 300 filtered or failed before tree estimation | SNP calling and filtering are documented; the partition segment names include `filtP2` and `min4`. However, neither the paper nor README identifies 300 exclusions, and no retained-window list or inference logs are deposited. | Plausible, unverified. Do not label 300 particular intervals as failures. |
| B. Partitions and trees are from different stages/versions | The partition names encode a filtered, repeat-masked chromosome FASTA; the README links those chromosome FASTAs to the GF-tree workflow, and the paper describes 3,623 GF alignments. No version history or alternate counts explain 3,923 versus 3,623. | Possible, but unsupported as the exact cause. |
| C. Newick record order implicitly follows chromosome/window order | The README says trees derive from every other GF, but does not state serialization order. Newick records have no identifiers. Partition lines are lexical rather than coordinate ordered, and the count mismatch prevents a positional join even under a putative order. | Not independently verifiable; rejected as a mapping method. |
| D. Partition rows are not one tree each | Confirmed: each line maps a physical 50-kb segment to a site interval in a chromosome FASTA; it does not name a Newick record. The README says the FASTAs were used to make trees, but does not provide the one-to-one retained-tree subset. | This corrects the earlier interpretation. The 3,923 entries are FASTA alignment partitions, not a tree-coordinate map. |
| E. Parsing, duplicate rows, or coordinate convention caused the difference | All 3,623 Newick records parse; no malformed or extra records were found. Partition coordinate keys are unique; their site intervals tile each concatenated alignment. Physical intervals have consistent inclusive 50-kb lengths. | Not explained by parser failure, duplicate partition rows, or a simple coordinate convention. |

The publisher's supplementary-information index lists Methods S1–S4 (sampling, network estimation, D statistics, and ancestral-range estimation), not a genomic-fragment generation workflow. A direct retrieval attempt for the listed 6.7-MB PDF returned HTTP 403. Searches for the distinctive Dryad filenames and a linked public tree-generation repository did not locate scripts or logs. The extra all-species ASTRAL tree is a species-tree result, not a coordinate manifest.

The cross-chromosome distribution of the 300 unmatched alignments/records is **unknown**: there is no chromosome label in any Newick record. Counts of partition entries by chromosome are available in `partition_coordinate_summary.tsv`, but those are not per-chromosome tree counts.

## S-locus representation

The reference used for read mapping is the chromosome-scale *P. veris* maternal/haploid assembly that contains the S-locus. Potente et al. describe a roughly 260-kb hemizygous region in chr1's approximate 27.43–27.70 Mb pericentromeric bracket. The Stubbs Dryad README separately states that the nuclear genome FASTAs exclude the hemizygous S-locus, while the S-locus alignment/tree is a separate product made from S-morph accessions plus an outgroup. Therefore reference-coordinate presence must not be confused with sequence representation in the 3,623 genome-fragment gene trees.

Three partition entries fall within the approximate chr1 bracket (27,450,001–27,500,000; 27,550,001–27,600,000; 27,650,001–27,700,000). Their site-block lengths are 53, 55, and 228 columns, respectively, in the chromosome FASTA partition metadata. Those counts do not establish that S-locus bases are present; the README says the S-locus is excluded. Nor can any of the three blocks be joined to a Newick record. The separate 255,577-site S-locus alignment supports its published phylogeny, not a local sliding-window series.

The chromosome FASTAs are listed by Dryad but are large (chr001 is about 131 MB). They were not downloaded: they cannot restore the missing association between bare tree records and blocks, and the README already documents S-locus exclusion. Recovering the mapping from sequences alone would require re-estimating/comparing trees or access to original per-GF files/logs, beyond a justified pilot and not warranted for the S-locus MSRC question.

## Determination and next evidence needed

The exact source of the 300-count discrepancy remains unresolved, and the exact tree-to-coordinate map cannot be recovered reliably from the deposited small metadata/tree files. `tree_coordinate_mapping.tsv` contains one row for each tree record with chromosome/start/end explicitly `NA`; ordinal is only a locator within the Newick archive, not a source gene-tree ID. No candidate interval is assigned to a tree.

A future mapping attempt would require the original retained-GF list/per-window IQ-TREE outputs, a source script that reproduces and records the order and exclusions, or an author-supplied crosswalk. For this MSRC application, however, the nuclear GF FASTAs' documented S-locus exclusion means the local gene-tree series does not provide direct genealogies of the hemizygous supergene. Even a recovered coordinate map would mainly describe nearby background windows, not the supergene genealogy itself.

**Conclusion: NO-GO for Primula as an empirical spatial MSRC example using this archive.** Keep the separate S-locus versus nuclear-tree discordance as biological context only; it does not demonstrate a spatial local-tree shift or genome-wide inference bias.

## Sources

- Stubbs et al. (2023), [paper and methods](https://doi.org/10.1111/nph.18525), [Dryad deposit and file listing](https://doi.org/10.5061/dryad.b2rbnzsgz).
- Potente et al. (2022), [*P. veris* chromosome-scale assembly and S-locus analysis](https://doi.org/10.1093/molbev/msac035).
