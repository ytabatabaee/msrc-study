# Primula Stage 0 data audit

## Decision

**NO-GO for spatial MSRC inclusion using this deposit.** All 3,623 Newick records are usable and have complete branch lengths, but the 3,923 entries in chromosome partition metadata are alignment blocks, not tree identifiers. No validated tree-to-coordinate crosswalk or explanation for the 300-record count difference is deposited. More decisively for this use case, the Dryad README says the nuclear genome FASTAs exclude the hemizygous S-locus; its sequence/tree is analyzed separately.

## Gene-tree archive

`PrimPrim_phylogeny_allGF.trees` contains 3,623 valid Newick records. Each has the same 106 taxa and 209 numeric non-root edge lengths. There are no malformed records or non-root polytomies; each is written with an unrooted trifurcating Newick root. Across 757,207 non-root lengths, min=0, median=0.00946494, mean=0.0190781, max=1.631064; 40 are exactly zero and none negative. The fitted IQ-TREE branch lengths are substitutions per site, not CU.

No record has a source filename, coordinate, chromosome, or GF identifier. The file is bare concatenated Newick. Record ordinal is only an archive locator.

## What the partition files represent

Each row has syntax like `DNA, phylo2_noreps_filtP2_pve_haplotypeT_001:10050001-10100000.min4.fasta = 1-2410`. The filename embeds a 50-kb chromosome interval; the right-hand range gives the inclusive columns assigned to that interval in the chromosome-concatenated FASTA. These are alignment partition definitions, not gene-tree IDs. The Dryad README says chromosome FASTAs/partitions were used to generate trees from every other GF.

Across all eleven files there are 3,923 rows, all with unique chromosome/coordinate keys and inclusive 50-kb physical intervals. The site ranges tile each chromosome FASTA when sorted by site start. Coordinate line order is lexical, not numerical/physical. Physical start spacing is commonly 100 kb, consistent with 50-kb GFs sampled every other interval. The paper reports 3,623 selected alignments/trees, but no source list indicates which 300 partition entries did not yield an archived tree, whether they were filtered/failed, or whether file versions differ. This difference cannot be distributed among chromosomes from the unlabeled trees. It is not a Newick parse error or duplicate partition-key issue.

The hypotheses, evidence, and limitations are detailed in `mapping_reconciliation_report.md`. `tree_coordinate_mapping.tsv` now includes all 3,623 Newick ordinals with missing coordinates marked `NA`; no positional guess was made.

## S-locus and reference

Stubbs et al. mapped reads to the chromosome-scale haploid *P. veris* reference that contains the S-locus. Potente et al. locate an approximately 260-kb hemizygous region in the chr1 approximate 27.43–27.70 Mb bracket. The Dryad README states the nuclear genome FASTAs exclude the hemizygous S-locus and lists the S-locus alignment/tree separately. Thus the reference contains the locus, but the genome-fragment tree alignments do not represent its sequence according to the deposited README.

Three chromosome-1 partition blocks fall in the broad published bracket. Their alignment block lengths are 53, 55, and 228 columns, but they are not evidence of S-locus bases and cannot be linked to Newick records. Coverage and mapping statuses are tabulated in `s_locus_coverage_audit.tsv`. The S-locus tree's discordance from the nuclear phylogeny is a separate whole-locus result and not evidence of a spatial shift among the genomic-fragment trees.

## Stage 0 answers

- All 3,623 deposited gene trees parse and have complete substitution/site branch lengths: **yes**.
- Are 3,923 partition rows 3,923 trees? **No.** They are coordinate/site-block entries in the chromosome FASTA partition metadata.
- Why do counts differ by 300? **Not exactly recoverable from this deposit.** Filtering/failed inference or version differences are possible; no retained-window manifest/log resolves it.
- Can each Newick record be mapped to coordinates? **No.**
- Does the mapping count difference imply 300 missing biological windows? **No such conclusion is supported.**
- Is the reference genome's S-locus present? **Yes.**
- Do deposited GF trees directly include S-locus sequence? **The README says the nuclear FASTAs exclude it.**
- Which GF trees genuinely overlap/include the S-locus? **No record can be assigned; direct S-locus sequence is excluded by the documented nuclear alignment construction.**
- Proceed with spatial Stage 1? **No.** The direct locus sequence is absent from the genome-fragment series and coordinates are unlinked.
