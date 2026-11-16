# Primula Stage 0 data audit

## Decision

**NO-GO for empirical inclusion at this time; CONDITIONAL GO for acquisition and feasibility work.** Dryad metadata and the article methods verify that the intended dataset exists and that its genomic fragments are 50-kb, nonoverlapping, every-other-fragment samples (3,623 trees). The actual tree and partition payloads could not be retrieved in this execution: Dryad file-stream endpoints returned HTTP 403. Consequently, tree usability, branch lengths, tree-to-coordinate mapping, S-locus overlap, local genealogies, and all downstream statistics remain unverified. No numerical inference is made from unavailable files.

## Verified deposit metadata

Dryad lists `PrimPrim_phylogeny_allGF.trees` (19.59 MB), `PrimPrim_phylogeny_astral.tre` (896 B), `PrimPrim_phylogeny_astarl_allspecies.tre` (4.62 KB; spelling as deposited), `PrimPrim_phylogeny_slocus_iqtree.tre` (2.56 KB), `README.md` (2.03 KB), and chromosome partition files `PrimPrim_phylogeny_chr001.partition.txt` through `...chr011.partition.txt` (27.71–43.42 KB each). The deposit also lists large FASTA alignments, which were not requested or downloaded. These are listing sizes, not local file sizes. Checksums are unavailable because no payload was retrieved.

The Dryad listing says the collection includes complete nuclear, plastid, and S-locus alignments; it does not provide direct evidence here that the specific 3,623 tree records retain usable branch lengths or that the partition files encode physical coordinates.

## Sampling design and coordinate caution

The paper reports 50-kb nonoverlapping genome fragments from 11 chromosomes of the *P. veris* reference assembly, with every other fragment selected to create 3,623 alignments. Thus the selected locus windows are separated by unsampled 50-kb fragments under the described design. No interpolation across these gaps is appropriate. Partition filenames being chromosome-specific are not, alone, a verified mapping from each Newick record to a chromosome and coordinate interval.

## Questions not yet answerable from retrieved evidence

- Number of Newick records successfully parsed; branch-length completeness and units; taxa, roots, polytomies, and malformed records.
- Whether partition metadata link each tree to a chromosome, start, and end.
- Whether any local tree overlaps the S-locus or its documented recombination-suppressed interval.
- Whether the S-locus-associated discordance is visible in local trees or affects genome-wide inference.
- All requested quartet frequencies, spatial autocorrelation, run lengths, branch-length summaries, region comparisons, and figures.

## Acquisition and analysis blockers

1. Dryad lists the target files at https://datadryad.org/dataset/doi:10.5061/dryad.b2rbnzsgz, but file-stream downloads returned HTTP 403 from the execution environment.
2. The exact partition-file schema and whether it encodes physical positions cannot be inspected until the files are acquired.
3. S-locus annotation must be reconciled to the *P. veris* assembly coordinates before any overlap determination. Published *P. vulgaris* locus sizes cannot be substituted as *P. veris* coordinates.

No Stage 1 analysis or figure is supported yet.
