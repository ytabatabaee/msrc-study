# Primula Stage 1 feasibility reassessment

## Decision

**NO-GO: do not include Primula as a spatial MSRC example from these deposited genomic-fragment trees.** The gene trees themselves are parseable with complete branch lengths, but the tree archive has no locus identifiers and cannot be joined to the FASTA partition metadata. The separate issue is decisive: the Dryad README says the nuclear genome FASTAs exclude the hemizygous S-locus. The deposited local-tree series therefore does not directly sample the supergene sequence whose recombination suppression motivates this analysis.

## Answers

1. **Why are there 3,923 partition entries and 3,623 trees?** The partition files are coordinate-to-site-block definitions for chromosome-concatenated FASTAs. They are not per-tree mapping tables. The Newick archive has 3,623 records, the number reported in the paper. The deposit has no retained GF list, failed/filtered tree list, tree-run logs, generating script, or order manifest to explain the count difference exactly. Filtering, failed inference, or file-version differences remain hypotheses; the 300 cannot be apportioned by chromosome.
2. **Can the exact mapping be recovered?** No, not from the deposited small files. The Newick records have no embedded identifier; the partition entries map FASTA coordinate blocks to alignment columns, not Newick records. The README's statement that trees were made from every other GF is not a serialized order specification, and the 300 count difference prevents a one-to-one ordinal join. See `mapping_reconciliation_report.md` and `tree_coordinate_mapping.tsv`.
3. **Does the reference genome represent the S-locus?** Yes. Stubbs et al. say the reads were mapped to the chromosome-scale haploid *P. veris* reference containing the S-locus. Potente et al. identify the roughly 260-kb locus in chr1's approximate 27.43–27.70 Mb region.
4. **Do the genomic-fragment alignments/trees contain that locus?** Dryad's README says nuclear genome FASTAs exclude the hemizygous S-locus and treats the S-locus alignment separately. Three chr1 partition rows fall in the broad coordinate bracket, but this is not proof that S-locus bases occur in those blocks; they may represent non-locus/flanking variants. No local Newick record can be tied to them. The number of verified tree overlaps is zero; the actual coordinate overlap count among unlabeled records is unknown.
5. **Is there enough spatial resolution?** The intended design samples 50-kb GFs every other interval, roughly 100-kb start spacing. A 260-kb locus would span only about two to three sampling intervals even if it were included. It is excluded from the nuclear FASTAs per the README, so this local-tree set cannot support inside-versus-outside topology, A(w)/B(w), local quartet branch lengths, autocorrelation, run lengths, or a locus leave-out species-tree sensitivity analysis. Window counts within 0.5/1/2/5 Mb in tree units are unknowable without a crosswalk.
6. **Does the S-locus tree differ from the nuclear history?** Yes, Stubbs et al. report shallow concordance but hard discordance at deep backbone relationships. That compares a separate S-morph-only locus phylogeny against nuclear phylogenies; it is not evidence of a localized shift in the genome-fragment series.
7. **Is there evidence that the S-locus changed genome-wide species-tree inference?** No. ASTRAL used all 3,623 GF trees, but none can be identified as locus-specific and the S-locus is a separate, excluded sequence alignment. The published discordance does not demonstrate inference bias.
8. **Could Primula meet the four empirical inclusion requirements?** Recombination suppression is documented via hemizygosity. But direct local trees are excluded from the deposit's nuclear GF FASTAs, the tree-coordinate mapping is not recoverable from available metadata, and the locus is small relative to sampling spacing. No substantial spatial topology or branch-length transition has been shown. It fails requirements 2–4 for this data product.

## Scope and next step

No Stage 1 spatial statistics or figures were generated because there are no defensible local genealogies to analyze. A full mapping exercise would require the original retained-GF file list and per-window IQ-TREE outputs/logs or an author-supplied crosswalk. Reconstructing thousands of trees or downloading large chromosome FASTAs is not justified for this MSRC question, since the documented nuclear alignments exclude the S-locus itself.

Primula can remain in background discussion as an example of discordance between an independently inferred supergene phylogeny and a nuclear phylogeny. It should not be presented as evidence for a spatially coherent genealogy shift or phylogenomic bias without a different dataset that retains the S-locus together with coordinate-verified local genealogies.

## Sources

- Stubbs et al. (2023), [paper](https://doi.org/10.1111/nph.18525), [Dryad deposit](https://doi.org/10.5061/dryad.b2rbnzsgz).
- Potente et al. (2022), [*P. veris* genome and S-locus](https://doi.org/10.1093/molbev/msac035).
