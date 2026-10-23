# Full 363-taxon comparison with the published PNAS trees

This audit parses the committed 363/48-taxon ASTRAL outputs and the cached
versioned PNAS repository archive. It does not invoke ASTRAL or regenerate
gene trees. All RF values are symmetric differences of canonical unrooted
internal bipartitions.

## Published reference files

The cached archive is the Zenodo `chr4avian-v1.0.0` repository snapshot
(0dc24ebde4dd904115bd1e56735550e792d5cf8f7af437dbadf809108c313763) at
https://zenodo.org/api/records/10699424/files/smirarab/chr4avian-v1.0.0.zip/content, corresponding to repository commit
`17b408a622c9589fe9f6ca8a434c150974e39cd4` (https://github.com/smirarab/chr4avian/tree/17b408a622c9589fe9f6ca8a434c150974e39cd4). The paper is
https://pmc.ncbi.nlm.nih.gov/articles/PMC11009670/. The source repository identifies `score-columbea/` as
the location of the species trees shown in Figure S9.

| Reference | Archive path | Taxa | SHA-256 |
|---|---|---:|---|
| Published full-data tree/reference artifact | `score-columbea/score.tre`; Git blob `aedc2d4890a7bff26e779724a05beb35777d6a0b` | 363 | `064e34ad6132d8c043650bf5aa5ae21528bb4b60e825afc37d91fb17db34a6d8` |
| Published outlier-removed tree/reference artifact | `score-columbea/score-post-astral-63430_filter_mRNA_random_region.gene.trees.aLRT0.05-by-63430_filter_mRNA_random_region.gene.trees.aLRT0.05.chr4outlierremoved.gene.trees.tre`; Git blob `f88f18219e915ddbed71161b0a25c542c72a3370` | 363 | `9357df58f4e25a693c9b30655a09e93badf4cf402209ba372e23a52ee29d5477` |
| Stiller `63K.tre` reference | `empirical/neoaves_chr4/external/stiller2024/63K.tre` | 363 | `d3b4a9f202bb5337338f6814ea4b64757d5ffa6547b8ff9e50faaa9260ed076f` |

The repository README documents `score-columbea` as containing the species
trees shown in Figure S9. The two files are therefore treated as published
species-tree reference artifacts. Their names and annotations show that they
were used in post-ASTRAL scoring, but the cached repository does not establish
whether each file is the original inference output or a scored representation
of that output. No ASTRAL command line or version identifier is included, so
their exact inference settings are not independently recoverable.

## Result

Our full-data `T0` is identical to both cached published full-data references
at the unrooted topology level (RF 0). The cached published full-data and
outlier-removed trees are also identical (RF 0). Our reconstructed `T_PNAS`
differs from the cached published outlier-removed tree by RF 2: one internal
split is removed and one is added. Therefore the apparent full-taxon
disagreement is real for the current reconstructed input, but it is not yet a
contradiction of the PNAS inference because our `T_PNAS` is based on
`named_recollapsed.tre`, whereas the published analysis used the original
aLRT-collapsed collection.

The changed split is the N62 relationship. In our T0 and in the published
full/outlier-removed references, the 9-taxon Columbimorphae split is present
and the N62 focal state is q1: Columbiformes and Other Columbimorphae are
together relative to Otidimorphae. In our reconstructed T_PNAS, that split is
absent and the 18-taxon Columbiformes + Otidimorphae split is present; the N62
state is q2. This is an unrooted split statement. It does not assert a rooted
clade without specifying a root.

## Interpretation of the PNAS comparison

The PNAS paper states that the 48-taxon analysis with all loci recovered the
Columbea/J2014-side topology, while removing the outlier region produced a
topology similar to S2024 but with Otidimorphae sister to doves, breaking
Columbimorphae. It also states that with the increased Stiller taxon sampling,
the S2024 topology was recovered with or without the outlier region. The exact
full-taxon references cached in `score-columbea` support that latter statement:
both published 363-taxon trees have the same N62 q1 state.

Our 363-taxon T0 agrees exactly with the published full-data topology. Our
reconstructed T_PNAS does not reproduce the published outlier-removed tree.
Thus the current data show a preprocessing/input reproduction discrepancy,
not a demonstrated contradiction of the published taxon-sampling result.

## Remaining explanation

The existing source audits establish that 63,113 of 63,430 collapsed records
have strongly supported header-derived locus correspondences; 317 remain
problematic, including 200 whose taxon sets are absent from the named archive.
The current T_PNAS has the same expected 61,999-locus count and the frozen
1,431-locus PNAS exclusion, so count and mask membership alone do not establish
identity. Ordinary aLRT branch collapsing cannot change leaf sets. The
remaining source/preprocessing discrepancies can therefore plausibly affect
the outlier-removed topology, but this audit does not demonstrate that the 317
records alone cause RF 2.

Our runs use ASTRAL-IV v1.25.4.8 with `-u 2`, `--genelength 1000`, and
`--sulengthtype 1`, as recorded in the run provenance. The cached PNAS archive
does not document a comparable ASTRAL version or command, so a version-only
explanation is unsupported. A controlled follow-up with the exact published
collapsed trees and the current ASTRAL binary is the minimal useful next test;
an additional inference using the current reconstructed input is not needed.

## Files

- `full363_pnas_tree_comparison.tsv`: RF and taxon-set comparisons.
- `full363_focal_split_comparison.tsv`: direct unrooted focal split and state calls.
- `full363_focal_branch_lengths.tsv`: CU, SU, localPP, and quartet annotations for homologous focal states.
- `full363_pnas_source_manifest.json`: source paths and checksums.
- `figures/neoaves_pnas_replication_comparison.*`: verified-results figure.
