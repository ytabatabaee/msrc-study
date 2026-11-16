# Primula Stage 0 data audit

## Decision

**NO-GO for spatial MSRC inclusion from this deposit as currently cross-referenced; CONDITIONAL GO for a topology-only, genome-wide gene-tree comparison.** All 3,623 deposited gene trees are usable with complete branch lengths, and all 11 chromosome partition files are present. However, the partition metadata enumerates 3,923 genomic windows—300 more than the Newick records—and the Newick records have no GF identifiers. Neither the README nor the deposited metadata explains which 300 partition windows lack trees. Therefore tree index cannot be assigned to physical coordinates without an unsupported ordering assumption. No S-locus-overlapping tree can be identified from the current files.

## Local gene-tree inventory

`PrimPrim_phylogeny_allGF.trees` contains 3,623 Newick records, all usable, with one invariant 106-tip taxon set. There are no malformed records, no missing taxa, and no non-root polytomies. Every record is represented as an unrooted tree with a trifurcating Newick root (the root degree is not counted as a biological polytomy). Each tree contains 209 non-root edges, and all 209 have numeric branch lengths.

Across the 757,207 non-root edge lengths, the minimum is 0, median 0.00946494, mean 0.0190781, and maximum 1.631064; 40 edges are exactly zero and none are negative. These are IQ-TREE maximum-likelihood branch lengths in substitutions per nucleotide site, not species-tree coalescent units. The paper reports IQ-TREE v2.1.2 with ModelFinder and ascertainment-bias correction for SNP-only alignments. These branch lengths should not be called divergence times or CU values.

The exact taxon labels are in `tree_taxon_sets.tsv`; the per-tree parse and branch-completeness inventory is `gene_tree_inventory.tsv`; aggregate tree-quality numbers are in `tree_quality_summary.tsv`.

## Partition and coordinate audit

All 11 chromosome partition files are present. Their 3,923 rows each contain an explicit 50-kb physical interval and its alignment-site range. Every chromosome's alignment-site ranges are contiguous in partition-file order, and every row is 50 kb. The sum is 3,923 partition windows, versus 3,623 Newick records. This is a difference of 300, not an explained one-to-one map. The metadata does not label which windows have corresponding Newick records, and the Newick file itself has no record IDs. We do not infer that trees 1–3,623 correspond to the first 3,623 partition rows.

The partition-file line order is not physical-coordinate order. Coordinates must be sorted by numeric start after a verified join. Among the partition intervals, 100-kb start spacing is common but some consecutive intervals are farther apart; these gaps must not be interpolated. The per-chromosome audit and complete partition-window inventory are in `partition_coordinate_summary.tsv` and `partition_window_inventory.tsv`. Rows in that inventory deliberately have `tree_index=NA` because tree association is unresolved. `mapping_reconciliation.tsv` records the 300-row mismatch, and `tree_coordinate_mapping.tsv` remains header-only.

## S-locus location and candidate windows

Potente et al. (2022) describe the *P. veris* S-locus as an approximately 260-kb region on chromosome 1, present only in the maternal (S-morph) haplotype and therefore hemizygous; the paper places it within a putative pericentromeric chromosome-1 bracket around 27.43–27.70 Mb. That bracket is regional context, not a verified exact S-locus start/end. The Stubbs et al. partition metadata contains three chr001 50-kb intervals inside that bracket: 27,450,001–27,500,000; 27,550,001–27,600,000; and 27,650,001–27,700,000. These are **candidate partition windows**, not verified local-tree records. Until the 300 missing tree-to-partition associations are resolved, we cannot determine whether any of the 3,623 local gene trees overlap the locus.

The five reported S-haplotype genes are CCM^T, GLO^T, CYP^T, PUM^T, and KFB^T. Hemizygosity is the documented recombination-suppression mechanism. We found no source establishing an inversion spanning this locus in the Stubbs sample panel. The S-locus tree is built from S-morph samples plus one outgroup, whereas the nuclear GF trees contain 106 accessions including both floral morphs. Homology/synteny across all sampled species is not established merely by reads mapped to the *P. veris* reference.

## Sampling scheme and tree limitations

Stubbs et al. report 50-kb nonoverlapping genomic fragments from 11 chromosomes, selecting every other GF and estimating 3,623 trees with IQ-TREE v2.1.2; ModelFinder used ascertainment-bias adjustment and the analyses used 1,000 ultrafast bootstrap replicates. The partition files instead enumerate 3,923 coordinate-bearing windows, and the deposit provides no explicit crosswalk for the 300-window mismatch. This prevents a defensible tree-to-coordinate mapping. Even when mapped, the 50-kb intervals and unsampled gaps must remain discontinuous in spatial plots.

The tree file supports downstream topology and branch-distance calculations in aggregate. Individual gene-tree branch lengths are fitted substitution distances and are not CU estimates. ASTRAL-III v5.7.3 was used by the paper for its species-tree estimate after collapsing branches with <10 bootstrap support.

## Exact source verification

`source_manifest.tsv` records Dryad URLs, listing sizes, local byte sizes, and local SHA-256 checksums for every requested file. Dryad does not advertise checksums in the listing, so the hashes are local provenance checks, not provider comparisons. The 19.59-MB tree file, 896-B ASTRAL tree, 2.56-KB S-locus tree, 2.03-KB README, and all 11 small partition files were used; no large FASTA or raw sequencing files were downloaded. `.DS_Store` is an operating-system artifact and was not analyzed.

## Stage 0 answers

- **3,623 deposited trees usable, with branch lengths?** Yes: all parse; 106 taxa; all non-root edges have numeric branch lengths in the IQ-TREE substitution/site scale.
- **Can every tree be mapped to coordinates?** No. 3,923 coordinate-bearing partitions versus 3,623 trees, with no crosswalk or tree IDs.
- **Which trees overlap the S-locus?** Unknown. Three partition windows fall in the published chr1 regional bracket, but their tree associations are unknown.
- **Can spatial Stage 1 proceed?** No. Resolving the 300-window mismatch requires the authors' original window-to-tree crosswalk, a deposited generation script/log, or an equivalent verified mapping source.
