# Primula Stage 1 analysis report

## Decision

**NO-GO for spatial MSRC inclusion from the deposited files as currently cross-referenced.** The gene trees themselves are suitable, but the deposit's coordinate partitions enumerate 3,923 windows for 3,623 unlabeled Newick records. No verified mapping connects the records to physical windows. The analysis stops before spatial topology/branch-length or regional inference; this is a data crosswalk limitation, not evidence for or against a localized S-locus effect. A topology-only comparison of the deposited S-locus tree with the species tree remains possible, as the paper already reports.

## Answers to the project questions

1. **Are the 3,623 trees usable, with branch lengths?** Yes. All 3,623 Newick records parse; all have the same 106 tips and complete numeric lengths on every non-root edge. They have unrooted trifurcating Newick roots and no non-root polytomies. IQ-TREE v2.1.2 DNA branch lengths are substitutions per nucleotide site, not coalescent units. Mean edge length is 0.0190781, median 0.00946494, range 0–1.631064; 40 edges are zero. See `tree_quality_summary.tsv`.
2. **Can every tree be mapped to genomic coordinates?** No. All 11 partition files contain 3,923 50-kb physical windows with alignment ranges, while the tree file contains 3,623 unlabeled records. The 300-window difference is unresolved. We do not assign tree indices by order.
3. **Which trees overlap the S-locus?** No local tree overlap can be verified. The *P. veris* S-locus is about 260 kb on chromosome 1 and hemizygous in the maternal/S-morph haplotype. Potente et al. place it within a putative pericentromeric region around 27.43–27.70 Mb; this is not an exact locus boundary. Three chr001 partition windows fall in that bracket, but the tree-to-window crosswalk is absent.
4. **Does the S-locus-associated phylogeny differ from genome-wide history?** Yes, as reported in the study: the S-locus and nuclear species-tree backbones show hard discordance, while shallow relationships are mostly concordant. Examples include *P. megaseifolia* and the *P. vulgaris*/*P. veris* relationships. ILS, S-locus introgression, and gene flow elsewhere are discussed as explanations; topology discordance alone does not distinguish them.
5. **Is that difference visible in the local-tree series?** Unknown. The local trees are usable, but their physical coordinates cannot be established from the deposit's current metadata.
6. **Is there a spatially coherent topology or branch-length shift?** Not determined. q1/q2/q3 spatial frequencies, D(w), run lengths, autocorrelation, A(w), B(w), and regional branch lengths were not calculated because tree-to-window identity is unresolved.
7. **Does recombination suppression coincide with a measurable genealogical signal?** Recombination suppression within the S-locus is biologically supported by hemizygosity. Coincidence with a measurable local-tree signal is not testable from the current crosswalk. No surrounding-chromosome shift or inversion is established here.
8. **Could this affect species-tree inference?** The paper's ASTRAL tree used all 3,623 GF trees, but we cannot identify the locus-overlapping trees for a leave-region-out comparison. The S-locus itself is a single linked unit; its published phylogeny difference is not evidence that it changed the genome-wide species-tree estimate.
9. **Does Primula satisfy MSRC empirical inclusion requirements?** No, pending a verified mapping. It remains a biologically interesting candidate, not a confirmed spatial MSRC example.

## Biological context and cross-system comparison

The S-morph (thrum) carries the S-locus; the L-morph (pin) lacks it. In *P. veris*, the S-locus is an approximately 260-kb maternal-haplotype region on chromosome 1 with five S-haplotype genes (CCM^T, GLO^T, CYP^T, PUM^T, KFB^T). Hemizygosity is the mechanism preventing recombination within the supergene. This differs from systems where recombination suppression is driven by a large inversion. The article mapped sampled genomes to the *P. veris* reference, but that does not independently establish full synteny of the locus across all taxa.

Existing MSRC comparisons provide design context only. *Neurospora* has a multi-megabase mating-type interval with physical local-tree tracks; *Anopheles* 2La is a large inversion with many mapped windows; fire-ant chromosome-16 analyses use local genealogies in a distinct supergene context; and *Microbotryum* comparisons rely more on evolutionary strata than continuous chromosome order. Primula's expected 260-kb locus and hemizygous architecture are biologically distinct. No cross-system figure is supported by the current Primula crosswalk.

## Required next input

The needed item is a verified crosswalk from each Newick record (1–3,623) to its chromosome and 50-kb interval, or the author scripts/logs that unambiguously reproduce this association and explain the 300 extra partition windows. The partition-site ranges alone do not identify which Newick record corresponds to which window. Once resolved, the S-locus overlap can be checked and Stage 1 resumed only if informative windows exist.

## Sources

- Stubbs et al. (2023), [paper](https://doi.org/10.1111/nph.18525) and [Dryad deposit](https://doi.org/10.5061/dryad.b2rbnzsgz).
- Potente et al. (2022), [*P. veris* chromosome-scale assembly and S-locus analysis](https://doi.org/10.1093/molbev/msac035); the paper reports the approximately 260-kb maternal-haplotype locus and chr1 pericentromeric location.
- IQ-TREE official [branch-length FAQ](https://www.iqtree.org/doc/Frequently-Asked-Questions): standard DNA-model branches are in substitutions per site.
