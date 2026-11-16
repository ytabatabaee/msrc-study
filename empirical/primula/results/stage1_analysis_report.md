# Primula Stage 1 analysis report

## Decision

**NO-GO for Primula as a confirmed MSRC empirical example on the evidence currently available.** Conditional resumption remains reasonable after the listed small Dryad data files can be downloaded. The required primary inputs were unavailable locally because Dryad download endpoints returned HTTP 403. Stage 1 statistics and figures were not produced. This is an acquisition and verification failure, not evidence that the locus lacks a signal.

## Answers to the project questions

1. **Are the 3,623 trees usable, with branch lengths?** The article reports 3,623 IQ-TREE gene trees from 3,623 50-kb genomic-fragment alignments, but payload inspection was blocked. Parseability, exact record count, taxon coverage, branch-length completeness, and branch-length units remain unverified. IQ-TREE branch lengths would ordinarily be in the fitted substitution-distance scale; this cannot be assumed for these files without inspecting the records and deposited README.
2. **Can every tree be mapped to genomic coordinates?** Not yet established. The deposit lists 11 chromosome partition files, but their contents were not retrievable. No physical positions were inferred from tree order.
3. **Which trees overlap the S-locus?** Unknown. The study mapped resequencing reads to the chromosome-scale haploid *P. veris* assembly described by Potente et al. (2022), and the Dryad README describes complete nuclear and S-locus alignments. The overlap between the sampled GF tree intervals and S-locus coordinates still requires an explicit coordinate map and reference-matched locus annotation.
4. **Does the S-locus phylogeny differ from genome-wide history?** Yes, the paper reports hard incongruence at deep backbone relationships, while shallow relationships are mostly concordant. Examples include the relationships of *P. megaseifolia*, *P. vulgaris*, and *P. veris*. The paper discusses ILS and S-locus introgression as possible explanations; the topology difference alone does not establish either mechanism.
5. **Is that difference visible in the local-tree series?** Unknown, because the local tree records could not be inspected and their genomic positions were not recovered.
6. **Is there a spatially coherent topology or branch-length shift?** Not established. No q1/q2/q3 frequencies, D(w), run lengths, autocorrelation, patristic contrasts, or internal branch statistics were calculated.
7. **Does recombination suppression coincide with a measurable genealogical signal?** Not established. A well-characterized *P. vulgaris* S-locus is a 278-kb thrum-specific hemizygous interval with five genes; this is evidence of hemizygosity and suppressed recombination within the S region, not a measurement of a wider surrounding genealogy shift in the *P. veris* local-tree series. The published 278-kb value must not be substituted for a verified *P. veris* interval.
8. **Could this affect species-tree inference?** In principle a small region can affect a genome-wide summary only in proportion to its sampled weight and linkage structure. Here neither the number of overlapping GF trees nor their quartet assignments is known, so influence cannot be quantified. The reported S-locus discordance should not be conflated with demonstrated genome-wide species-tree sensitivity.
9. **Does Primula satisfy MSRC empirical inclusion requirements?** No, not yet. Spatial coordinates, local-tree overlap, gene-tree statistics, branch metrics, and sensitivity analyses are prerequisites for inclusion as a regional MSRC example.

## Verified study design and biological context

Stubbs et al. state that they partitioned 11 chromosomes in the *P. veris* reference assembly into nonoverlapping 50-kb fragments, selected every other fragment, and estimated trees from the resulting 3,623 alignments using IQ-TREE 2.1.2. They used ModelFinder with ascertainment-bias adjustment and 1,000 ultrafast bootstrap replicates. The described sampling leaves intervening 50-kb fragments unsampled, so a plotted series must contain gaps and avoid interpolating across them. Their ASTRAL-III species tree used all 3,623 GF trees (ASTRAL-III v5.7.3).

The article reports strong backbone discordance between the S-locus and nuclear genome phylogenies. It states the S-locus topology places *P. vulgaris* as sister to *P. elatior* plus Caucasus species (excluding *P. el. pseudoelatior*), while the nuclear topology places *P. vulgaris* with *P. juliae* and Caucasus taxa, then *P. veris*. These are published whole-locus versus whole-nuclear-tree comparisons; neither establishes a spatial transition in surrounding local trees.

The floral polymorphism is conventionally described as L-morph/pin and S-morph/thrum. Comparative genomic work in *P. vulgaris* describes a roughly 278-kb region containing five S-haplotype genes (CCM^T, GLO^T, CYP^T, PUM^T, KFB^T), hemizygous in thrum plants and absent in pins. This structural hemizygosity provides a mechanism for low/no recombination within the region. Whether the same coordinate boundaries and structural state apply across every sampled species must be confirmed from the *P. veris*-referenced alignments and source data; shared gene labels or alignment to one reference do not alone prove conserved synteny across all taxa. The retrieved sources do not establish a specific inversion in the studied panel.

## Cross-system comparison (interpretive context only)

- **Neurospora:** existing MSRC work uses physically ordered local genealogies across a large mating-type suppressed-recombination interval and flanks. Its multi-megabase scale offers many more local windows than a 278-kb plant locus could provide under 100-kb start-to-start GF sampling.
- **Anopheles 2La:** a large inversion provides a coordinate-defined inside/outside comparison with many windows; its suppressed recombination is tied to an inversion and should not be assumed to match Primula hemizygosity.
- **Fire ants:** the existing chromosome-16 analyses use multiple local gene trees and supergene/mating-type-associated states; structural and taxonomic context differs from Primula.
- **Microbotryum:** existing local trees support comparisons among evolutionary strata, but physical ordering is less complete; strata-level contrast is not equivalent to a continuous neighborhood track.

All comparisons above are design context, not a claim that results or mechanisms transfer across taxa. Primula's expected narrow footprint, hemizygosity, potentially ancient persistence, and whole-locus phylogenetic discordance make it a biologically distinct candidate.

## Required next steps

1. Retrieve only the listed small tree, README, species-tree, S-locus-tree, and partition files.
2. Reconcile all 3,623 Newick records with partition metadata and independently validate chromosome, start, end, and 50-kb gap structure.
3. Obtain a reference-matched *P. veris* S-locus coordinate annotation and document assembly accession/version; define flanks before looking at local topology results.
4. Only if mapping and overlap are adequate, freeze all biologically supported quartets, compute local topology and branch metrics, then conduct region and ASTRAL sensitivity analyses.
5. Generate figures only from verified mapped observations; visually inspect PNG/PDF/SVG outputs.

## Sources

- Stubbs et al. (2023), [paper](https://doi.org/10.1111/nph.18525), methods and results.
- [Dryad deposit](https://doi.org/10.5061/dryad.b2rbnzsgz), published 2022-10-13.
- Potente et al. (2022), [Comparative Genomics Elucidates the Origin of a Supergene Controlling Floral Heteromorphism](https://doi.org/10.1093/molbev/msac035), *P. veris* assembly and comparative S-locus analysis.
- Li et al. (2016), [Genetic architecture and evolution of the S locus supergene in *Primula vulgaris*](https://doi.org/10.1038/nplants.2016.188).
