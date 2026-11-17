# Sunflower Stage 1 spatial genealogical audit

## Decision

**NO-GO for claiming a defensible recombination-suppression-associated spatial genealogical effect or a sunflower MSRC success from these inputs.** The deposited Newicks support an arrangement-blind descriptive pilot, but the published regions are rounded to whole Mbp, observed windows are sparse across the very large blocks, flank/control replication is weak or unavailable, and linked windows are not independent replicates. Results below are descriptive; no significance claims are made.

## Data, definitions, and fixed quartets

The analysis re-parsed 3,242 coordinate-indexed Newicks from the local published tree table and joined called-base counts from the local window manifest. Each tree had 22 unique tips, finite nonnegative numeric branch lengths, and both fixed four-species sets. Two quartet sets were frozen: ANN (H. annuus), ARG (H. argophyllus), PET (H. petiolaris), plus either DEB (H. debilis) or PAR (H. paradoxus). All 48 observed individual combinations are scored separately within each species set and averaged to one window-level observation; they are not counted as independent loci. For both sets, q1 is the independently supported H. annuus + H. argophyllus sister pair; q2 and q3 are fixed alternatives ann-pet and arg-pet. The PAR set includes a hybrid-origin species; species histories in this radiation are reticulate and introgression is widespread, so q1 is a reference contrast, not an unambiguous bifurcating history [Owens et al. 2023](https://doi.org/10.1093/molbev/msad013).

Todesco et al. define the reported intervals as haploblock regions; their table notes that some regions comprise rearranged pieces. Support codes `g` mean recombination suppression in genetic maps, `G` means orientation differences between genetic maps, and `H`/`h` refer to Hi-C comparisons. These are not interchangeable: the table endpoints are not fine-mapped inversion breakpoints. The focal ann05.01 carries A,g (assembly differences plus map suppression), and ann13.01 h,g (reference-relative Hi-C plus map suppression) [Todesco et al. 2020, Extended Data Table 1](https://doi.org/10.1038/s41586-020-2467-6). Published integer-Mbp endpoints were treated as ±0.5-Mb uncertainty; the primary interior excludes 0.5 Mb at every edge. Arrangement genotypes for the exact tree tips remain unavailable, so interpretation is arrangement-blind.

## Coverage and coordinate audit

- **ann05.01:** 57 tree windows overlap the rounded interval; 57 lie wholly inside the conservative interior; 1 touch endpoint uncertainty bands; left/right 5-Mb flanks contain 8/10 windows. The union of observed tree intervals spans 1.435% of the reported physical interval.
- **ann13.01:** 97 tree windows overlap the rounded interval; 94 lie wholly inside the conservative interior; 3 touch endpoint uncertainty bands; left/right 5-Mb flanks contain 2/3 windows. The union of observed tree intervals spans 0.689% of the reported physical interval.

Validation found 3242 unique tree IDs, 1 distinct tip counts and 1 distinct branch counts. There were 2 overlapping tree-interval pairs across chromosomes. The manifest contains 3273 windows with ≥2,000 called bases; 31 of these lack a deposited Newick and are listed explicitly in `missing_tree_manifest_windows.tsv`.
Window spacing is highly uneven. Trees are discrete observations separated by unobserved genomic sequence; no interpolation is used. Called-base counts are retained per tree in the track table. The tree source's 3,242 records are used; no missing trees are inferred from the paper's 3,267 count.

## Spatial topology and branch lengths

- **ann05.01:** interior n=57, flanks n=8 left and 10 right. Mean q1/q2/q3 in the interior = 0.5205/0.2946/0.1849; mean D=q2−q1 = -0.2259. Mean branch summaries: A-like fixed species-pair distance contrast 0.0010 subs/site, B mean pairwise distance 0.0116 subs/site, and mean quartet internal branch 0.0021 subs/site.
- **ann13.01:** interior n=94, flanks n=2 left and 3 right. Mean q1/q2/q3 in the interior = 0.8063/0.0771/0.1166; mean D=q2−q1 = -0.7292. Mean branch summaries: A-like fixed species-pair distance contrast 0.0030 subs/site, B mean pairwise distance 0.0145 subs/site, and mean quartet internal branch 0.0033 subs/site.

Here `A_species_pair_distance_q2_minus_q1` is an explicitly arrangement-blind contrast of fixed species-pair distances; it is not the arrangement-defined A(w) used in MSRC analyses. `B_mean_pairwise_distance` is the mean patristic distance across the six pairs in each fixed four-tip combination, averaged over combinations. Quartet internal branches and all distances are substitutions/site from local ML gene trees, not time or coalescent units.

Endpoint sensitivity (0/0.5/1-Mb shrinkage) gives ann05.01 n=57/57/57 and mean D=-0.226/-0.226/-0.226; ann13.01 n=97/94/92 and mean D=-0.731/-0.729/-0.723. Equal-weighted nonempty 1/5/10-Mb bin means range -0.273 to -0.230 in ann05.01 and -0.767 to -0.719 in ann13.01; empty bins are omitted, not filled. These are sensitivity summaries, not independent replications.

Conservative adjacent-window transitions per observed adjacent Mb are listed separately. A run absent under the 20-kb adjacency threshold does not indicate a lack of change across intervening unobserved sequence.

Run detection treats two tree windows as adjacent only if the inter-tree gap is ≤20 kb. This conservative threshold makes almost all widely spaced tree observations isolated; an absence of adjacent runs is a lack of sampled adjacency, not evidence of genealogical continuity. Spatial autocorrelation is reported only with pair counts and is descriptive because pairs overlap and loci are linked.

## Controls and genome-wide quartet sensitivity

- **ann05.01 matched controls:** not_estimable; 0 eligible placements (20 required), matched on same chromosome and reported physical length with tree count within 20%, while excluding g-marked rounded intervals. No matched-control effect or null distribution is reported.
- **ann13.01 matched controls:** not_estimable; 0 eligible placements (20 required), matched on same chromosome and reported physical length with tree count within 20%, while excluding g-marked rounded intervals. No matched-control effect or null distribution is reported.

The fixed q1 fraction is 0.6495 with all 3,242 trees, 0.6518 after excluding ann05.01, 0.6446 after excluding ann13.01, and 0.6469 after excluding both. These small descriptive changes do not suggest a large impact on this selected quartet summary, but no uncertainty or independent matched-removal baseline is available. The exclusion analysis recomputes a fixed four-species quartet concordance summary; it is not a full ASTRAL species-tree inference, cannot estimate coalescent-unit branches, and counts linked windows only as input records. Because suppression-associated arrangement genotypes are unavailable and reticulation violates a simple MSC interpretation, a change in these descriptive frequencies is not called phylogenomic bias. Matched-control removals could not be constructed under the predeclared physical-length and same-chromosome constraints.

The ann13.01 q1 enrichment is similar across the two fourth-species choices: 0.806 with DEB versus 0.779 with PAR, compared with their respective genome-wide 0.650 and 0.550 baselines. Ann05.01 shows q1=0.520 (DEB set) and 0.488 (PAR set), with discordant alternatives carrying more support. The direction of ann13's q1 enrichment is therefore insensitive to this prespecified fourth-species choice, while ann05's specific discordant topology depends on which quartet set is used; neither result has matched controls.

## Answers to the Stage 1 questions

1. **Sustained topology shift?** No sustained alternative-topology replacement is established. Ann13.01 has q1=0.806 in its interior versus 0.650 genome-wide; ann05.01 is more discordant (q1=0.520, q2=0.295). Conservative adjacency yields only 6 and 7 adjacent interior tree pairs, respectively, and matched controls are unavailable, so these regional compositions cannot yet be attributed to recombination suppression.
2. **Branch lengths inside versus outside?** The report provides descriptive interior/flank distances in substitutions/site, but does not claim a robust regional difference from sparse flanks.
3. **Exceptional versus matched genomic controls?** Not estimable: eligible same-chromosome, same-length, similar-tree-count controls are insufficient.
4. **Independent suppression evidence?** Both focal rows are marked `g` for genetic-map recombination suppression. ann05.01 additionally has `A`; ann13.01 has `h`. The interval coordinates describe published haploblocks, not exact inversion breakpoints.
5. **Consistency across quartets?** Ann13 q1 enrichment is directionally consistent with either DEB or PAR as the fourth species; ann05 discordance is less stable in which alternative topology is favored. Both sets include all 48 sampled individual combinations averaged within each tree, not 96 independent replicates.
6. **Species-tree topology or CU lengths after exclusion?** Only fixed quartet concordance sensitivity is estimable here. No species-tree topology was re-inferred and CU branch lengths are unavailable.
7. **Arrangement-specific interpretation?** No. The exact 22 tree tips cannot be assigned to arrangements with available local metadata.
8. **Suppression, introgression, linked selection?** These data do not distinguish mechanisms. Introgression is independently documented in this annual sunflower radiation; suppression and linked selection remain plausible co-occurring processes.
9. **Strong enough plant dataset for MSRC?** No-go on current evidence. Reconsider only after exact boundaries, arrangement crosswalks, and adequately sampled matched controls become available.

## Outputs and next step

Coordinate and coverage audits, missing manifest windows, validation summaries, per-window quartet/distance tracks, region summaries, run lengths, transition counts, boundary/bin-size and quartet-selection sensitivity, distance-bin autocorrelation, control feasibility, and fixed-quartet exclusion sensitivity are stored in the adjacent TSVs. Six figures are provided as PNG/PDF/SVG and show only observed points; connecting segments appear only across gaps ≤20 kb. No raw inputs were downloaded or added to Git. The limiting next step is sample-to-arrangement crosswalk recovery plus finer independent interval definitions, not expansion to more haploblocks.
