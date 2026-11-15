# Atlantic cod nucleotide Patterson's D: results and interpretation

## Data and estimator

The public Matschiner et al. (2022) Zenodo record 4560275 provides the four inversion SNAPP VCFs, four Dsuite VCFs and an outside-supergene SNAPP VCF. The analysis downloaded and verified the selected files by record MD5 and gzip integrity. All sample calls use the gadMor2 reference and the frozen linkage-group coordinates. The original feasibility audit's network blocker is superseded: the host was reachable after network permission was granted. No whole-genome reads or 3.36-GB outside Dsuite file were downloaded. See `data_feasibility_audit.md`, `nucleotide_source_inventory.tsv`, `sample_mapping_audit.tsv`, `vcf_sample_inventory.tsv`, and `download_validation.tsv`.

For each biallelic SNP within a published 250-kb tree window, diploid ALT frequencies in the three frozen ingroup populations produce ABBA and BABA contributions. The homozygous outgroup call defines the ancestral allele at that site; fixed-ALT outgroup calls are correctly repolarized, and missing or heterozygous outgroup genotypes are excluded. Population calls are complete in the primary SNAPP VCF sets. Regional D is `(sum ABBA - sum BABA)/(sum ABBA + sum BABA)` after pooling contributions, not an average of window D values.

## Primary results

| Linkage group | Inside D | Outside D | Inside − outside | Informative SNPs inside / outside | 1-Mb block-jackknife 95% CI for difference | Circular-shift p |
|---|---:|---:|---:|---:|---:|---:|
| LG01 | +0.9371 | −0.0362 | +0.9733 | 802 / 681 | [+0.887, +1.060] | 0.0175 |
| LG02 | +0.8282 | −0.0225 | +0.8507 | 297 / 1,245 | [+0.780, +0.922] | 0.0206 |
| LG07 | −0.9430 | −0.0107 | −0.9323 | 1,002 / 1,549 | [−0.989, −0.876] | 0.0238 |
| LG12 | +0.8496 | +0.0347 | +0.8149 | 564 / 928 | [+0.711, +0.919] | 0.1636 |

Across all matched windows (including windows that overlap a frozen boundary),
the primary VCF scan retained 19,560 / 16,490 / 25,613 / 20,295 biallelic SNP
records for LG01 / LG02 / LG07 / LG12. Homozygous-outgroup callable counts
were 19,178 / 16,167 / 25,075 / 19,852; outgroup missing/heterozygote
polarization failures were 382 / 323 / 538 / 443. ABBA/BABA-informative site
counts were 1,491 / 1,553 / 2,575 / 1,493. These totals are window-matched;
they exclude variants outside the exact published tree-window spans.

The estimates are extreme within each interval, and all four differ strongly in direction/magnitude from their same-linkage-group complements. LG01/LG02/LG12 favor ABBA for their frozen quartet orientation inside the interval; LG07 favors BABA. The physical circular-shift validation is nominally unusual for LG01, LG02 and LG07, while LG12 is not. This is descriptive with only 113, 96, 125 and 109 valid circular shifts, respectively, and should not be read as a multi-linkage-group causal test. The block-jackknife intervals address spatial clustering but do not remove population-history confounding.

The LG12/Bornholm-associated primary design has inside D=+0.850 and outside D=+0.035. Its nearest 1-Mb left and right flanks have only 2 and 11 informative sites, respectively, so their D values (−1.00 and +0.186) are too sparse to establish a sharp local transition or distinguish Bornholm-specific history from the broader interval signal.

The equal-direction statement is **not** true across all four groups: LG07 has the opposite sign. That signal is stable to transversions, reducing each ingroup to one individual, the alternative Gadoga outgroup, and the denser Dsuite inversion VCF at a 50% per-population call threshold. Inside-minus-outside D ranges are +0.986 to +1.000 (LG01), +0.849 to +0.881 (LG02), −0.969 to −0.871 (LG07), and +0.806 to +0.839 (LG12). The outside estimates remain the same complete-call SNAPP background in these contrasts; Dsuite is an inside-region variant/callability sensitivity, not a fully matched Dsuite outside comparison. See `robustness_analysis.tsv`.

Only the primary population contrast per LG was estimated, so robustness to
alternative population triples is not established. Geography/ecotype remains
confounded, particularly because the LG07 P1/P3 populations share a migratory
ecotype. Gadcha and Gadmac were inventoried but not scanned as additional
outgroups. Reference bias and SNP ascertainment were not separately quantified
beyond comparing the deposited SNAPP and Dsuite variant sets and using
transversions-only sensitivity.

Boundary shifts of ±50 and ±100 kb, 50/100/250/500-kb aggregations, per-window records, and 250/500/1,000-kb block deletion results are supplied as separate tables. D is pooled from site contributions at each aggregation; no window D average is used. The inversion definitions remain those in the frozen repository manifest.

## Spatial relation to published trees

The 426 published 250-kb SNAPP population trees have the 12 cod-population tips but do not include the VCF outgroup. Consequently, a four-taxon induced quartet and exact `D_topology=q3-q2` are **not identifiable** from those trees. The matched output reports the rooted three-population sister-pair contrast separately and never calls it exact D_topology. Per-window nucleotide D correlates descriptively with that rooted-triple contrast (Pearson r: LG01 +0.717, LG02 +0.524, LG07 +0.505, LG12 +0.528; 1-Mb block-bootstrap intervals are in `spatial_correlations.tsv`). D versus the existing arrangement-minus-baseline signal is strong on LG01 (+0.887), LG02 (+0.646), and LG12 (+0.692), and negligible on LG07 (+0.071). This recreates the expected broad distinction between the large prior arrangement topology shifts and weak LG07 arrangement contrast, but D on LG07 is itself extremely strong and negative. The allele-sharing signal therefore does not simply reproduce arrangement concordance.

Correlation with the available SNAPP divergence-time contrast is positive on LG01 (+0.842), LG02 (+0.816), and LG12 (+0.805), and negative on LG07 (−0.842). SNAPP branch lengths are population divergence-time units; they are not substitution-unit gene-tree lengths or coalescent-unit branch lengths. These associations are spatially descriptive and not independent-window significance tests.

## Interpretation and comparison

The defensible conclusion is that the four inversion intervals contain strong, spatially concentrated nucleotide site-pattern asymmetry for these pre-frozen population quartets, and that the sign is quartet-specific. An introgressed or otherwise alternative genealogy is one possible explanation, alongside incomplete lineage sorting, ancestral structure, demographic admixture, linked selection, population structure, SNP ascertainment, and polarization or reference effects. The estimates do not establish introgression uniquely and do not demonstrate that a species-tree estimator changes topology or coalescent-unit branch lengths. An exact four-taxon genealogical comparison is blocked by the missing outgroup tip in the published window trees; the separate three-population concordance is not an estimator of species-tree bias.

Relative to the repository's other empirical examples, the cod shifts are much more extreme than the prior Neoaves chromosome-4 results, where Columbea and N62 showed essentially no inside/outside change (−0.0004 and +0.0009), while N61 shifted by −0.1776. For *Anopheles gambiae* 2La, the primary inside D values are approximately −0.66 to −0.73 for A/C-designs and +0.81 to +0.84 for B/D-designs, against near-zero outside values: the interval-specific magnitude and design-dependent sign are qualitatively similar to cod. For the house-mouse t-complex corrected quartet analyses, inversion-union D was +0.24 to +0.39 against +0.014 to +0.044 on chromosome complements, a more moderate but clear spatial contrast. These datasets differ in taxa, ploidy/call filtering, physical unit, and quartet definitions, so this is contextual comparison, not a pooled cross-system effect estimate.

**Recommendation: CONDITIONAL GO** for a manuscript result describing a strong and directionally heterogeneous cod nucleotide site-pattern signal across inversion regions. Before making a causal claim or asserting an estimator bias, prioritize fully matched outgroup-supported local genealogies, more independent individuals per outgroup, demographic/geographic controls, and an explicit species-tree/coalescent-unit estimation sensitivity using these genotypes. The measured D is an allele-sharing statistic, not a unique signature of MSRC or introgression.
