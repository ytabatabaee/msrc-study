# Neurospora nucleotide Patterson's D: results and interpretation

## Decision

**CONDITIONAL GO** for reporting a rooted, nucleotide-based spatial D analysis on the deposited L10 variable-site panel. A checksum-matched *N. discreta* alignment supplied locally makes outgroup polarization and coordinate mapping possible. This remains conditional because only the filtered 50-SNP-per-window L10 matrix was analyzed, four mating-type labels were excluded, and the published local trees do not include *N. discreta*.

## Primary estimate

The frozen quartet was `((L10 mat-a, L10 mat-A), N. hispaniola), N. discreta)`. Verified mat-a group: 10 isolates; verified mat-A group: 8; *N. hispaniola* P3: `8817@999`; *N. discreta* O: `Ndis@999`. Per-site allele-frequency ABBA/BABA contributions used haploid calls, a 50% minimum called fraction in each multi-isolate ingroup, a called singleton P3, and an unambiguous outgroup matching one of exactly two ingroup alleles. Fixed-ALT outgroup sites are properly polarized; gaps, ambiguous bases and outgroups outside the two ingroup alleles are excluded. The primary quartet and sample membership were frozen before the D scan.

Across 2,265,369 variable-site rows, 646,088 were polarizable and 110,775 informative. Pooled all-genome D was **+0.5497** (ABBA 52,851.51; BABA 15,356.64). This pooled result is dominated by LGI and should not be read as a genome-wide estimate from a random set of callable bases.

On LGI, the frozen interval 1.5–8.0 Mb has D **+0.8712** (ABBA 37,024.50; BABA 2,547.90; 39,726 informative sites). Outside it, D is **+0.4460** (ABBA 3,698.94; BABA 1,417.32; 8,822 informative sites). Inside-minus-outside is **+0.4253**. Genomic block jackknife 95% intervals for the contrast are +0.084 to +0.766 (250-kb; 12 outside blocks), +0.126 to +0.725 (500-kb; 6 outside blocks), and +0.185 to +0.666 (1-Mb; 4 outside blocks). The outside estimate itself is imprecise (250-kb block interval +0.106 to +0.786; 12 blocks). A circular-shift test of the 100-kb spatial track gives p=0.0105 (95 nonidentity shifts; fixed tract length/boundaries). Boundary sensitivity using 1.2–8.2 Mb retains a strong inside/outside difference (+0.8179; see summary table). The 250-kb left and right flanks have D +0.7191 and +0.7869, respectively, so the elevated asymmetry is not confined to the inversion interior.

The two documented opposite-mating-type pairs combined yield LGI D +0.8748 inside and +0.4809 outside. Transversion-only D is +0.8709 inside and +0.4719 outside. These sensitivities preserve direction and magnitude; they do not resolve the high outside baseline. Remaining linkage groups' pooled D are LGII +0.0229, LGIII −0.0041, LGIV +0.0715, LGV −0.0049, LGVI −0.0016, and LGVII +0.0069. The filtered variable-site matrix is enriched for variants and includes donor taxa; these are not neutral genome-wide baselines.

## Genealogy comparison and limits

At matched 100-kb LGI bins, nucleotide D is descriptively associated with published local-tree TWISST q3 (Pearson r=+0.949, Spearman rho=+0.796, 93 bins) and inversely with q1 (r=−0.940, rho=−0.745). These q weights are not nucleotide pattern probabilities. More critically, the tree set lacks *N. discreta*, so it cannot supply the exact rooted quartet splits needed for `D_topology=q3−q2`; no exact D-vs-D_topology comparison is claimed. Available branch metrics are published substitution-distance/divergence proxies, not coalescent-unit gene-tree branches. Correlations with those metrics are descriptive, linked-bin associations only.

The very high D outside the frozen introgressed interval and high values in both flanks indicate that an interval-only “introgression causes D” interpretation is not justified. Possible contributors include other linked genealogy, the SNP ascertainment and filtering scheme, population history/geography, reference/alignment effects, and incomplete sampling. The D statistic detects asymmetry in derived-allele sharing under the chosen rooted quartet; it is not unique evidence for introgression and does not uniquely support MSRC.

## Scope and next steps

This is a L10 analysis, not an analysis across all seven Neurospora lineages. CJ07, CJ11, CJ12, and CJ18 remain excluded because no defensible mating-type assignments were found. Per-window tables retain site/filter counts. The outgroup sequence is aligned to the 2509 reference coordinate axis; all 2,265,369 genotype reference alleles match the corresponding alignment column exactly. Alignment gaps/ambiguities are substantial and removed sitewise. No raw reads were downloaded or processed.

Before manuscript-level causal claims, an independent callable/ascertainment-aware variant panel, additional independently rooted population contrasts, and exact matched outgroup-bearing local trees (or an explicitly separate genealogy estimator) would be needed. Reproduce from repository root with `python3 empirical/neurospora/analysis/patterson_d/scripts/run_patterson_d.py`, then `python3 empirical/neurospora/analysis/patterson_d/scripts/make_figures.py`.
