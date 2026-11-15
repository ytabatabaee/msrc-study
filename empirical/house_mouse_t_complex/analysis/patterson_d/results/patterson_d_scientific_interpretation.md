# House mouse chromosome 17 nucleotide Patterson D — interpretation

## Feasibility and frozen population designs

**Stage 1 passed (READY_TO_RUN), and Stage 2 completed.** The local Kelemen–Vicoso archive contains 4,046 inferred ML gene trees and position-only annotations, not genotypes. The public Harr et al. AllMouse PASS VCF was verified as a tabix-indexed nucleotide source and queried remotely without downloading its 55.7-GB compressed file. It contains 67 samples and spans mm10/GRCm38 chr17:1–94,987,271; all 55 selected local tree-tip sample IDs are present. The complete nucleotide scan covers 18,998 adjacent 5-kb bins across chr17. The inversion-focused segment is 5.0–42.05 Mb; exact tree overlap is 5.0–37.055 Mb.

Archived sample records identify 15 heterozygous t-designated mice (Mmd 4, Mmm 8, Mmc 3), 32 standard mice, and 8 *M. spretus*. Three rooted contrasts were fixed in the runner before D calculation and are explicitly tabulated in `config/frozen_quartet_members.tsv`: ((Mmd_t,Mmd_std),Mmm_std),O; ((Mmm_t,Mmm_std),Mmd_std),O; and ((Mmc_t,Mmc_std),Mmm_std),O. *M. spretus* is a close, separate outgroup used by the source study. These are exploratory arrangement-associated population contrasts rather than interchangeable conventional species-tree null quartets. The VCF is unphased, so carrier diploids do not identify the t-bearing chromosome.

## Inversion-union contrasts

For each biallelic single-base A/C/G/T SNP, the population-frequency estimator is ABBA=(1-p1)p2p3 and BABA=p1(1-p2)p3, where p1–p3 are ALT frequencies and the eight outgroup samples must all be called and fixed for REF. ABBA/BABA are weighted contributions; D is calculated from pooled contributions. The VCF is SNP-only, so invariant callable-base denominators and mutation-rate interpretations are unavailable. For inversion contrasts, a 5-kb bin is assigned inside if it overlaps the published interval; the inversion-union uses the union of the four validated intervals. “Outside” below means the full chr17 complement, not a presumed neutral control. Estimates use 500-kb physical-block delete-one jackknives.

- **D1_Mmd_t_vs_standard:** inversion-union D=-0.2467 (ABBA 3133.1, BABA 5185.2; 55626 informative SNPs), full-chromosome complement D=-0.0670 (ABBA 7123.2, BABA 8146.4; 84231 informative SNPs). Signed difference=-0.1797 (500-kb block-JK SE 0.0502); absolute-D difference=+0.1797 (SE 0.0502). 190 physical blocks overall; informative in/out blocks 68/124.
- **D2_Mmm_t_vs_standard:** inversion-union D=+0.0375 (ABBA 2555.6, BABA 2370.7; 52910 informative SNPs), full-chromosome complement D=+0.0149 (ABBA 2753.4, BABA 2672.4; 71649 informative SNPs). Signed difference=+0.0226 (500-kb block-JK SE 0.0608); absolute-D difference=+0.0226 (SE 0.0608). 190 physical blocks overall; informative in/out blocks 68/124.
- **D3_Mmc_t_vs_standard:** inversion-union D=+0.2142 (ABBA 10627.9, BABA 6878.3; 99805 informative SNPs), full-chromosome complement D=+0.0095 (ABBA 29036.0, BABA 28490.1; 260557 informative SNPs). Signed difference=+0.2047 (500-kb block-JK SE 0.0245); absolute-D difference=+0.2047 (SE 0.0245). 190 physical blocks overall; informative in/out blocks 68/124.

By absolute D, asymmetry is stronger inside the inversion union for D1 (+0.180, SE 0.050) and D3 (+0.205, SE 0.024), while D2 is close to zero (+0.023, SE 0.061). D1 and D3 have opposite signed interpretations because the arrangement/population design differs. These are exploratory effects; the 500-kb SEs do not correct for design/interval multiplicity or population confounding.

## Individual inversions, collinear gaps, and flanks

Individual interval estimates are provided in `patterson_d_jackknife.tsv` at 100-kb, 250-kb, 500-kb, and 1-Mb block sizes. D1 shows larger absolute D inside In(17)1–4 (+0.188, +0.246, +0.288, +0.059). D2 shows increases in In(17)1–2 (+0.331, +0.297), but little change in In(17)3–4 (+0.046, −0.019). D3 has positive absolute-D differences in all four intervals (+0.256, +0.094, +0.190, +0.194). The comparisons are design-dependent and individual proximal intervals have only three informative 500-kb blocks.

D is also elevated in collinear gap 3 (chr17:18.05–19.40 Mb; D1 +0.615, D2 +0.732, D3 +0.510), based on only three informative 500-kb blocks and with high uncertainty for some designs. Other gaps and flanks vary by contrast. Thus asymmetry is not confined to inversion spans, and inversion-versus-complement differences cannot establish recombination suppression as its cause. Each segmented interval and its full-contig complement are tabulated.

## Gene trees and other project comparisons

Exact-start matching yields 4,046 nucleotide/tree windows. For ((P1,P2),P3),O, ABBA reflects derived sharing by P2+P3 and BABA by P1+P3. The existing full-tip q_t_alt/q_species and tree D tracks do not use these exact population subsets, so they are not equivalent estimands. Across 58 matched nonempty 500-kb blocks, D1 nucleotide D versus mean q_t_alt has descriptive Pearson r=-0.097; overall co-location is weak. Branch-distance and recombination classifications are joined by exact window start. No formal spatial permutation test was run.

Compared qualitatively with the local Neoaves report, where the interval effect was strong in one design but nearly absent in two, house mouse effects also vary by quartet and interval; here D1/D3 show stronger absolute D in the union while D2 does not, and collinear gap 3 also has high D. The local Anopheles 2La report found large, consistent inversion-localized effects across four merus designs, unlike this design-dependent mouse pattern. These studies differ in samples, alignment/VCF ascertainment, and quartet estimands; numerical effect sizes are not directly comparable.

## Robustness and interpretation

Resolution pooling is reported at 5, 50, 100, 250, and 500 kb from summed ABBA/BABA contributions; empty and low-information bins remain present. Synthetic tests use the exact estimator and outgroup-polarization helpers and pass for pure ABBA, pure BABA, balanced, uninformative, missing-outgroup, and derived-outgroup examples. Structural checks passed for 56,994 5-kb rows, D bounds, 50-kb count conservation, 4,046 exact tree matches, and 15 nonempty PNG/PDF/SVG exports.

Transversion-only, alternative-outgroup, and alternate SNP-filter scans were not run. Partial ingroup missing calls are omitted from population-frequency means; sites with an entirely missing group or any missing outgroup are tallied. Per-sample call rates and invariant callability are not available from this SNP-only export. Geography/subspecies structure, admixture or introgression, ancestral population structure/polymorphism, incomplete lineage sorting, SNP ascertainment, and polarization remain alternatives. Nonzero D does not uniquely support introgression or MSRC.

No formal test establishes that the observed nucleotide asymmetry is caused by suppressed recombination or matches genealogy shifts causally. See the window, multiresolution, inside/outside, jackknife, topology, robustness, and figure outputs for complete records.
