# Patterson's D robustness validation

## Scope and estimator

This is a direct genotype-level validation of the Atlantic cod Stage-4 Patterson's D analysis. It does not infer nucleotide patterns from local trees. All population quartets and sample IDs were frozen before the additional scans. The 12 added designs are deterministic alternatives drawn from Stage-3 outside-inversion baseline sister relationships, with two individuals per ingroup population and the predeclared Borsai outgroup. The designs are in `../config/frozen_validation_quartets.tsv`; their source candidate, baseline split, arrangement states, geography/ecotype, and rationale are recorded there.

The direct scanner recalculates diploid allele frequencies and ABBA/BABA contributions at each PASS (or `.`) biallelic SNP. It requires complete calls in the six frozen ingroup individuals and a fixed homozygous outgroup. Fixed-ALT outgroups are recoded correctly; heterozygous or missing outgroups are excluded. D is calculated from pooled ABBA and BABA sums, not averaged across windows.

## Primary result reproducibility

The independent parser exactly reproduces the previous primary SNAPP results to floating-point precision. The absolute difference from the prior calculation is below 3e-11 for all inside and outside estimates. The inside estimates are LG01 +0.9371, LG02 +0.8282, LG07 −0.9430, and LG12 +0.8496. Same-linkage-group outside estimates are −0.0362, −0.0225, −0.0107, and +0.0347, respectively. The inside-minus-outside contrasts are +0.9733, +0.8507, −0.9323, and +0.8149.

| Linkage group | Source | Inside D | Outside D | Inside − outside | Polarized sites inside / outside |
|---|---|---:|---:|---:|---:|
| LG01 | SNAPP complete-call | +0.9371 | −0.0362 | +0.9733 | 11,110 / 7,886 |
| LG01 | Dsuite missing-call-allowed | +0.9361 | −0.0232 | +0.9593 | 103,909 / 72,905 |
| LG02 | SNAPP complete-call | +0.8282 | −0.0225 | +0.8507 | 2,998 / 13,019 |
| LG02 | Dsuite missing-call-allowed | +0.8461 | −0.0260 | +0.8721 | 27,806 / 123,305 |
| LG07 | SNAPP complete-call | −0.9430 | −0.0107 | −0.9323 | 9,767 / 15,024 |
| LG07 | Dsuite missing-call-allowed | −0.9371 | +0.0013 | −0.9384 | 79,205 / 136,810 |
| LG12 | SNAPP complete-call | +0.8496 | +0.0347 | +0.8149 | 10,015 / 9,762 |
| LG12 | Dsuite missing-call-allowed | +0.8584 | +0.0214 | +0.8370 | 151,010 / 156,338 |

Boundary-overlapping published 250-kb windows were excluded (one on LG01, LG02, and LG12; two on LG07). Inside and outside are same-linkage-group comparisons, not neutral controls. The Dsuite `.full` data retain missing calls at the cohort level but the estimator still requires all six selected ingroup calls and a fixed outgroup at each site.

The missing-call-allowed source increases polarized sites by roughly 8–16-fold inside and 9–13-fold outside, while preserving each primary inside direction and the large absolute inside/outside contrast. Its inside/outside contrasts differ from SNAPP by −0.0139 (LG01), +0.0213 (LG02), −0.0061 (LG07), and +0.0220 (LG12). This is strong robustness to the complete-call site filter, not an independent biological replication. Every SNAPP regional position is also found in its corresponding Dsuite regional VCF; the additional Dsuite records reflect the less restrictive missing-call rule. The Dsuite roster lacks Gadmor_lfo2, so the frozen two-individual LFO sample is represented by Gadmor_lfo1 plus Gadmor_lfo3 in LG01/LG02 (P2) and LG07 (P1); LG12 is unchanged. The Dsuite contrast therefore combines a callability change with this explicit one-individual substitution. Exact sample IDs are listed in `source_sample_mapping.tsv`. The cohort-wide minor-allele-count filter and upstream SNP quality/masking pipeline remain shared ascertainment limits.

## Outgroup and population sensitivity

All four verified external taxa yield the same inside sign for each linkage group. Across Borsai, Gadoga, Gadcha, and Gadmac, inside D spans +0.932 to +0.937 on LG01, +0.828 to +0.862 on LG02, −0.959 to −0.943 on LG07, and +0.846 to +0.855 on LG12. Outside estimates remain close to zero (LG01 −0.072 to −0.036; LG02 −0.023 to −0.017; LG07 −0.017 to −0.008; LG12 +0.016 to +0.035). Sitewise outgroup heterozygosity and missingness are explicitly tallied in `outgroup_sensitivity.tsv`; polarization uses only fixed homozygous outgroup genotypes.

The three alternative frozen population designs show a coherent but nonuniform pattern. The alternative q3 arrangement-sharing design gives strong same-direction inside D on LG01 (+0.935), LG02 (+0.799), and LG07 (+0.980), while LG12 is near zero and negative (−0.066). The alternative q2 arrangement-sharing design is negative on LG01 (−0.880), LG02 (−0.941), and LG07 (−0.083), but near zero on LG12 (−0.004). Baseline-pair controls are smaller in absolute magnitude (0.337, −0.064, −0.161, and −0.012). These are one quartet per class and linkage group; shared populations mean these are descriptive sensitivity checks, not independent replicates or a multiplicity-tested survey.

The LG07 primary signal is large and robust across source and outgroup, despite the earlier Stage-4B arrangement-concordance contrast being approximately zero (−0.006). This mismatch is scientifically informative: the nucleotide asymmetry is not a mechanical restatement of the earlier arrangement topology statistic. It could reflect population history or other sources of gene-tree/site-pattern discordance.

## Spatial uncertainty and remaining limitations

For the primary design, block-jackknife uncertainty for inside-minus-outside D is stable across 250-kb, 500-kb, and 1-Mb physical block sizes. At 1 Mb, delete-block SEs are 0.044 (LG01), 0.036 (LG02), 0.029 (LG07), and 0.053 (LG12) for SNAPP, and 0.028, 0.027, 0.011, and 0.028 for Dsuite. All primary contrast intervals remain away from zero in both sources at these resolutions. These intervals summarize genomic block variation under this design; they do not correct demographic or SNP-ascertainment biases.

The Relate source is not a valid external-outgroup D dataset: it contains only 24 Atlantic cod individuals, no Borsai/Gadoga/Gadcha/Gadmac, and its deposited colinear VCF has no records on LG01/LG02/LG07/LG12. The Relate within-region site sets overlap the Dsuite VCFs but cannot independently polarize the four-population patterns. No less-ascertained, outgroup-bearing whole-genome callset was identified in the local/public deposits used here; raw reads were not downloaded.

Previously reported spatial comparison against population-tree features remains relevant: the published SNAPP trees lack the external outgroup, so an exact four-tip q1/q2/q3 or exact D_topology cannot be derived. Earlier arrangement-associated tree shifts correlate with D spatially on LG01 (+0.887), LG02 (+0.646), and LG12 (+0.692), but not LG07 (+0.071). The LG07 D result therefore reinforces the need to keep site-pattern asymmetry separate from topology-frequency contrasts. Neither nonzero D nor spatial coincidence uniquely establishes introgression or supports MSRC over ancestral structure, incomplete lineage sorting, or other demographic alternatives.

## Conclusion

The validation supports a CONDITIONAL GO for describing strong, spatially localized allele-sharing asymmetry across the four inversion regions. The magnitudes and signs are stable to direct estimator reimplementation, four external outgroups, and the missing-call-allowed Dsuite source. LG07 is strongly negative despite its weak prior arrangement-concordance signal, and LG12's alternative population designs are near zero. Results should be framed as polarized site-pattern asymmetry with population-geographic and ascertainment caveats, not as unique proof of introgression, species-tree estimation bias, or MSRC.
