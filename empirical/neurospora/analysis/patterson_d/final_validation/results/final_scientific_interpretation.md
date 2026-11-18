# Final scientific interpretation: Neurospora L10 nucleotide D

## Decision: CONDITIONAL GO

The analysis is suitable for a carefully qualified manuscript result about allele-sharing asymmetry **conditional on the deposited, ascertained variable-site panel**. It is not an unbiased genome-wide callable-site analysis. The rooted nucleotide estimator reproduces a strong LGI contrast in the verified L10 samples, but the outside estimate is itself strongly positive and the immediate flanks are high. The signal therefore is not confined to the published 1.5–8.0 Mb interval. The data do not uniquely attribute the pattern to the MSRC, recombination suppression, introgression, or any one mechanism.

## Four claims kept separate

**H1 — D differs from zero.** Pooled D is strongly positive inside (+0.87123) and positive in the LGI complement (+0.44595) in the deposited variant panel. Physical-block intervals for the outside region are broad and based on 12/6/4 blocks at 250 kb/500 kb/1 Mb. Because P1 and P2 are mating-type sample groups from one recombining L10 population outside the suppressed tract, “D=0 under no gene flow” requires additional assumptions about demographic structure, exchangeability, and ascertainment. The nucleotide D estimate is valid as the specified sample contrast; a simple species-tree null is conditional.

**H2 — D changes spatially along LGI.** The pooled inside-minus-complement contrast is +0.42527. Existing 250/500/1000-kb physical-block intervals are +0.084 to +0.766, +0.126 to +0.725, and +0.185 to +0.666. The contrast persists under common support in source-window span and raw callability. But the flanks are elevated: at 250 kb, left/right flank D is +0.719/+0.787; at 500 kb, +0.700/+0.602; at 1 Mb, +0.574/+0.455. D declines toward distant sequence. The focal interval is the strongest part of a broader spatial profile, not an isolated island with an obvious neutral baseline.

**H3 — D corresponds spatially with published genealogy weights.** Matched 100-kb D and TWISST q3 have r=+0.949 across 93 bins, +0.454 inside, and +0.843 outside. After demeaning within focal/complement classes, r=+0.785 (500-kb block-bootstrap 95% interval +0.574 to +0.916). This shows some association beyond the difference in regional means. The trees lack Ndis, so this is not exact rooted nucleotide `D_topology=q3−q2`; it is a descriptive association between distinct estimands from related selected data.

**H4 — The joint pattern uniquely identifies MSRC.** Not established. Prior work documents introgression involving L10 and *N. hispaniola*, while this dataset also has plausible ancestral structure, within-lineage/population effects, linked genealogy heterogeneity, branch-length/substitution effects, and SNP/outgroup ascertainment. The current D/topology analyses do not separate these causes and do not prove either that introgression occurred in every high-D segment or that no gene flow occurred.

## Computational validation and contradictions

An independent scalar reference counter exactly reproduced all 265 deterministic raw-site checks spanning seven linkage groups, focal/flank/outside regions, outgroup states, and low/high group-callability strata. No status, ABBA, or BABA discrepancy was observed. The original coordinate map was already validated at all 2,265,369 input SNP positions and was not repeated.

All 45,303 deposited 50-SNP windows map to 50 genotype rows. The 219 remaining rows are included in the original whole-matrix summaries but not in source-window tables. They contribute no information on LGI and only eight informative patterns on LGIV/LGV combined, negligibly changing chromosome D values. This explains the table-scope difference and is not an estimator error.

Ascertainment common-support filters based on physical span, Ndis callability, and P1/P2 callability retain the regional contrast. This shows the measured covariates alone do not explain it, but cannot correct unknown upstream variant discovery or restore invariant sites. The source panel has variable physical spans and no callable-base denominator.

The strongest contradictory result is sample sensitivity on LGIV: full-group D=+0.0715, while the previously reported matched-pair set gives about −0.740. Further decomposition shows the paired signal is dominated by CJ15(mat-A)/CJ16(mat-a), whereas CJ03/CJ04 has only four informative sites; balanced subsets and leave-one-out estimates vary widely. This prevents any broad claim of sample-insensitive genomic behavior. LGI itself is stable to leave-one-out and balanced subsampling, particularly inside the focal interval.

The original circular-shift p=0.010526 is a reproducible output of 95 circular rotations including the observed alignment (94 nonidentity rotations). Circular chromosome stationarity is not a defensible biological null for this linear, unevenly ascertained LGI track. We retain that p-value only as a description of that algorithm and withdraw strong significance language. The physical-block contrast is the more transparent uncertainty summary, with the few outside blocks stated explicitly.

## Manuscript wording boundary

Report this as a rooted, frequency-weighted nucleotide D contrast on the deposited variable-site panel, elevated across much of LGI and highest in the published interval. State the positive outside baseline, high immediate flanks, ascertainment and outgroup limitations, sample-specific LGIV sensitivity, and the prior introgression evidence. Describe the local-tree result as TWISST q3 association, not exact D_topology. Do not claim a neutral outside baseline, unbiased genome-wide D, exact breakpoints, a uniquely causal MSRC, or proof for/against introgression.
