# Genotype filter sensitivity interpretation

The primary analysis remains DP≥10 and at least 50% call fraction per ingroup.
The additive grid crosses DP thresholds 8/10/12/15 with call fractions
0.50/0.60/0.70/0.80 for both frozen primary quartets. It is a predeclared
diagnostic grid, not a search for the largest D. See `filter_attrition.tsv` for
stage attrition and `genotype_filter_sensitivity.tsv` for pooled regional
estimates and informative-site counts.

The grid's `records` field counts mapped biallelic SNP records eligible for
genotype parsing. The complete source scan separately reports 2,650,251 VCF
records genome-wide and 311,890 records overlapping the 852 gene bodies; these
source totals are documented in the parent feasibility audit. Variant positions
overlapping multiple gene intervals are not synonymous with author-window
duplicates.

Across non-strict combinations, supergene D ranges +0.801 to +0.840 for Q1
(55–3,744 informative sites) and +0.833 to +0.847 for Q2 (137–3,588 sites).
Restricting to DP≤12 and call fractions ≤0.70 retains thousands of informative
supergene sites per design and gives Q1 D +0.818 to +0.826 and Q2 +0.837 to
+0.847. The supergene-minus-chr16A+B contrasts remain positive across these
moderate combinations (roughly +0.82 to +0.89). Their 250-kb block-bootstrap
95% percentile intervals are wholly positive: Q1 lower bounds 0.771–0.786
(upper bounds 0.866–0.876), and Q2 lower bounds 0.785–0.815 (upper bounds
0.913–0.926). Each has at least 34
informative inside blocks and 22 outside blocks. At DP15/80%, Q1 has zero
informative supergene sites, Q2 has two, and the comparison region has none;
among 91,927 mapped biallelic supergene records, only 61 Q1 and 222 Q2 pass
the strict ingroup call-fraction rule. Outgroup fixation then leaves 47 and 194,
respectively, before ABBA/BABA informativeness is assessed. At the primary DP10
threshold, variant-site call rates are lowest in Q1 P2 and Q2 P3 (both about
0.795), compared with 0.969 in the outgroup. Thus strict attrition is driven
primarily by the combined high-depth/high-callability requirement in the
ingroup, with a secondary loss at outgroup polarization; it is not evidence
that the primary signal is absent. The strict estimates are NA/underpowered for
the regional comparison, not D=0 and not confirmation. At the most
stringent nonzero combinations counts and block counts drop sharply, so the
sign alone is not evidence of robustness. Per-combination counts and intervals
are in `genotype_filter_sensitivity.tsv`; per-window contributions are in
`filter_window_sensitivity.tsv`.

The VCF already underwent upstream filtering, so changing these thresholds
cannot recover genotypes discarded before release.
