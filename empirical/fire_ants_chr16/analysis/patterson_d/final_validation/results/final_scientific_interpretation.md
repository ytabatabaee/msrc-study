# Final scientific interpretation

## Findings

The frozen Q1 and Q2 estimators are numerically well specified and the
original unit tests cover key polarization, call, and pooling cases. Q1 has
supergene D=+0.8249 versus chr16A+B D=−0.0035 (ΔD=+0.8284); Q2 has +0.8443
versus −0.0192 (ΔD=+0.8635). Chr1 controls are −0.0075 and −0.0452. The
transversion-only and existing geographic/outgroup sensitivities preserve the
direction. These values summarize the sampled BUSCO gene-body SNPs, not all
chromosome bases.

The predeclared moderate filter grid preserves the supergene contrast. For
DP≤12 and call fractions ≤0.70, Q1 D ranges +0.818 to +0.826 and Q2 +0.837 to
+0.847, with thousands of informative sites per estimate. Contrasts versus
chr16A+B remain strongly positive. The strict DP≥15/80% sensitivity is
inconclusive: Q1 has zero informative supergene sites, Q2 has two, and the
comparison region has none. It is not a passed robustness test. Filter-specific
250-kb physical-block intervals for DP≤12 and call fractions ≤0.70 remain
positive (Q1 lower limits +0.771 to +0.786; Q2 +0.785 to +0.815). Strict-filter
intervals are NA because there are too few informative blocks.

Excluding windows with variant-site call rate below 0.80 leaves supergene D at
+0.8266 (Q1) and +0.8439 (Q2). Deleting any single supergene window leaves D
+0.812 to +0.831 (Q1) or +0.837 to +0.851 (Q2). These checks reduce concern
about one low-callability window or one author window driving the regional
contrast, but do not replace individual- or gene-level influence analyses.

The strongest spatial contrast is supported by physical-block intervals across
250 kb, 500 kb, and 1 Mb resampling. The circular-shift p=0.0208 does not
preserve physical spacing or gene density and should remain descriptive. The
nearby window sampling is sparse and linked; chr16B has only two windows.

Exact D_topology=q3−q2 correspondence is positive: across chr16, r=0.838 (Q1)
and 0.762 (Q2); within the supergene, r=0.544 and 0.466, with 1-Mb block
intervals 0.158–0.752 and 0.249–0.668. The within-stratum-centered whole-axis
correlations are 0.441 and 0.356. This supports spatial correspondence but
not independent replication or a unique mechanism.

In(16)1 has 47 midpoint-assigned nucleotide windows and carries most of the
window count; In(16)2 has four, and In(16)3 has zero. Exact SNP-base inversion
overlap is unavailable. GT:DP cannot rule out reference-mapping bias, copy
number effects, or allele-specific dropout. Individual leave-one-out and
per-gene influence are not established from retained outputs.

## Hypotheses and decision

H1 (D differs from zero) and H2 (supergene exceeds chr16A+B) have strong
descriptive support under the primary estimator. H3 (site/genealogy spatial
correspondence) is supported descriptively. H4 (robustness to all genotype,
sample, and gene ascertainment choices) remains incomplete because strict
filters lose information, individual/gene influence is not tested, and
read-level mapping controls are unavailable. H5 (unique proof of introgression
or MSRC) is not established.

**CONDITIONAL GO** for a narrow report of a large, locus/window-specific
derived-allele-sharing asymmetry concordant with sampled quartet genealogies.
Do not claim genome-wide callable-base D, exclusion of mapping/reference bias,
an In(16)3 nucleotide effect, independent replication, or a uniquely identified
mechanism. Competing explanations include recurrent introgression, persistent
ancestral arrangements, ancestral structure, selection on alternative
backgrounds, gene conversion, geography, regions of suppressed recombination, copy-number
variation, reference bias, and genotype ascertainment. The causal order between
structural variation and reduced recombination remains open.
