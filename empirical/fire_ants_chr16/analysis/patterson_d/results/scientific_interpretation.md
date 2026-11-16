# Scientific interpretation: fire-ant supergene Patterson's D

## Main estimates

The primary quartets were frozen from metadata and the independent chromosome 1–15 baseline before nucleotide analysis. Q1 is (invicta/macdonaghi SB, invicta/macdonaghi Sb, richteri Sb, *S. saevissima*); Q2 is (richteri SB, richteri Sb, invicta/macdonaghi Sb, *S. saevissima*). All quartets use the same derived-allele orientation and one biological allele per haploid individual.

| Quartet | Region | D | ABBA | BABA | Informative SNPs |
|---|---|---:|---:|---:|---:|
| Q1 primary Argentina | chr16 supergene | +0.8249 | 1,395.49 | 133.91 | 3,696 |
| Q1 primary Argentina | chr16A + chr16B | −0.0035 | 106.77 | 107.52 | 1,584 |
| Q1 primary Argentina | chr1 control | −0.0075 | 285.24 | 289.52 | 4,785 |
| Q2 primary Argentina | chr16 supergene | +0.8443 | 1,431.84 | 120.88 | 3,540 |
| Q2 primary Argentina | chr16A + chr16B | −0.0192 | 60.83 | 63.21 | 1,312 |
| Q2 primary Argentina | chr1 control | −0.0452 | 171.92 | 188.19 | 3,921 |

The pooled supergene-minus-chr16A+B contrasts are +0.8284 (Q1) and +0.8635 (Q2); versus chr1 they are +0.8323 and +0.8895. Estimates use 52 supergene windows, 42 chr16A windows, 2 chr16B windows, and 117 chr1 windows. The two chr16B windows are too few to support a stand-alone regional inference. Dropping the first and last supergene windows or adding the nearest one window on either side changes the contrasts only slightly (Q1 +0.8286 trimmed / +0.8260 expanded; Q2 +0.8632 / +0.8624 versus chr16A+B).

Across 250-kb, 500-kb, and 1-Mb physical block resampling, the 95% block-bootstrap intervals for the inside-minus-chr16-outside D remain positive: Q1 +0.784 to +0.874, +0.779 to +0.878, and +0.782 to +0.874; Q2 +0.813 to +0.919, +0.804 to +0.927, and +0.805 to +0.926. Corresponding delete-one-block jackknife intervals also remain positive. These are genomic-block uncertainties, not independent-SNP standard errors. Circularly shifting the 52-window supergene mask over the 96 chr16 author windows gives an exact two-sided p=0.0208 for each primary contrast; this is a coarse, linked-window spatial null with only 96 rotations and is descriptive.

## Robustness and topology matching

The positive supergene signal is stable across the frozen sensitivities: all-geography, S. invicta-only focal samples, Brazil-only *S. saevissima* outgroup, and French-Guiana-only outgroup each give Q1 supergene D about +0.802 to +0.832, with chr16A+B close to zero and Delta D about +0.785 to +0.834. Q2 all-geography and S. invicta-only supergene estimates are +0.844 and +0.850. Transversions-only supergene D is +0.829 (Q1) and +0.852 (Q2), with transversion-only inside/outside contrasts +0.819 and +0.896. The primary orientation is therefore not driven by transitions or by the tested geographic subset. There is no clean alternate outgroup species in the Stolle release.

The strict DP>=15 / 80%-ingroup-call sensitivity is severely underpowered in these windows: Q1 has zero informative supergene sites, and Q2 has two informative supergene sites in the supergene with none in the chr16A+B comparison. The D estimates from that filter cannot be treated as an independent confirmation or contradiction; the strict sites do not provide enough information for a useful regional comparison. This is a material limit on claims of robustness to callability stringency.

All four sample groups, including the external outgroup, occur in the original 267-tip Newick trees. The mapping validates every frozen sample ID against local tips. For each tree, q1/q2/q3 were calculated exactly over the Cartesian product of all available frozen sample combinations using the unweighted four-point metric; all 1,704 quartet-window cases have all groups and are resolved. Thus `D_topology=q3-q2` is available for the matching frozen quartet groups, but remains a topology frequency contrast rather than nucleotide D.

The inside-window mean induced D_topology is +0.828 (Q1) and +0.746 (Q2), versus −0.006 and approximately 0.000 on chr16A and −0.007 and −0.002 on chr1. Nucleotide D correlates descriptively with exact D_topology across chr16 (Pearson r=+0.838 and +0.762) and within the supergene (r=+0.544 and +0.466), but has near-zero correlation outside (r=−0.072 and +0.142). One-Mb block-bootstrap intervals for the supergene correlations are +0.158 to +0.752 and +0.249 to +0.668. Within-supergene D also correlates with the original seven-group q_H−q_S track (r=+0.588 and +0.466), q_H (r=+0.550 and +0.438), and inversely with q_S (r=−0.605 and −0.477); q_3 correlations are small (+0.055 and +0.033). Across all chr16, raw D_site vs D_topology correlations are higher because both measures share a sharp inside/outside transition; after centering each track within inside/outside strata, residual Pearson correlations remain +0.441 and +0.356. Chr1 is a control with no similarly strong spatial relation. These are descriptive correlations across linked windows, not independent-window significance tests.

## What the result supports

Both reciprocal population-history/arrangement quartets show a large, localized excess of ABBA over BABA in the sampled four-gene supergene windows, while the chr16A+B and chr1 comparisons are near balance. The signal agrees spatially with the already observed haplotype-associated TWISST shift and with the exact induced quartet topology contrast. This adds nucleotide site-pattern evidence to the local-genealogy result.

The estimate does not uniquely establish introgression or the MSRC mechanism. The analyzed data are variant-only, restricted to four BUSCO genes per window, and mapped to one reference in a structurally divergent region. D cannot distinguish introgression from ancestral structure, selection on ancestral polymorphism, gene conversion, or systematic allele sharing induced by mapping/reference bias. It cannot establish transfer direction or reconcile the competing published histories by itself. PNAS and Stolle used different samples, assemblies, site sets, inversion boundaries, and hypotheses. The justified manuscript claim is a localized, robust asymmetry in derived-allele sharing that is concordant with the supergene genealogy shift, with mechanism stated as consistent with but not proof of introgression.

**Assessment: CONDITIONAL GO** for inclusion as carefully qualified, locus/window-specific nucleotide evidence. Do not describe these results as genome-wide callable-base estimates, a standalone causal test, or definitive proof of introgression/MSRC.
