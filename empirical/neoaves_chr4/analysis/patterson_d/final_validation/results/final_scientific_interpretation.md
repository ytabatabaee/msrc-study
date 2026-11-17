# Final scientific interpretation

## Assessment: CONDITIONAL GO

A carefully qualified report of the chromosome-4 nucleotide site-pattern results is supportable. The strongest defensible result is a replicated regional N61 contrast in the frozen quartet, stable to the *Meleagris* alternative outgroup, transversion-only filtering, omission of PNAS_2 or the largest interval, and ±250-kb to ±1-Mb boundary perturbations on the original source coordinate. This is not evidence that one evolutionary mechanism uniquely caused the pattern.

## Distinct claims

### A. Is D different from zero?

At the 1-Mb block-jackknife scale, the primary all-sites D estimates outside and inside the union of the three published intervals are both different from zero under the descriptive 95% intervals for all three frozen designs. This is not a valid conventional introgression test for every quartet: Columbea's rooted species-history relationship is disputed, and all three show nonzero outside baselines. In particular, a nonzero D outside the target does not localize an effect to a suppressed-recombination region.

### B. Does D change inside versus outside?

- **Columbea:** inside 0.2164 (SE 0.0107), outside 0.2168 (0.0052), ΔD −0.0004 (0.0122; 95% block interval −0.0243 to 0.0235). No pooled regional shift. PNAS_2 alone is higher but has only 776 informative sites and 43 loci; this does not overturn the pooled result.
- **N61:** inside −0.2582 (0.0098), outside −0.0806 (0.0040), ΔD −0.1776 (0.0106; interval −0.1984 to −0.1568). The sign and regional contrast persist with *Meleagris* (ΔD −0.1666), transitions only (−0.1424), transversions only (−0.3181), excluding PNAS_2 (−0.1717), excluding PNAS_3 (−0.1565), and boundary expansion/erosion sensitivities (−0.154 to −0.178). The interval-specific contrasts are negative in PNAS_1, PNAS_2 and PNAS_3; PNAS_2 remains less precise.
- **N62:** inside 0.0394 (0.0111), outside 0.0385 (0.0041), ΔD 0.0009 (0.0120; interval −0.0226 to 0.0243). No all-sites pooled regional shift. The transition-only ΔD is −0.0061 (0.0121), and transversion-only ΔD is +0.0385 (0.0225); both 95% block intervals cross zero. Alternative outgroup ΔD is +0.0162 (0.0098), also compatible with zero.

All SEs and intervals use genomic blocks; they are not site-independent standard errors. `D_outside=0` and `D_inside=D_outside` are separate hypotheses and are tabulated separately.

### C. Does nucleotide D correspond to local genealogy?

Exact four-tip matching was possible. Across exact-tree loci the genome-wide exact q3−q2 topology contrasts were +0.015 for Columbea, −0.230 for N61 and −0.016 for N62. Matched 1-Mb D_site vs D_topology correlations were, respectively: Columbea r=0.088 whole chr4 (0.445 inside, 0.052 outside; 0.078 after inside/outside demeaning); N61 r=0.914 (0.649 inside, 0.755 outside; 0.732 after demeaning); N62 r=0.355 (0.405 inside, 0.399 outside; 0.375 after demeaning). These are descriptive spatial correlations. Two-Mb block-bootstrap intervals are in `spatial_validation.tsv`.

The N61 nucleotide/site and exact topology contrasts agree in sign and co-vary strongly. Columbea has a large nucleotide D with a near-zero q3−q2 frequency contrast, and N62's topology contrast is also near zero. This separation is consistent with topology frequencies being only one component of site-pattern asymmetry; branch lengths, substitution processes, polarization, and gene-tree estimation also matter. The induced unrooted gene-tree internal branch-length summaries are machine-readable in `quartet_branch_lengths.tsv`; they are ML substitution distances, not coalescent times.

### D. Is a specific evolutionary mechanism established?

No. Nonzero D does not prove introgression. An inside/outside change does not uniquely demonstrate MSRC. The analyses do not distinguish introgression from ancestral structure, asymmetric ILS, finite-sites homoplasy, taxon/history effects, or other model violations. Structural boundaries remain proposed genomic annotations; they are not demonstrated causes of suppression in these tests.

## Why the backgrounds are nonzero

Columbea’s positive D persists with the alternative Galliform outgroup and among transversions (transition-only pooled ΔD is near zero), so Gallus-specific polarization and transitions alone do not explain its positive background. Yet the 48-taxon and 363-taxon references do not agree on the exact frozen quartet, and the exact local-tree q3−q2 contrast is near zero. The available data cannot assign one cause. For N61, the outside D is negative and the exact gene-tree q2 frequency exceeds q3; genome-wide q3−q2 and D_site are both negative. This offers a supported description of the background asymmetry, not a causal explanation. Its additional inside shift is reproducible, but not uniquely attributable to gene flow or suppressed recombination.

## Limits on manuscript wording

Describe N61 as a robust inside/outside nucleotide D contrast in this fixed taxon quartet over the published source-coordinate intervals. State the nonzero outside baseline, quartet uncertainty, and linked genomic sampling. Present Columbea and N62 pooled regional results as near-null. Do not translate the results into proof of introgression, timing/direction, an inversion, or the origin of a suppressed-recombination region.
