# Patterson's D design audit

The four candidate tests use the ordered population design `((P1,P2),P3),O`.
The primary designs use A. merus as outgroup; A. melas is retained as a
sensitivity outgroup. Groups are frozen from the existing sample-level Ag3.10
karyotype calls and Stage-1R primary homozygote table. Heterokaryotypes are
not used as ingroup P1/P2/P3 samples. Outgroup arrangement is not used to
polarize sites; only reliable outgroup nucleotide calls determine the
ancestral allele at each site.

The species-tree baseline is supported by Fontaine et al. (2015) and
Thawornwattana et al. (2018): A. gambiae and A. coluzzii are sister taxa and
A. merus is the earliest branching member of the complex. The design therefore
tests population-level asymmetry conditional on that baseline; it does not
assert that arrangement-defined samples are a species-tree quartet.

The main confounding is geographic structure: standard and inverted samples
are not geographically exchangeable within either species in this frozen
sample set. Ancestral arrangement polymorphism, incomplete lineage sorting,
and introgression can also produce asymmetry. Designs A and B are primary
exploratory tests because the gambiae groups have 8 and 12 samples and use
coluzzii as P3. Designs C and D are sensitivity tests; C_i has only 3 samples.

A. christyi is a more external outgroup in the published whole-genome
framework, but no A. christyi sample is present in the authenticated
`fontaine-2015-rebuild` sample set. It is therefore not downloaded or used in
this analysis. The existing sample metadata contain 9 merus and 4 melas
samples; their genotypes, rather than arrangement labels, determine per-site
eligibility. A fixed-outgroup estimator excludes sites where the selected
outgroup is missing or polymorphic.

The complete frozen table is `frozen_quartet_designs.tsv`. A nonzero result is
not evidence by itself for introgression or MSRC causation.
