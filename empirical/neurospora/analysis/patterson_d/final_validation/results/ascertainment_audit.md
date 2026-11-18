# Variable-site ascertainment audit

## What the deposited panel is

The input is a coordinate-indexed matrix of 2,265,369 filtered variable nucleotide rows, not a random sample of callable genomic positions. There are no invariant sites and no per-base callable denominator. The primary D therefore describes allele-sharing among the deposited, polarizable variable sites under the specified genotype and outgroup filters; it is not an unbiased genome-wide nucleotide D.

Martin & Van Belleghem (2017) used L10 samples chosen because this lineage had a previously reported strong *N. hispaniola*-associated signal. Their local genealogies use nonoverlapping windows of 50 SNPs and require every sample in a tree to be genotyped at ≥40 of the 50 SNPs. The deposited coordinate table contains 45,303 windows, each covering exactly 50 genotype-matrix rows; physical span varies widely. The genotype archive is labeled `min20MAC1`, but neither the source archive documentation available locally nor a deposited variant-generation script specifies the exact upstream SNP discovery, allele-count, per-sample inclusion, or missingness rules. We do not infer those details from the filename.

The 50-SNP unit is the genealogy/window sampling unit, not a guarantee that SNPs are independent, evenly spaced, or representative. Physical spans range from 78–121,395 bp on LGI and vary across LGs. LGI median span is 768 bp inside and 649 bp outside (50 SNPs correspond to approximately 65,104 vs 77,042 rows/Mb at those medians). Window spans overlap broadly across focal and complement regions. Per-window median Ndis callability is 0.98 inside and 0.94 outside; pooled row callability is lower (90.8% inside and 85.7% outside) because of the skewed distribution. Median P1/P2 callability is 1.0 in both regions; P3 callability medians are 0.96 inside and 0.94 outside. These are call rates at already-ascertained SNP rows only.

## D sensitivity to observed panel properties

`ascertainment_diagnostics.tsv` provides physical span, SNP density, group callability, outgroup callability, polarizable/informative counts, ABBA/BABA, and D for every source window. `ascertainment_common_support.tsv` reports pooled-D comparisons after restricting windows to overlapping ranges of physical span, Ndis callability, P1/P2 callability, and informative-site density, plus a D-blind nearest-neighbor match on span and raw callability. These filters retain the regional difference: common-support estimates remain near +0.871 inside and +0.446 outside; the covariate-matched subset remains strongly positive in both regions (inside about +0.84, outside about +0.45).

These common-support analyses show that the observed contrast is not explained solely by the marginal differences in source-window span or observed outgroup/ingroup callability. They do not correct the upstream variable-site ascertainment, recover missing invariant sites, or prove that the selection process is ignorable. Informative-site density is itself related to the D-generating patterns, so filters on that diagnostic are descriptive rather than a causal adjustment. No site-level reweighting is used to claim recovery of unobserved variants.

## Rows outside source tree windows

Of 2,265,369 matrix rows, 2,265,150 are assigned to the 45,303 recorded 50-SNP windows. The remaining 219 are included in the original D totals and physical bins, but have no source tree window. They contribute zero informative patterns on LGI; across LGII–VII only eight informative patterns (six LGIV, two LGV) occur. Removing them changes the displayed LGIV/LGV D by less than 0.0001 and does not change the regional inference. Details are in `unassigned_219_audit.tsv`.

## Conditional conclusions

Conditional on this deposited variant panel and its callability filters, the positive LGI site-pattern asymmetry, the higher focal-region D than distant LGI, and the elevated immediate flanks are all observed and robust to measured span/callability common-support restrictions. The panel cannot support a claim that these estimates are genome-wide unbiased, that the entire outside region is a neutral baseline, or that SNP ascertainment has no effect on magnitudes or spatial patterns. No independent unascertained nucleotide panel was present locally; no raw reads or large downloads were initiated.

## Source

Martin & Van Belleghem 2017, [published methods and Neurospora Figure 3](https://doi.org/10.1534/genetics.116.194720); data are described by the [Dryad record](https://doi.org/10.5061/dryad.4jm83). The upstream whole-genome comparison and evidence for introgression are in [Corcoran et al. 2016](https://doi.org/10.1101/gr.197244.115).
