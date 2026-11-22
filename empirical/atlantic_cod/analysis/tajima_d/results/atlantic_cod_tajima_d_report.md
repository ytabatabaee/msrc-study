# Atlantic cod Tajima's D — exploratory, inside-only population summaries

**Status: EXPLORATORY.** D was calculated for each available sampling population separately on the existing Relate SNP VCFs; populations were not pooled to increase sample size. The source script selects biallelic records from `concatenated.masked.vcf.gz`, excludes sites with any missing GT across the full source panel, and removes sites fixed across that panel. The VCF is derived from a filtered/imputed GATK callset. The deposited inputs do not establish that every segregating SNP was discovered and retained. Missing invariant sites do not invalidate classical D, but site filters can omit segregating sites, so these values are conditional on the retained panel.

## Samples, ploidy, assignments, and coordinates

The Relate VCFs contain 24 diploid fish columns, phased GT (`0|1` etc.). Phasing is not needed for site-frequency counts; each diploid genotype contributes two alleles. Population membership is read from the independently named sample prefixes (`avc`, `avo`, `bat`, `bor`, `icc`, `ico`, `kie`, `lfc`, `lfo`, `low`, `twc`, `two`). Eleven groups have at least two fish (four haplotypes); Iceland coastal (`icc`) has one fish and is not estimable. The sample-group manifest also records independently established population-level arrangement states and confidence for each linkage group. These are not individual arrangement assignments, so no pure arrangement haplotype contrast is made. LG12 arrangement orientation has demographic-inference rather than direct-outgroup confidence.

Coordinates are the VCF's linkage-group coordinates on the gadMor2 study axis, matched by exact LG and position to the published 250-kb window inventory. The existing frozen boundaries are LG01 9,114,741–26,192,386; LG02 18,489,307–24,048,607; LG07 13,610,591–23,019,113; LG12 638,100–14,327,837 (1-based source coordinates). The 250-kb windows follow the existing tree-window table. GT-only VCFs contain no depth/callability denominator; records retained by source filtering have no missing genotypes, but missingness before filtering is not recoverable from these files.

## Results and decisive data limitation

The four supergene Relate VCFs provided 530,017 retained records in the selected LG intervals. D is estimable in 1,914 population-window rows (11 populations × 174 existing inside windows). For the More-stationary (`avc`, two fish) population, mean D inside was LG01 **−0.143** (66 windows), LG02 **−0.087** (21), LG07 **−0.067** (36), and LG12 **−0.196** (51). The full per-population table is `atlantic_cod_tajima_d_inside_outside.tsv`.

Crucially, the existing `gadus_morhua_outside_supergenes_relate.vcf.gz` contains **zero records on LG01, LG02, LG07, or LG12** under the audit parser. Therefore there are no estimates for the outside windows on those four linkage groups, and no inside-minus-outside effect or uncertainty interval is reported. The source tree table does contain outside windows, but tree presence cannot supply nucleotide allele counts. The Dsuite and SNAPP VCFs use different sample sets/processing and were not mixed into the Relate comparison. Per-base pi/theta and callability are unavailable.

## Interpretation

These within-population inside-window SFS summaries add descriptive information beyond population-tree branch lengths, but the SNP filtering and lack of matched outside SNPs prevent a test of inversion-associated D shifts. They do not reproduce or refute the Neurospora pooled/background pattern, because no within-population individual arrangement assignments support pooled-versus-background D. No neutrality-test p-values are interpreted. The immediate missing input for regional comparison is an outside-2LG genotype VCF with the same 24 individuals, same call filters, and target-linkage-group flanks; unbiased quantitative interpretation additionally requires evidence on the pre-filter masked GATK segregating-site set.

Reproduce from repository root with `python3 empirical/atlantic_cod/analysis/tajima_d/scripts/run_cod.py`. See the window table, source audit, sample-group manifest, and PDF/PNG figure in this directory.
