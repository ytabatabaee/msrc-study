# SNP ascertainment and source filtering audit

## Documented processing

The deposited files are derived from the `mmatschiner/supergenes` cod phylogenomics workflow and Zenodo record 4560275. The upstream `apply_filter.slurm` workflow applies site quality filters (FS < 20, QD > 2, MQ > 20, depth within the 10th–90th percentile, InbreedingCoeff > -0.5, and read-position/mapping-rank filters), masks indels/repeats and other masked positions, then sets individual genotypes to missing when DP < 4 or GQ < 20. It removes unobserved alternate alleles and retains biallelic SNPs with both alleles observed.

`split_vcf.sh` derives the SNAPP and Dsuite site sets from filtered/masked VCFs using biallelic SNP filters and a minor allele count of at least one across the full input cohort. The SNAPP VCF additionally excludes any site with a missing genotype anywhere in its full sample roster; the Dsuite `.full` VCF permits missing genotypes. Thus the site spectrum is globally ascertained, and the two sources differ in complete-call filtering. The project does not claim these SNPs are an unbiased genome-wide sample.

The SNAPP preparation script declares a 250-kb window, a random maximum of 1,000 SNPs per window, and 50-bp minimum distance between selected SNPs. The nominal 500-site minimum is set in the shell wrapper but is not passed by the published SLURM command. These XML-preparation filters are distinct from the complete SNAPP VCF files used here for D; no random 1,000-SNP thinning or 50-bp thinning was applied by this validation.

The Relate input was Beagle-imputed from filtered cod genotypes with ancestor-heterozygous sites removed and a repeat/indel/ancestral-heterozygosity mask. It contains cod samples only, no external outgroup, and no ancestral-state field in the VCF. The deposited colinear Relate VCF is restricted to LG03–LG06, LG08–LG11, and LG13–LG23, so it supplies no sites in the four inversion linkage groups under study.

## Validation estimator

For each frozen population design, the scanner parses diploid genotypes directly, averages allele frequencies across the two preselected individuals per ingroup population, and requires complete genotypes in all six selected ingroup individuals. It includes a site only when the chosen outgroup is homozygous reference or homozygous alternate. Fixed ALT outgroups are correctly repolarized; heterozygous and missing outgroups are excluded. Multi-allelic and non-SNP records are excluded. REF is never treated as ancestral by default. Dsuite results are shown only where the deposited samples permit the frozen design and should be interpreted as source/filter sensitivity, not a new independent biological replication.

The positional overlap audit compares VCF records in the exact 250-kb published tree windows, excluding boundary-overlapping windows in the regional contrasts. D is computed from pooled ABBA/BABA contributions, never as a mean of window-level D values.
