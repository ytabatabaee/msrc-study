# Ascertainment, callability, and structural-bias audit

## Scope of the analyzed sequence

The local-tree grid contains 213 authors' windows (117 chr1, 42 chr16A, 2 chr16B, 52 chr16_supergene); each window is a four-BUSCO-gene set. The nucleotide calculations include only biallelic SNP records whose reference-contig positions fall inside those exact gene bodies. Across all target windows, 311,890 of 2,650,251 released VCF records map to the 852 genes. The primary quartets retain 91,927 biallelic SNP records in the 52 supergene windows and 49,954 + 1,969 in chr16A + chr16B. These are ascertained variable sites, not genome-wide invariant-plus-variable alignments. The target-gene body span is 1,751,361 bp in the supergene windows and 1,019,759 bp in chr16A+B; this is sequence-span coverage, **not** an estimate of callable bases.

Physical-bin tracks at 50 kb, 100 kb, 250 kb, 500 kb, and 1 Mb pool window-assigned site contributions using the frozen author coordinates and midpoints. Each row includes its contributing-window count, variant/SNP counts, and informative counts. They are sparse BUSCO-gene tracks. The gene identity bridge does not provide a continuous scaffold-to-chromosome coordinate transformation; no SNP is assigned an invented chromosome-16 base coordinate.

## Genotype and polarization filters

- Minimum genotype depth is DP 10, consistent with the released pipeline's hard-call mask. The final VCF has GT:DP only; GQ, AD/RO and per-read mapping quality are unavailable for this audit.
- A haploid individual's `0`, `1`, `0/0`, `1/1`, `0|0`, or `1|1` call contributes one allele. Discordant pseudo-diploid heterozygotes are excluded as ambiguous; they are already absent or masked in the supplied filtered VCF, so the observed ambiguous-call count is zero.
- Worker pools and source-pipeline random-allele samples were excluded from frozen membership. They are not treated as independent haploid individuals.
- Each ingroup population needs at least 50% of its frozen individuals called at a site. Outgroup calls must be unanimous among called individuals, with at least four calls (or all available individuals when the frozen outgroup subset has fewer than four; this applies only to the two-male French Guiana sensitivity).
- Fixed REF outgroup makes ALT derived; fixed ALT outgroup makes REF derived. Missing, polymorphic, or discordant outgroups are excluded.
- Indels and multiallelic records are excluded. All ABBA/BABA contributions are allele-frequency products, pooled before D is calculated.

## Observed variant-site callability

At biallelic variant positions, Q1 call rates are 0.830 in the supergene and 0.840 on chr1; Q2 rates are 0.870 and 0.881, respectively. Primary quartet call rates vary little between supergene and flanks; the Q1/Q2 call-rate differences mainly reflect their different sample membership and denominators. Per-region missing calls, group-callability exclusions, outgroup failures, polarized sites and informative counts are in `variant_diagnostics.tsv`. The VCF has no invariant-site genotype matrix, so these rates do not measure base-level mappability or callable sequence.

The deposited calls were mapped to a single Si_gnGA reference and then quality-filtered upstream. The Sb haplotype is known to differ structurally and in repeat/copy-number content. This local release cannot test allele-specific mapping bias, copy-number effects, or mapping-quality differences because it lacks reads, per-read alignments, and genotype likelihoods. Upstream filters did remove problematic high-coverage/repetitive regions, but they do not establish equal mapping behavior between SB and Sb. Thus the strong regional D is robust to the tested nucleotide substitution class and sample sets, but structural/reference bias remains a plausible alternative contributor.

## Sensitivity status

Transversions-only D, Argentina-only vs all-geography populations, S. invicta-only focal samples, and Brazil vs French Guiana *S. saevissima* outgroup subsets all retain the positive supergene D and near-zero chr16A+B background. A stricter sensitivity (DP>=15 and >=80% calls in each ingroup population) sharply reduces usable data: Q1 has no informative supergene sites; Q2 has only two informative supergene sites and none in chr16A+B. The Q2 strict-filter D is based on two sites and is not an adequately powered confirmation. This sensitivity does not falsify the primary result; it shows that the deposited BUSCO-window VCF cannot support a high-stringency estimate with useful precision.

A different external outgroup species, high-GQ sensitivity, genotype-likelihood estimator, read/mappability/CNV audit, leave-one-individual-out scan, and balanced subsampling design were not available or frozen: no clean alternative species outgroup was frozen, GQ/AD and read-level files are unavailable in this release, and balanced/LOO designs were not frozen in Stage 1. No new sample design was selected after viewing D.

The two chr16B windows are reported separately as a two-window flank summary; they do not provide an independent evolutionary replicate. Treat small-region D and spatial results as descriptive.
