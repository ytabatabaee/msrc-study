# Fire-ant chromosome 16 Patterson's D: Stage 1 data feasibility

**Audit outcome: feasible for a second-stage nucleotide analysis, conditional on explicit filtering and using the frozen designs. No D values were calculated in Stage 1.**

## Data actually available

The Stolle et al. release provides a real, sample-level VCF at the Wurm Lab data page. The nested download URL is `https://github.com/wurmlab/wurmlab.github.io/raw/master/data/supergene_introgression/gt.vcf.gz/gt.vcf.gz`; the shallower URL is only a Git LFS pointer. The local file is 1,503,782,377 bytes; SHA256 `e782b2e523357ee09ffe81ce0cc2c6b25f733159ba965c638cc03c6785306775`; it matches the Git LFS object ID and passes `gzip -t`. The VCF is VCFv4.2, called with FreeBayes, has 69,592 contigs and 368 individual sample columns, and contains `GT:DP` calls. Its reference header is `gng20170922wFex.fa`.

The local Supplementary Data 1 has 386 sample rows. The released VCF matches 368 of them exactly; the 18 nonmatching rows agree with the paper's stated removal of nine high-missingness and nine suspected diploid males. Every released VCF sample has a metadata match. All 267 frozen local-tree tips match exactly to both VCF and metadata. The source paper reports the 368-sample filtered release and links the genotype file and sequence archive ([Stolle et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8917144/); [Wurm Lab data](https://wurmlab.com/data/supergene_introgression/)).

## What is and is not accessible

| Requirement | Finding |
|---|---|
| Individual VCF genotypes | Available locally; 368 samples; GT:DP |
| Population frequencies | Not supplied as a separate product; can be estimated from individual genotypes and the frozen metadata |
| Nucleotide SNP matrix | Available as VCF variant records; not an alignment of invariant sites |
| Sample identifiers and population metadata | Available; exact matching for every VCF sample |
| Outgroup genotypes | Available for *S. saevissima*; six clean individual males retained in the local tree and VCF |
| Reference and genomic positions | VCF reference and scaffold positions available |
| Random-allele/pool caveats | Documented in original filtering scripts; flagged samples and pools are excluded from primary quartets |
| Tabix/CSI index | Not supplied at the public endpoint (TBI URL returned 404); no local index or `bcftools`/`tabix` executable was present during audit |
| Whole-genome raw reads | Not downloaded |

The full VCF was streamed locally because there is no index for efficient range queries. This scan was completed and checksummed. Do not treat raw REF as ancestral: future D computation must polarize using the frozen outgroup calls and exclude unsuitable outgroup sites under a predeclared rule.

## Genotype quality and ploidy caveats

The published pipeline aligned reads to Si_gnGA and used FreeBayes with ploidy 2, including haploid males. Its downstream genotype filter resolves many apparent heterozygotes using reference-read fractions, masks low-depth and low-quality calls, removes true heterozygotes from ordinary haploid samples, and applies site-level filters. The released VCF contains only GT and DP, so the original read-count and genotype-quality evidence cannot be independently reconstructed from this release. The source pipeline documents three worker pools represented by majority alleles and seven outgroup males given random-allele treatment at heterozygous sites. These are excluded from the primary sample designs. One random-treatment identifier (`AR142-1-bigB-m`) is absent from the current VCF and supplement release; the other six listed source IDs resolve. A source-script intermediate count of 370 differs from the deposited 368; the exact current deposited roster is used here, with no undocumented reconciliation.

## Coordinate compatibility

The 213 frozen windows consist of four BUSCO genes each (852 unique genes). All 852 IDs match rows in Stolle et al. Supplementary Data 2, and their 41 scaffolds occur in the VCF contig header. The VCF scan found 2,650,251 variant records; 311,890 records overlap the exact gene-body intervals of those 852 window genes, including 3,424 positions overlapping multiple gene intervals. This supports exact **window membership** by joining VCF contig/position to the BUSCO interval and gene ID. The upstream external scaffold-to-chromosome linkage map is absent from this repository snapshot, so these results do not claim a continuous scaffold-to-chromosome coordinate lift or assign SNPs to invented base-pair positions on chr16. The 1-kb flanks used by some upstream steps are outside the exact gene-body membership used here.

## Stage 1 decision

**CONDITIONAL GO** to implement nucleotide D on exact gene-body variants assigned to frozen windows. Conditions: preserve the documented haploid/pool/random-allele exclusions, use outgroup-based polarization rather than REF-as-ancestral, report variable-site ascertainment and callable-site limitations, and keep scaffold/window coordinates distinct from chromosome coordinates. The analysis should not claim genome-wide invariant-site density from a variant-only VCF.

Machine-readable inventories: `nucleotide_source_inventory.tsv`, `sample_mapping_audit.tsv`, `window_gene_coordinate_mapping.tsv`, `vcf_record_scan_by_contig.tsv`, and `window_vcf_coverage.tsv`.
