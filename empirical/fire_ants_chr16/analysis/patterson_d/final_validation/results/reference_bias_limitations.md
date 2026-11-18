# Reference and structural bias limitations

| Concern | Status | Available evidence / limit |
|---|---|---|
| SB/Sb variant-site callability differences | Tested and supported as measurable | Population-group DP≥10 call rates are tabulated in `structural_bias_diagnostics.tsv`; variant-only genotypes cannot measure invariant-site dropout. |
| GT:DP depth differences | Tested and supported as measurable at retained variants | DP thresholds are evaluated in the filter grid. No GQ, allele depths, mapping quality, or read-level evidence accompanies the release. |
| Fixed REF versus fixed ALT polarization | Tested and supported as correctly handled | VCF REF is not assumed ancestral; the fixed outgroup determines orientation. A REF/ALT class association alone does not prove mapping bias. |
| Missingness conditional on genotype/arrangement | Partly tested and supported at segregating VCF sites | Genotype classes absent from the VCF or filtered upstream cannot be recovered. |
| Low-callability BUSCO/window influence | Tested and supported at window level | Window call rate and D are available. Gene-level reference behavior cannot be disentangled from gene ascertainment without per-read data. |
| Allele-sharing SNP depth/callability | Partly tested and supported for retained biallelic records | The source VCF is already filtered and lacks read evidence; no causal interpretation of DP association is possible. |
| SNP-specific ABBA/BABA contribution versus DP/callability | Not tested at individual-site level | The additive scan retains threshold and window aggregates, not per-site allele-sharing contributions. |
| Alignment/reference mapping bias | Not testable with this release | No BAM/CRAM/FASTQ, allele depths, mapping qualities, or genotype likelihoods. |
| Copy-number and multimapping effects | Not testable directly | No depth-normalized copy-number or uniquely mapped-read analysis is supplied. |

Consequently, mapping bias is neither ruled out nor quantifiable here. The most
defensible statement is that GT:DP and variant-site callability can be audited,
while allele-specific read mapping and structural copy-number effects remain
unresolved.
