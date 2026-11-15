# Methodological audit before matched inference

## Sample identities and group sizes

The exact frozen member file is treated as authoritative. The archived source mapping and VCF header both verify the same 55 exact IDs; each is present in the tree-tip union and the indexed VCF header (67 total VCF samples). Per design, group sizes are:

| Design | P1 | P2 | P3 | Outgroup |
|---|---:|---:|---:|---:|
| D1 Mmd t vs Mmd standard | 4 | 12 | 13 | 8 |
| D2 Mmm t vs Mmm standard | 8 | 13 | 12 | 8 |
| D3 Mmc t vs Mmc standard | 3 | 7 | 13 | 8 |

GT coding is VCF diploid `0=REF`, `1=ALT`; `0/1` contributes ALT frequency 0.5, and incomplete diploid calls are missing. Population frequencies average called individuals. The outgroup identities are the eight frozen `Ms_SPRE...` samples.

## Coordinates, trees and filters

The source README states that local-tree window labels are chr17 coordinates minus 5,000,000 bp. Tree windows were restored to mm10/GRCm38 chr17 coordinates and matched to VCF windows by exact 5-kb start. The source archive has 4,046 usable trees spanning chr17:5,000,000–37,055,000; no gene-tree inference is made beyond that coverage. VCF POS is 1-based and is assigned to 0-based half-open 5-kb windows using `(POS−1)//5000`.

The nucleotide input is Harr et al.'s indexed AllMouse PASS VCF. The scan retains only biallelic one-base A/C/G/T SNPs. The original implementation required all outgroup calls fixed REF; the extension fixes that polarization omission by including fixed ALT as the opposite ancestral-state case while continuing to exclude polymorphic or missing outgroups. No minimum derived allele count is imposed. The published four inversion intervals and existing collinear/flank masks are reused unchanged.

## Material comparability limitation

The sample identities match, but the carrier sequence tips in the Kelemen & Vicoso trees are reconstructed pseudo-t-haplotype sequences made from t-specific SNP profiles, whereas the nucleotide scan uses unphased diploid carrier genotypes. Thus exact population membership is matched, but the chromosome copies represented by the two measurements are not identical. Results are interpreted as population/individual-level correspondence, not a chromosome-phased t-haplotype comparison.

Arrangement assignments come from the frozen source-backed sample metadata and `_tHaplSubset` designation. No arrangement label is inferred from local tree topology. The M. spretus cohort is the only outgroup present in the frozen source design; no alternative outgroup was selected.
