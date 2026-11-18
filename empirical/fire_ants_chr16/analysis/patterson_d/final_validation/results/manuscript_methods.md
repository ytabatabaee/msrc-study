# Manuscript methods draft

We estimated frequency-weighted Patterson's D from the Stolle et al. variant-only
GT:DP VCF using the frozen Q1 and Q2 male-individual quartets. VCF biallelic SNP
records were assigned by scaffold and position to the exact four BUSCO gene
bodies composing each author-defined window. Haploid calls and homozygous
pseudo-diploid calls were represented as one allele; heterozygous, missing, or
below-depth calls were excluded. The primary filter required DP≥10 and at least
50% called individuals in each ingroup population. Outgroup alleles had to be
unanimous among at least four called *S. saevissima* individuals. Fixed-REF and
fixed-ALT outgroups separately determined the derived allele. Site ABBA/BABA
frequency products were pooled before computing D.

Q1 contrasted invicta/macdonaghi SB and Sb against richteri Sb, rooted with
*S. saevissima*. Q2 contrasted richteri SB and Sb against invicta/macdonaghi Sb,
using the same outgroup. Membership was fixed from sample metadata and the
independent chr1–15 species history before nucleotide results. Regional
uncertainty used physical genomic-block resampling at 250 kb, 500 kb, and 1 Mb.
Exact induced quartet topology frequencies used matching frozen sample groups
and D_topology=q3−q2; these are treated as distinct from nucleotide D.

Coordinates shown for chr16 are author-window coordinates on the Stolle axis.
SNP coordinates remain on VCF scaffolds; no continuous SNP lift is implied.
Inversion window classes use Stolle's source midpoint rule and existing mapped
outer envelopes. The site set is restricted to selected BUSCO genes and
variant-only records, so results are conditional on this ascertainment.
