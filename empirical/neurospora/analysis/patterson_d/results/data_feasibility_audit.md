# Nucleotide data feasibility audit — completed follow-up

## Decision

**GO for a conditional rooted nucleotide D calculation**, after the user supplied the checksum-matching Dryad outgroup ZIP. This is a local analysis; no fresh download or raw-read processing was used. The deposited Martin & Van Belleghem archive and user-provided outgroup ZIP are present under ignored `data/neurospora/raw/`.

## Source inventory and access status

- Martin & Van Belleghem (2017), Zenodo record 4946041: the local archive is 15,577,951 bytes, MD5 `571324702709eb1c2732a3fd993edc93`. It contains the L10 haploid nucleotide matrix, coordinate-indexed 50-variable-site windows, trees and topology weights. The genotype matrix has 2,265,369 rows and 26 sample columns; the source is filtered variable SNPs rather than invariant/callable sequence. See `nucleotide_source_inventory.tsv`.
- Corcoran et al. (2016) Dryad DOI 10.5061/dryad.162mh: exact Ndis alignment file 11939, 63,382,228 bytes, MD5 `95a3eb66abc7d0c44cda5cb7953b2a4f`. The ZIP was user supplied. It contains seven chromosome-group FASTAs, including `Ndis@999` and the gap-free `2509@6` reference row. Its row count/lengths and sequence classes are recorded in `outgroup_alignment_inventory.tsv`.
- Exact coordinate check: genotype POS is 1-based; alignment columns are 0-based, and `2509@6` has no gaps. At all 2,265,369 genotype positions, the matrix reference base equals the corresponding alignment base under column `POS-1`; zero mismatches. The competing `POS+1` mapping does not match. Ndis gaps and ambiguous symbols are excluded per site; callability counts are in `outgroup_callability.tsv`.
- Available nucleotide formats: haploid base matrix and multi-sequence nucleotide FASTA alignment. No VCF/BCF, population allele-frequency file, or REF/ALT ancestral annotations are used. Sample IDs, genomic linkage-group coordinates, and source reference sequence are present. Ndis independently supplies polarization. Individual-to-population membership was checked against the mating-type source audit, not inferred from local trees.

## Remaining limits

The four labels CJ07, CJ11, CJ12 and CJ18 have no independently supported mating-type assignments and are excluded. Ten verified mat-a and eight verified mat-A individuals remain. The L10 trees lack Ndis, so exact four-taxon D_topology comparison remains unavailable. The variable-site panel's ascertainment means reported D is conditional on its SNP filtering and is not a genome-wide callable-site estimate. See `scientific_interpretation.md` and `run_manifest.json`.
