# N. discreta alignment and coordinate validation

## Inventory and taxa

The user-provided ZIP matches Dryad's recorded size (63,382,228 bytes) and MD5 (`95a3eb66abc7d0c44cda5cb7953b2a4f`). It contains one FASTA file for each of LGI–LGVII. Each file has six same-length sequences: `2508@6`, `2509@6`, `Ncra@999`, `Ndis@999`, `8817@999`, and `1135@999`. Thus the *N. discreta* row is actually present. The `Ndis@999` sequence is the outgroup record in this alignment; the paper's phylogeny identifies the *N. discreta* strain as 8879. Full file sizes, checksums, record lengths, and sequence-character counts are in `outgroup_mapping_validation.json` and `outgroup_alignment_inventory.tsv`.

## Coordinate mapping

Dryad metadata says alignment coordinates follow the 2509 *N. tetrasperma* mat-a reference. In all seven files, the `2509@6` row has no gap characters and its sequence length is the linkage-group coordinate length. The deposited genotype matrix also has a `2509@6` sample column. At every genotype coordinate on every linkage group, the aligned 2509 base at column `POS-1` exactly equals the matrix `2509@6` call: 2,265,369/2,265,369 matches and zero mismatches. The alternative `POS+1` alignment-column interpretation matches only about 26.6–26.9% of compared bases. This validates the 1-based coordinate convention and the direct coordinate-to-column mapping against the matching reference sample across all SNP positions.

Because `2509@6` contains no alignment gaps, alignment column `POS-1` maps directly to reference position `POS`; gaps in the `Ndis@999` row represent an unavailable/unaligned outgroup base at that reference coordinate and are excluded. The coordinate-axis mapping works on all seven LGs, including 301,014 callable *N. discreta* bases among 331,711 genotype rows inside LGI 1.5–8.0 Mb and 149,665 among 174,635 outside. LGII–LGVII together have 1,474,160 callable bases among 1,759,023 genotype rows. Remaining sites are split into alignment gaps and ambiguous characters in `outgroup_callability.tsv`.

The validation is anchored to the 2509 reference row and to the matching deposited `2509@6` genotype column, not to coincidentally similar numeric positions. No separate de novo alignment was inferred. The deposited alignment is the source of homology; within it, any Ndis gap or ambiguous symbol is excluded. The 2509 reference row has no gaps, so there are no reference insertion columns requiring a nontrivial indel lift. Inversion/structural blocks are represented in the source alignment; no cross-chromosome or name-only joins were used. The high Ndis gap rate is explicit missing alignment coverage, not a callable ancestral base.

## Callability summary

| Region | Genotype SNP rows | Ndis A/C/G/T | Alignment gap | Ambiguous |
|---|---:|---:|---:|---:|
| LGI inside 1.5–8.0 Mb | 331,711 | 301,014 | 30,635 | 62 |
| LGI outside 1.5–8.0 Mb | 174,635 | 149,665 | 24,944 | 26 |
| LGII–LGVII | 1,759,023 | 1,474,160 | 241,802 | 43,061 |

The counts concern already-ascertained genotype SNP rows only, not invariant callable genomic bases. They do not constitute an unbiased genome-wide missingness estimate.
