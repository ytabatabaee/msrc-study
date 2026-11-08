# Microbotryum Stage 2 coordinate validation

## Mapping result

The deposited two-column `gene_trees` file contains 260 records with prefixes such as `black-g2118`; the Newick tips retain only six taxon labels. No orthogroup identifier is retained in these records. The local `Orthogroups_filtered.txt` contains full sample-specific gene identifiers and orthogroup rows, but there is no deposited table linking the `gNNNN` tree prefix to an `OG` row. Therefore orthogroup recovery is **not verified**.

The numeric suffixes in all 260 tree records map to an `Msup_*_gNNN` identifier in `MPbarcode50_tsebra.gtf`; this is an identifier-suffix match, not proof of orthology or sample identity. Mapping counts are: 260 contig-local provisional matches, 0 chromosome-scale coordinates, 0 orthogroup-resolved records, and 0 unmapped records.

## Assembly/sample validation

The paper identifies the focal phased assemblies as M. superbum strain 1065 and M. shykoffianum strain c212. The local README describes `MPbarcode50.fa` only as a haploid Flye ONT assembly; its anonymous `contig_*` headers and `Msup_tig...` annotation namespace are compatible with an M. superbum annotation, but do not identify the isolate or mating type. The article describes a submitted de novo annotation for M. superbum strain 6P, whereas the gene-tree taxa are Mv-sup-1065-A1/A2 and Mv-cat-C212-A1/A2. This prevents treating MPbarcode50 as the 1065 tree assembly without an explicit submitted-file cross-reference.

## Physical order and boundaries

The publication uses M. lagerheimii as a proxy for ancestral mating-type-chromosome gene order and documents the PR chromosome fusion, edge PARs, and colored strata. It does not provide a chromosome-coordinate transform in the local deposited files. The current MPbarcode50 matches span several anonymous contigs, including multiple strata sharing contigs. No chromosome-scale order, validated gene-order rank, or base-pair stratum boundary was recovered. No spatial genealogy tracks are generated and no contigs are concatenated.

## Consequence

The Stage 2 spatial question remains blocked. The valid empirical comparison is stratum-level, with the deposited stratum label as the independent regime annotation and the PAR category as a categorical control. The complete one-row-per-tree mapping table is retained for future resolution.

References: [Lucotte et al. 2025](https://doi.org/10.1038/s41467-025-60222-5), [Dryad](https://doi.org/10.5061/dryad.0k6djhb89), [authors' repository](https://github.com/eliselucotte/MicrobotryumOnDianthus).
