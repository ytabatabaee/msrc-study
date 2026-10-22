# PNAS 48-taxon source reconstruction

Date: 2026-10-22. No ASTRAL-IV inference or SLURM job was run.

## Determination

The original collapsed archive can be assigned to genomic loci with a strong
but incomplete provenance chain. The 63,430-row published header supplies the
`Gene` index and genomic coordinates. Its 63,430 derived locus IDs are all
present exactly once in `63430.named.gene.trees.gz`. Comparing collapsed record
`i` with the named tree for header `Gene=i` gives identical taxon sets for
**63,113** records. This explains why the previous
same-row comparison was misleading: the named archive changes order after the
first block of loci, producing 31,066 apparent disagreements when the two
files are compared by physical line number.

The remaining **317** records have different taxon
sets from their header-derived named locus. Their identities are retained as
unresolved. No locus ID is assigned from topology similarity or from a guessed
row correspondence. Therefore the complete published collapsed collection is
not an exact verifiable input for locus-filtered 48-taxon treatments.

## Source evidence

| source | checksum / version |
|---|---|
| `63430.aLRT-0.95-collapsed.gene.trees.gz` | `5b058690d31531b3eb4f287556597bf5aa7db380ddf74c84bc53f71b512497df` |
| `63430.named.gene.trees.gz` | `cbf5a4f631fa92ef6da65c48a5aaf39b424fce94d2812c526eaed904dc61392b` |
| `master_table_gene_trees.txt` | `e6d01f4ecd9811f525da07a8702f811804998f135342d94a46a2f670247332f8` |
| PNAS repository | `17b408a622c9589fe9f6ca8a434c150974e39cd4` |
| PNAS Zenodo archive | `chr4avian-v1.0.0.zip`, `0dc24ebde4dd904115bd1e56735550e792d5cf8f7af437dbadf809108c313763` |
| Stiller Figshare | `10.6084/m9.figshare.25825408.v1` |
| Stiller ERDA | `https://sid.erda.dk/sharelink/ENhZODU9YE` |

The Stiller README explicitly identifies the collapsed aLRT-0.95 archive as
the ASTRAL input and describes the named archive as locus-name plus Newick
records. The PNAS `delta_quartet.R` script independently establishes the
inclusive GalGal4 window-start predicate used for outlier selection.

The exact producer script for the 63K collapsed archive is not included in the
available source snapshot. The archive itself uses decimal support labels and
multifurcations, but that does not establish branch-length or support
serialization details well enough to regenerate the file from the named
archive. No threshold-tuning reconstruction was used.

## Counts and exclusions

- Collapsed records: **63,430**.
- Unique header-derived loci: **63,430**.
- Taxon-set matches: **63,113**.
- Unresolved collapsed records: **317**.
- Published outlier loci under the frozen predicate: **1,431**.
- Existing J48 inputs preserved: **7**.

The 1,431-locus outlier mask is verified against the header-derived locus
coordinates and remains unchanged. This does not prove that the same 1,431
IDs can be selected from the unresolved collapsed records.

## Dataset decision

No `PNAS48_T0` or outlier-removed exact input is created. Doing so would
silently turn the 317 unresolved records into genomic loci. Existing `J48_*`
files remain the reconstructed named-tree re-collapse datasets. They are
appropriate for the MSRC extension after being labeled as reconstructed
design inputs, but they should not be described as exact PNAS input
reproductions.

The source archive is usable for a future **partial verified** dataset made of
the 63,113 matched records, but that dataset would omit
317 loci and must be labeled as such. It is not produced by this audit.

The 48-taxon mapping is the existing frozen `jarvis48_taxa.tsv` mapping with
48 source labels and 48 Stiller labels. The old `J48_*` counts and checksums
are recorded in `pnas48_dataset_comparison.tsv`; no corrected counterpart is
claimed because the unresolved records prevent exact locus filtering.

## Published species-tree outputs

The available `score-columbea/` files are Figure S9 scoring outputs and the
available `score.tre` is a full published tree; the archive does not identify a
48-taxon full-data or outlier-removed species tree. No new species tree was
inferred here. The 363-taxon published reference remains a reference output,
not validation of a 48-taxon reproduction.

The source reconstruction regression suite passes 4/4 tests. It checks the
63,430 unique header locus IDs, both source checksums, the 1,431-locus frozen
outlier mask, and preservation of the seven existing J48 manifest entries.

## Recommended next steps

Use the existing reconstructed-design `J48_T0`, `J48_T_PNAS`,
`J48_T_STRUCT_PRIMARY`, `J48_T_BLOCK_500K`, and `J48_NONCHR4` inputs for the
MSRC extension only, with their current checksums. For a strict PNAS replication,
first obtain a source-side file or author confirmation that links all 63,430
collapsed records to locus IDs, or explicitly accept the 317-locus omission as
a partial verified sensitivity. Then prepare Bridges-2 single-treatment jobs
without changing the 363-species inputs; do not submit them from this audit.
