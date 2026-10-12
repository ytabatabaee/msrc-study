# Final targeted fire-ant recombination recovery

## Classification

`QUANTITATIVE_RECOMBINATION_TRACK_NOT_RECOVERED`

The final targeted attempt recovered a validated historical assembly placement
but no validated marker-level physical coordinates on the frozen `Si_gnGA`
axis. No cM/Mb estimate, interval-level rate, window-level rate, association
test, or recombination figure is therefore produced.

## Frozen analysis audit

The prior recovery state was audited at commit `ed101dcc5040ca3ab3fed5281b49a74bda9fe5cb`.

- Reference: Stolle `Si_gnGA`, `gng20170922wFex.fa`; target `chr16`.
- Frozen supergene: `chr16:11,680,438–27,917,498`.
- Existing windows: 213; chromosome-16 windows: 96.
- Wang marker rows retained from Stage 7: 27934; total marker/pairwise rows: 27937; target chr16 mappings: 0.
- Existing Stage 4A–6C topology, discordance, and branch-length outputs were not modified.

## What the targeted correspondence recovered

Pracana et al. (2017) used seven Wang linkage maps to order and orient
scaffolds from the older `Si_gnG`/`GCA_000188075.1` reference with ALLMAPS.
The cited public supplement contains methods and population-genetic tables,
not the ALLMAPS input files, scaffold-order table, or an AGP. The cited Dryad
record is `10.5061/dryad.js509`.

Stolle et al. (2019) provide the closest recoverable bridge: the archived
`S.inv.AGP.txt` orders and orients 31 `GCF_000188075.1` scaffold
components on an optical `chr16` pseudomolecule. This is preserved in
`data/fire_ants_chr16/processed/final_historical_chr16_scaffold_mapping.tsv`.
The table deliberately leaves Si_gnGA coordinates blank. Its scaffold names
are accessioned `NW_...` identifiers, whereas every Wang row uses a
`Si_gnF.scaffold..._nt...` identifier. The available records do not establish
that identifier correspondence, and no sequence-based transfer to
`gng20170922wFex.fa` was available. Optical chromosome `chr16` was therefore
not silently treated as the current Si_gnGA `chr16`.

## Wang et al. (2013)

The Wang tables contain linkage-group labels, marker identifiers, original
Si_gnF scaffold/offset strings, and family-specific cM positions. They contain
no marker sequences. The historical AGP can establish an optical linkage-group
context, but it cannot lift the Wang Si_gnF marker offsets without a validated
Si_gnF-to-GCF/Si_gnGA correspondence. Consequently, zero markers were newly
mapped to target chr16, including zero inside and zero outside the supergene.
The earlier zero-marker result was therefore an assembly/identifier
correspondence failure, rather than evidence that Wang measured no chr16
markers.

## Ross and Shoemaker (2018)

Ross and Shoemaker did not publish a 2,796-marker RAD-seq map. The 2,796
RAD-marker figure belongs to Wang's P034 linkage map. Ross and Shoemaker used
a 15-locus pedigree panel and reported pairwise segregation/recombination
measurements for the three supergene-linked loci `C294`, `Gp-9`, and `i_126`
among polygyne queen progenies. These measurements remain pairwise-only:
their physical positions were not validated on Si_gnGA and the public raw
workbook could not be retrieved from Dryad in this environment. They cannot
be converted to cM/Mb or assigned to windows without unsupported assumptions.

## Quantitative estimates and integration

No adjacent marker pair has both validated Si_gnGA physical coordinates and
genetic positions. Therefore the rate formula was not applied. No marker
spacing, interval length, confidence bound, genotype-specific map, or
SB/SB-versus-SB/Sb comparison is reported. All 213 window-level
recombination values remain missing; 96 chr16 windows have
zero supporting markers. Quartet discordance, topology support, and branch
lengths are retained as existing measurements but are not tested against
recombination because the required spatially resolved quantitative covariate
is absent. No continuous track or association figure is scientifically
supported.

## Remaining blockers

1. Original Wang marker sequences or a validated Si_gnF-to-GCF/Si_gnGA
   scaffold correspondence is unavailable.
2. Pracana's described ALLMAPS inputs and intermediate scaffold-order files
   are not present in the cited supplement/data record recovered here.
3. Ross's pedigree panel supplies useful pairwise evidence but not a
   validated physical map on the frozen axis; its raw workbook was not
   accessible.
4. The historical optical chr16 AGP is an assembly correspondence artifact,
   not a direct meiotic recombination map.

The final-recovery source manifest, historical mapping table, explicit
missing-value window table, and regression checks are self-contained in this
repository. Execution: `python3 empirical/fire_ants_chr16/scripts/08_final_recombination_recovery.py`.
