# PNAS 48-taxon data audit

Date: 2026-10-22. This is a data and provenance audit; no ASTRAL-IV inference
was run and no cluster job was submitted.

## Answers to the required questions

1. **Are the current J48 inputs identical to the published PNAS inputs?** No.
They are a reconstructed published-design dataset. The PNAS/Stiller README
identifies `63430.aLRT-0.95-collapsed.gene.trees.gz` as the ASTRAL input,
whereas the current Stage 4E inputs come from the ID-preserving
`named_recollapsed.tre` fallback. A provenance-aware comparison resolves the
archive-order problem behind most of the 31,066 naive same-row disagreements:
63,113/63,430 collapsed records match the header-derived named-locus taxon set,
while 317 records remain unresolved.

2. **What differs?** The input gene-tree collection/preprocessing is the
material discrepancy. Current pruning retains the 48 mapped Stiller labels,
omits trees with fewer than four retained taxa, and preserves branch lengths;
317 collapsed records still lack a verified locus assignment.

3. **Did we identify original published 48-taxon gene-tree files?** No.
The available PNAS archive contains scripts, support data, and species-tree
outputs, but no authoritative 48-taxon gene-tree file with locus IDs. The
published method describes restricting the 63,430 Stiller loci to the 48
Jarvis species. A strong partial mapping exists, but exact complete 48-taxon
files remain unavailable.

4. **Are PNAS excluded intervals and locus IDs correct?** The interval
predicate and local locus membership are exactly reproduced: 1431
loci, with zero implementation-only or published-only loci in the local
manifest. This validates the mask against the published script, not the
unavailable original 48-taxon file. The three intervals are GalGal4 `chr4`
window-start intervals and are inclusive.

5. **What was corrected?** No raw `.tre` was overwritten and no corrected
published-equivalent dataset was fabricated. The audit records
header-derived locus assignments and retains the 317 unresolved records
explicitly; current J48 inputs remain explicitly reconstructed inputs.

6. **What remains unverifiable?** Exact collapsed-tree-to-locus identity,
published 48-taxon gene-tree topology/branch/support equality, the original
48-taxon ASTRAL command/version/output, and a published 48-taxon full-data or
outlier-removed species tree.

7. **Which datasets should run on Bridges-2?** Run the reconstructed-design
baseline and MSRC treatments only after labeling them as reconstructed:
`J48_T0`, `J48_T_PNAS`, `J48_T_STRUCT_PRIMARY`, `J48_T_BLOCK_500K`, and
`J48_NONCHR4`. Keep the current files for provenance; do not call these exact
published-input reproductions.

## Taxa and filtering

The frozen mapping contains 48 source 2014 labels and 48 Stiller labels. The
only non-exact mapping status is recorded in the mapping file (`Manacus` was
reidentified in the Stiller source); no label was invented by this audit.
Current Stage 4E counts and checksums are listed below. Fewer-than-four-taxon
trees are omitted after pruning, as implemented by `05_stage4e_prepare_aster.py`.

| treatment | loci | omitted <4 | SHA-256 |
|---|---:|---:|---|
| J48_T0 | 62945 | 485 | `cb5195cc449582084a0afc1ce11c87cc8eaa176863f4682f78ad4a29eb22dc14` |
| J48_T_PNAS | 61524 | 475 | `6c60dce91603d46f5d0af75ec46093065cb605bcc8a8ce26e112d420482953b8` |
| J48_T_STRUCT_PRIMARY | 62116 | 479 | `000213a5e5628e1498a0dd8aebca889c056d789cde0d743c5a858d059ff15ef3` |
| J48_T_BLOCK_250K | 57084 | 451 | `f1e187276b8bf3783b15c947b74f604efaa32efa0056923d1234dc86a52056fa` |
| J48_T_BLOCK_500K | 56905 | 450 | `d9d6004ae80a6b17289ac57afe7941fc539d041eb76576bd0d6a27879d405357` |
| J48_T_BLOCK_1MB | 56817 | 448 | `a01a90932b564afbcdfc68f2e4e42343d79a2c3019178ed29805ed765a1cc4e9` |
| J48_NONCHR4 | 56721 | 447 | `cfecbc68da62fc6bd03b971f09bd95b34fa5efa5ff247cef3672ae90f6667760` |

## Replication-fidelity classification

| component | determination |
|---|---|
| Species selection | reconstructed match with documented label mapping |
| Source gene trees | mismatch / exact identity unverifiable |
| Gene-tree preprocessing | mismatch: current ID-preserving re-collapse vs published collapsed archive |
| Locus filtering | verified for the local reconstruction |
| PNAS outlier exclusions | exact published predicate and local locus membership verified |
| Species-tree inference settings | unverifiable for the 48-taxon comparison |
| Published output topology | unverifiable for a 48-taxon output |

This therefore qualifies as a **reconstructed published-design reproduction**.
The structural masks, block controls, and CU/SU comparisons remain **MSRC
extensions**.

## Sources and versions

- `pnas_repository`: https://github.com/smirarab/chr4avian/tree/master
- `pnas_repository_commit`: 17b408a622c9589fe9f6ca8a434c150974e39cd4
- `pnas_zenodo`: https://doi.org/10.5281/zenodo.10699423
- `pnas_zenodo_downloaded_record`: https://zenodo.org/api/records/10699424/files/smirarab/chr4avian-v1.0.0.zip/content
- `stiller_figshare`: https://doi.org/10.6084/m9.figshare.25825408.v1
- `stiller_erda`: https://sid.erda.dk/sharelink/ENhZODU9YE

Selected downloaded source checksums:

- `63430.aLRT-0.95-collapsed.gene.trees.gz`: `5b058690d31531b3eb4f287556597bf5aa7db380ddf74c84bc53f71b512497df`
- `63430.named.gene.trees.gz`: `cbf5a4f631fa92ef6da65c48a5aaf39b424fce94d2812c526eaed904dc61392b`
- `chr4avian-v1.0.0.zip`: `0dc24ebde4dd904115bd1e56735550e792d5cf8f7af437dbadf809108c313763`
- `delta_quartet.R`: `35e0e7094e01a957643efc001bdd5820320502ab4ffa8dd74485c2a67dfdab37`
- `genetreesupport-README.md`: `89a5166ed9abb1e4cb5c232f0f3a0ac049b8a696466c846c21000bb7141f937a`
- `master_table_gene_trees.txt`: `e6d01f4ecd9811f525da07a8702f811804998f135342d94a46a2f670247332f8`


Local reconstruction checksums:

- `data/stage4d/jarvis48_taxa.tsv`: `5bd5c663c65c5b1bad797f2d337e97a29b6400184c32a2d1e0892a012da37f64`
- `data/stage4d/named_recollapsed.tre`: `c0a0ec9b432405da548a6cdb37460ccf1c0d9d1b5a5211d74d5f397c70e0f943`
- `results/stage4d_locus_manifest.tsv`: `8290505033761276c458399ac2821773a77cb13a39d9679fe53c00f1605a23ac`

Published interval coordinates are GalGal4 window-start intervals
25,030,000–32,670,000, 33,510,000–34,470,000, and 44,130,000–56,810,000;
the paper's GalGal6 table gives corresponding approximate ranges
25,555,000–33,202,000, 34,230,000–34,999,000, and 44,690,000–57,179,000.
The local locus mask uses the published GalGal4 predicate because that is the
coordinate system used by the authors' script.

The local Stiller archive manifest records the ERDA archive files and checksums;
the local PNAS download manifest records the GitHub commit and Zenodo ZIP
checksum. Large raw archives remain external and are not duplicated by this
audit. No published 48-taxon input was found in the available archive inventory.

## Recommended command preparation

Prepare, but do not submit, one treatment per job using the existing single-
treatment workflow after transferring the current Stage 4E inputs:

```bash
for T in J48_T0 J48_T_PNAS J48_T_STRUCT_PRIMARY J48_T_BLOCK_500K J48_NONCHR4; do
  sbatch -A bio240015p --qos=rm -p RM-shared --cpus-per-task=32 \
    --time=16:00:00 --export=ALL,TREATMENT="$T" \
    empirical/neoaves_chr4/cluster/stage4e_aster/01_single_treatment.sbatch
done
```

These commands are intentionally recommendations only; no command was
submitted by this audit.
