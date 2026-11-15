# Stage 0 data audit — Channel Island monkeyflowers

## Data used

The locally supplied files were found in `empirical/monkeyflowers/data/raw/`; no files were redownloaded. The available coordinate-bearing files are the 100-kb window metadata TSV, the combined window/weight CSV, and the nonoverlapping coordinate grid. A fourth file contains the 15 TWISST topology counts but no coordinates. Exact byte sizes and SHA-256 checksums are in `input_manifest.tsv`.

The source study is Short & Streisfeld (2025), [Molecular Ecology 34:e17778](https://doi.org/10.1111/mec.17778). Their methods state that reads were aligned to the *M. aurantiacus* reference assembly from Stankowski et al. (2019). Window labels are `LG1`–`LG10`; the supplied start/end values form 1-based inclusive, nonoverlapping 100-kb intervals. The coordinate grid, metadata TSV and combined weight export have identical unique `(scaffold,start,end)` key sets (1,948 each), with no within-scaffold coordinate gaps. The combined file's `row_number` is 1–1,948 and its metadata columns agree with the TSV by coordinate key.

The combined file contains normalized weights and coordinates. The count-only weight file has 1,948 rows and no coordinate/window ID fields. Each count row, after normalization, matches the corresponding combined export weights across all 15 classes; this check necessarily uses the count file's row order because that file has no keys. Analyses use the combined coordinate-bearing export after this validation, and independently join that export to the metadata and coordinate grid using genomic keys.

All 1,948 windows have positive site counts (1,688–12,840), finite weights, and weights summing to one (maximum absolute sum error <1.6e-15). The raw count table and combined weight table are duplicates of the same topology information in different encodings; the combined export additionally joins coordinates and window metadata. The coordinate-grid file supplies the same scaffold/start/end keys and is retained for independent coordinate validation.

## Taxa and topology mapping

The tree weights refer to four ingroup taxa and *M. clevelandii* as outgroup: `ari` = *aridus*; `par` = *parviflorus*; `isl_lon` = Island *longiflorus*; `lon` = mainland *longiflorus*; `cle` = *M. clevelandii*. The TWISST file explicitly labels 15 rooted topologies. Following the article's definitions:

- Species history: `topo3`, `(cle,((ari,par),(isl_lon,lon)))`.
- Recent-introgression set: `topo9`, `topo12`, and `topo13`, each with a `par` + `isl_lon` sister pair. The set weight is their sum.
- Ancient-introgression topology: `topo10`, `(cle,(ari,(par,(isl_lon,lon))))`.
- The other 11 rooted topology classes remain separately available in `topology_weight_tracks.tsv`; their total is reported as `other_weight`, not redistributed into the focal three classes.

The authors describe the recent history as gene flow between Island *longiflorus* and *parviflorus* after mainland/Island *longiflorus* split; the ancient history is gene flow between the ancestor of *parviflorus* and the ancestor of clades C/D. Their interpretation is based on simulation comparisons and other evidence, and does not make every local alternative topology uniquely attributable to introgression.

## Local trees and scope

The supplied local Newick trees are not present in the repository's monkeyflower data or analysis folders. Consequently, this audit verifies topology-weight coordinates, not a one-tree-per-window set of Newicks. The topology weights cannot provide A(w), B(w), quartet internal branch lengths, substitution-unit branch metrics, CU estimates, or independent gene-tree species-tree inference. The completed work is a topology-only spatial analysis.

## Coordinate coverage

| Scaffold | Windows | Coordinate range (bp) |
|---|---:|---:|
| LG1 | 205 | 1–20,500,000 |
| LG2 | 243 | 1–24,300,000 |
| LG3 | 207 | 1–20,700,000 |
| LG4 | 185 | 1–18,500,000 |
| LG5 | 255 | 1–25,500,000 |
| LG6 | 157 | 1–15,700,000 |
| LG7 | 175 | 1–17,500,000 |
| LG8 | 202 | 1–20,200,000 |
| LG9 | 145 | 1–14,500,000 |
| LG10 | 174 | 1–17,400,000 |

No chromosome boundaries are crossed by any spatial runs. The primary chromosome 5 comparison uses the two independent 100-kb windows reported at 24.7–24.9 Mb; because the paper calls the candidate region approximately 300 kb without nucleotide-precise endpoints in its text, 24.7–25.0 Mb is an explicitly labeled sensitivity interval, not a newly asserted published breakpoint.
