# Stage 0 — Neurospora source and tree audit

## Provenance

The repository started clean on `main`. `data/neurospora/raw/neurospora.tar.gz` was absent, so the public Zenodo mirror was downloaded from the [record for Martin & Van Belleghem's TWISST data](https://zenodo.org/records/4946041). The deposited archive is 15,577,951 bytes and its MD5, `571324702709eb1c2732a3fd993edc93`, matches the record. Raw and extracted files are excluded from Git. Compressed-file sizes and checksums are in `source_manifest.tsv`.

The archive README-equivalent description and actual contents identify genotype data, per-window Newick trees, window-coordinate metadata, four-lineage quartet-resolution definitions, and one published three-resolution TWISST weight file. The main files are:

| File | Observed contents |
|---|---|
| `neurospora/trees/*.data.tsv.gz` | One row per 50-variable-site window: linkage group, start, end, midpoint, and site count |
| `neurospora/trees/*.trees.gz` | 45,303 window records: 32,472 Newick trees with branch lengths and 12,831 literal `NA` records |
| `neurospora/cras_his_L10.*.weights.tsv.gz` | Three topology subtree counts per window; complete counts normalize by their row sum (120) |
| `neurospora/cra_his_L10.topos` | q1 species-associated, q2 remaining, q3 N. hispaniola + L10 mat-A topology definitions |
| `neurospora/geno/*.geno.gz` | Genotype matrix whose header includes 22 L10 CJ samples and four additional reference samples |

The tree and coordinate files each contain 45,303 records and align one-to-one. There are 32,472 actual Newick trees; 12,831 corresponding rows are literal `NA` in both trees and weights. The 10,126 LGI records span coordinates 752–9,464,227 bp; the remaining 35,177 records cover LGII–LGVII. Coordinate intervals do not overlap within linkage groups. Each window represents 50 variable sites; physical widths vary. Within the approximately 1.5–8.0 Mb LGI interval, 6,634 window records are represented and 5,767 have both a local tree and valid weights.

## Taxa and branch lengths

The trees use CJ strain labels for the L10 *N. tetrasperma* samples and labels `Ncra@999`, `8817@999`, `1135@999`, and `2509@6` for the additional samples. The TWISST topology definitions explicitly group the L10 samples as `L10_A` and `L10_a`, with *N. crassa* and *N. hispaniola* as the other two quartet lineages. Martin & Van Belleghem describe 10 mat-A and 12 mat-a L10 samples; Corcoran et al. identify CJ03/CJ04 and CJ15/CJ16 as matched mat-A/mat-a pairs. The trees themselves do not encode mating-type names on each CJ tip, so individual assignments are not inferred from labels alone. Topology estimates use the authors' deposited group assignments; branch-length contrasts use those independently documented matched pairs.

The trees are BIONJ output from PhyML GTR windows, treated as unrooted for quartet analysis. Their branch lengths are substitution distance per site, not coalescent units. The per-tree inventory records taxon-set IDs, missing non-root branch lengths, and multifurcations; the sole 26-tip set is listed once in `tree_taxon_sets.tsv`. The trees are haplotype/strain genealogies. Treating each of 22 L10 samples as an independent species would be biologically incorrect; the analysis therefore compares four lineage groups via deposited quartet weights rather than presenting an individual-tip ASTRAL tree as a species tree.

## Biological boundary and history

Martin & Van Belleghem report that q3, grouping L10 mat A with *N. hispaniola*, replaces q1 across the central approximately 7 Mb of LGI. They interpret this in light of Corcoran et al.'s population-genomic evidence for introgression into L10 mat-A chromosomes and strong recombination suppression on the mating-type chromosome. The primary interval used here is approximately 1.5–8.0 Mb on the study coordinate axis; the 250-kb flanks are a sensitivity comparison, not independently asserted breakpoints. Independent work documents large inversions on the mating chromosome, but recombination suppression is not reducible to structural heterozygosity alone. Introgression and suppressed recombination are treated as interacting processes.

Primary sources: [Corcoran et al. 2016, Genome Research](https://pmc.ncbi.nlm.nih.gov/articles/PMC4817772/); [Martin & Van Belleghem 2017, Genetics](https://pmc.ncbi.nlm.nih.gov/articles/PMC5419486/); [Corcoran et al. 2017, Nature Communications](https://pmc.ncbi.nlm.nih.gov/articles/PMC5658415/).
