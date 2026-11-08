# Microbotryum Stage 0 feasibility audit

## Decision

**NO-GO for this run; conditional follow-up only after Dryad file access is restored.**

The public metadata identifies the requested `gene_trees` (59,490 bytes) and
`speciestree.tre` (1,619 bytes), but the actual files could not be downloaded
in this environment. The Dryad API metadata request succeeded; the file
download endpoint returned `401 Unauthorized: must have current bearer token`,
and the landing-page `file_stream` endpoint returned HTTP 403. Therefore no
gene-tree count, taxon inventory, branch-length assessment, stratum count, or
quartet result is claimed.

## Sources inspected

- Lucotte et al. (2025), *Nature Communications*, DOI
  [10.1038/s41467-025-60222-5](https://doi.org/10.1038/s41467-025-60222-5).
- Dryad DOI [10.5061/dryad.0k6djhb89](https://doi.org/10.5061/dryad.0k6djhb89).
- Authors' repository [MicrobotryumOnDianthus](https://github.com/eliselucotte/MicrobotryumOnDianthus).
- Authors' code archive [Zenodo 10.5281/zenodo.14918233](https://doi.org/10.5281/zenodo.14918233).

Dryad metadata confirms that the deposited record contains `gene_trees`,
`speciestree.tre`, README material, assemblies, annotations, and feature-count
files. The small gene-tree and species-tree file checksums are recorded in
`data/microbotryum/metadata/source_manifest.tsv`.

## Coordinate and stratum recoverability

The authors' code was inspected. `TP/2a.make_code_to_group_trees.py` reads a
separate stratum-dating table and maps orthogroups to strata, while
`TP/3.Awk_command_parse_genetrees.sh` parses IQ-TREE outputs and rewrites tip
labels. The code therefore documents an orthogroup-to-stratum route, but it
does not by itself establish physical chromosome coordinates for every
deposited tree. The deposited gene-tree file and the source stratum table are
required before coordinate recovery can be classified as verified.

The audit deliberately records physical coordinates, gene order, and stratum
assignments as unresolved. No position is inferred from topology or from a
stratum label.

## Recombination suppression

The publication describes mating-type chromosome linkage and evolutionary
strata as biological evidence for historical recombination suppression. That
is independent structural/evolutionary evidence, not a direct local crossover
rate measurement. Region boundaries and control regions cannot be audited
from the inaccessible gene-tree and annotation resources. The existing
reporting files therefore do not assert coordinates or quantitative rates.

## Exploratory quartet analysis

Not run. The intended fixed quartet is
`S1_a1,S1_a2,S2_a1,S2_a2`, with q1 = species clustering, q2 = mating-type
clustering, and q3 = the third resolution. The parser and fixed-label tests
are implemented, but there are zero verified trees in this run.

## What is needed for a conditional GO

Restore authenticated or otherwise permitted access to the two small Dryad
files, then rerun `stage0_audit.py`. A complete feasibility decision requires:

1. actual Newick record counts and taxon/mating-type labels;
2. branch-length and resolution completeness;
3. exact orthogroup-to-stratum mapping;
4. independently verified gene order or physical coordinates; and
5. trees from recombining control regions.

Until those checks are complete, this dataset should not be added as a sixth
empirical example.
