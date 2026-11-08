# Microbotryum coordinate validation

## Finding

All 260 unique gene identifiers in the 260 deposited gene-tree records map numerically to `MPbarcode50_tsebra.gtf` (14876 unique GTF identifiers). This is useful identifier continuity, but it is not by itself proof that the annotation and the deposited trees use the same isolate, mating type, assembly, or orthology release.

The current files support **contig-local order only**. No chromosome-scale placement of MPbarcode50 contigs, no validated mapping to the M. superbum or M. shykoffianum mating-type chromosomes, and no independently verified base-pair boundaries for each color stratum were recovered from the local Dryad files or the inspected author scripts. No spatial tracks are therefore generated.

## Assembly and sample identity

`MPbarcode50.fa` has anonymous `contig_*` headers. Its annotation uses `Msup_tig..._g...` identifiers, indicating an M. superbum-oriented annotation namespace, but the local README does not state the exact isolate or mating type. The tree tips include `Mv-sup-1065-A1/A2`, `Mv-cat-C212-A1/A2`, and `Mv-lag-1253-A1/A2`; the MPbarcode50 sample correspondence remains unverified. The separate Mv-sup and Mv-cat FASTA resources use different sample-specific contig namespaces, so numerical gene-ID matching must not be treated as a chromosome mapping.

The published study states that mating-type chromosome contigs were identified using PR/HD homology and ortholog links, and that PARs occur at chromosome edges. It also states that gene trees were inferred for one-to-one orthologs and deposited with the species tree. Those statements independently support the stratum labels and biological interpretation, but do not supply a validated MPbarcode50 contig order in the files currently available here.

## Mapping confidence

| component | status | confidence |
|---|---|---|
| gene-tree Newick parsing | all 260 parse with six taxa and branch lengths | high |
| tree gene identifier to MPbarcode50 GTF identifier | 260/260 unique tree IDs represented | moderate for identifier continuity; low for biological sample correspondence |
| physical coordinates within an MPbarcode50 contig | available from GTF | moderate, annotation-dependent |
| chromosome-scale ordering | not recovered | none |
| stratum-to-boundary base-pair mapping | not independently verified | none |

## Consequence

The valid analysis unit is evolutionary stratum, not a continuous chromosome window. The PARs category is used as a categorical recombining control. Contig-local coordinates are retained in `coordinate_mapping_audit.tsv` but are not concatenated or plotted as a chromosome.

References: [Lucotte et al. 2025](https://doi.org/10.1038/s41467-025-60222-5), [Dryad](https://doi.org/10.5061/dryad.0k6djhb89), [authors' repository](https://github.com/eliselucotte/MicrobotryumOnDianthus).
