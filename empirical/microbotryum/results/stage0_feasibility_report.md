# Microbotryum Stage 0 feasibility audit

Sources: [Lucotte et al. 2025](https://doi.org/10.1038/s41467-025-60222-5), [Dryad dataset](https://doi.org/10.5061/dryad.0k6djhb89), and the [authors' analysis repository](https://github.com/eliselucotte/MicrobotryumOnDianthus).

## Decision

**CONDITIONAL GO:** the deposited gene trees support a stratum-level genealogical analysis, with a usable PAR control category; physical positions are only verified on an MPbarcode50 scaffold and no quantitative crossover-rate track is available.

## Verified data

The local Dryad files contain 260 parseable gene trees. All usable trees have six taxa, complete branch lengths, and resolved bifurcating topology under the parser. The first column is a combined `stratum-gene` identifier. Deposited strata and counts are: `HD_RR` (9), `PARs` (29), `black` (134), `blue` (17), `dark_ruby` (25), `in_HD` (2), `light_ruby` (13), `orange` (4), `purple` (3), `red` (2), `turquoise` (22).

Taxa are `Mv-cat-C212-A1/A2`, `Mv-sup-1065-A1/A2`, and `Mv-lag-1253-A1/A2`; the first four provide the focal quartet. The deposited A1/A2 labels are retained exactly and are not inferred arrangement states.

## Quartet results

q1 is species clustering `(Mv-cat-A1,Mv-cat-A2)|(Mv-sup-A1,Mv-sup-A2)`, q2 is mating-type clustering `(cat-A1,sup-A1)|(cat-A2,sup-A2)`, and q3 is the third resolution. Labels were fixed before examining stratum frequencies.

| stratum | trees | q1 | q2 | q3 | B | internal branch | A |
|---|---:|---:|---:|---:|---:|---:|---:|
| HD_RR | 9 | 1.000 | 0.000 | 0.000 | 0.03056 | 0.04360 | 0.04360 |
| PARs | 29 | 1.000 | 0.000 | 0.000 | 0.07030 | 0.10437 | 0.10437 |
| black | 134 | 0.119 | 0.873 | 0.007 | 0.14770 | 0.05926 | -0.02324 |
| blue | 17 | 1.000 | 0.000 | 0.000 | 0.28713 | 0.15573 | 0.15573 |
| dark_ruby | 25 | 0.600 | 0.400 | 0.000 | 0.25593 | 0.31392 | 0.29695 |
| in_HD | 2 | 0.000 | 1.000 | 0.000 | 0.08562 | 0.04156 | -0.02078 |
| light_ruby | 13 | 1.000 | 0.000 | 0.000 | 0.26830 | 0.37704 | 0.37704 |
| orange | 4 | 0.250 | 0.750 | 0.000 | 0.10283 | 0.04266 | -0.01242 |
| purple | 3 | 0.000 | 1.000 | 0.000 | 0.15853 | 0.09364 | -0.04682 |
| red | 2 | 1.000 | 0.000 | 0.000 | 1.77099 | 0.84642 | 0.84642 |
| turquoise | 22 | 0.955 | 0.045 | 0.000 | 0.04997 | 0.02626 | 0.02436 |

Frequencies are across genes within each stratum, not probabilities inferred from one tree. Distances are substitution-tree units, not coalescence times. These results establish a stratum-dependent genealogical pattern but do not establish MSRC as a unique explanation.

## Older versus younger strata

The published study identifies the black region as an older trans-specific nonrecombining stratum and turquoise as a younger, species-specific stratum. In these trees, black has q2 = 0.873, B = 0.14770, internal branch = 0.05926, and A = -0.02324; turquoise has q2 = 0.045, B = 0.04997, internal branch = 0.02626, and A = 0.02436. This is a strong stratum-dependent contrast, but the differences are descriptive and may reflect mutation rate, alignment, or gene-tree estimation effects.

## Coordinate and control-region audit

The gene identifiers map to 14876 unique genes in `MPbarcode50_tsebra.gtf`; all 260 unique tree gene identifiers are represented. The mappings provide contig, start, and end coordinates on MPbarcode50; they are scaffold/contig coordinates, not chromosome-level coordinates. Gene order is therefore available within contigs but not as a validated chromosome-wide track. The identifier match does not independently prove that the MPbarcode50 annotation and the Mv-sup-1065 gene-tree reference use the same sample/assembly, so these coordinates remain provisional.

The deposited `PARs` category contains 29 usable trees and is the available pseudo-autosomal/recombining control category. It is a categorical control, not a quantitative crossover-rate track. Stratum labels and PAR assignment are independent evidence of recombination regime, not direct local rate measurements.

## Recommendation

The dataset is suitable for a complete stratum-level exploratory analysis and for testing species-versus-mating-type genealogy shifts. It remains **CONDITIONAL GO** rather than a sixth coordinate-indexed MSRC example until chromosome-level gene order and the MPbarcode50-to-tree sample correspondence are independently recovered. A follow-up can use the PARs category as a control and should validate ordering of MPbarcode50 contigs on the mating-type chromosome.

Raw input directory: `/Users/ytabatabaee/Desktop/msrc-paper/data/microbotryum/raw`.
