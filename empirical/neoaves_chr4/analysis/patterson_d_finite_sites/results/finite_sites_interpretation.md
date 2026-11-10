# Finite-sites interpretation

## Numerical result

| Region | observed D | JC69 | HKY85 | GTR+Gamma | observed ABBA/BABA | GTR expected ABBA/BABA |
|---|---:|---:|---:|---:|---:|---:|
| all PNAS | -0.25819389863468784 | -0.3957154161839579 | -0.34750840847103087 | -0.291682924993958 | 9155.0/15528.0 | 2385.5/4350.2 |
| outside | -0.0805788766970934 | -0.09353633294231861 | -0.09170065822412289 | -0.0929047555333484 | 33894.0/39835.0 | 8888.3/10709.0 |
| inside minus outside | -0.177615 | -0.302179 | -0.255808 | -0.198778 | | |

The near−1 infinite-sites predictions arise from the extreme inferred discordant-topology composition inside the PNAS mask. Complete finite-sites pruning, including terminal branches, recurrent substitutions, and the empirical base composition, substantially attenuates this ratio. GTR+Gamma predicts the N61 shift at approximately −0.1988, versus observed −0.1776; JC69 and HKY85 predict −0.3022 and −0.2558. Thus terminal-branch/homoplasy effects explain much of the mismatch, while model choice changes the residual modestly.

Expected absolute counts remain lower than observed counts, so agreement in D is not a full likelihood fit. The trees and alignments are not independent: the published gene trees were estimated from these alignments. The original inference model and a faithful parametric re-estimation experiment were not recovered, so gene-tree error cannot be quantified here.

Columbea remains positive and N62 near zero, unlike N61. These results are compatible with N61-specific genealogy asymmetry, but do not distinguish MSRC from introgression, ancestral structure, or other causes. PNAS intervals are the published GalGal4-derived mask applied to the Stiller GalGal6 axis; exact boundary lifting remains unresolved.

See `statistical_validation.tsv` for 500-kb, 1-Mb, and 2-Mb block bootstrap intervals.
