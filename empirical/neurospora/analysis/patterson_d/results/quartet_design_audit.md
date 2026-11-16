# Quartet design audit

## Frozen primary design

`((P1=L10 mat-a, P2=L10 mat-A), P3=N. hispaniola), O=N. discreta`. P1 contains the 10 independently verified mat-a L10 individuals and P2 the 8 verified mat-A individuals; unresolved CJ07, CJ11, CJ12 and CJ18 are excluded. P3 is `8817@999`; O is `Ndis@999`. The exact sample list, groups, predeclared baseline and candidate alternative are in `config/frozen_quartets.tsv` and `config/frozen_sample_membership.tsv`.

The population-history baseline groups the two L10 mating-type backgrounds as sisters relative to the donor-associated *N. hispaniola* lineage. The alternative of interest is increased derived-allele sharing between mat-A and *N. hispaniola*; this hypothesis did not determine the outgroup choice. The external *N. discreta* outgroup is present in the chromosome alignments and is coordinate validated against the genotype reference at all 2,265,369 positions.

## Sensitivity and caveats

CJ03(A)/CJ04(a) and CJ15(A)/CJ16(a) are two documented matched isolate pairs. Their combined set is an explicitly limited sampling sensitivity, not a balanced full population design. Four CJ labels remain unresolved. Alignment gaps, ambiguity, outgroup polarization failure, multi-allelic patterns and insufficient calls are handled by the frozen site filters. The four-tip outgroup design cannot be exactly induced from the L10 published local trees because Ndis is absent; comparisons to TWISST q weights and available branch metrics are reported only as distinct, related-estimand spatial summaries.
