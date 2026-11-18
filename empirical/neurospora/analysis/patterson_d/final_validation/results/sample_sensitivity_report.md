# Isolate-level sensitivity

## LGI primary result

Leave-one-isolate-out estimates across the 10 mat-a and 8 mat-A samples are tightly clustered: D_inside ranges from +0.87116 to +0.87198 and D_outside from +0.43992 to +0.45795. The inside-minus-outside direction is unchanged in every leave-one-out run. The outside estimate is somewhat more influential than the inside estimate, but no single isolate explains the LGI contrast.

Balanced random subsets were selected without using D values, with seed 20261118. Twenty independent draws were made for each 2/2, 4/4, 6/6, and 8/8 design on LGI. Inside D remained approximately +0.870–+0.872 at all sizes. Outside D was more variable for small groups: ranges were approximately +0.336 to +0.525 for 2/2, +0.402 to +0.484 for 4/4, +0.428 to +0.468 for 6/6, and +0.441 to +0.452 for 8/8. This supports a robust LGI regional contrast in the verified sample set, while showing expected small-sample noise in the lower-coverage complement.

## LGIV signal and matched-pair discrepancy

LGIV is not sample-robust in the same way. The full 10/8 groups have D=+0.0715, but the two verified opposite-mating-type pairs combined have D≈−0.740. In 10 balanced random draws at each size, LGIV values span roughly −0.69 to +0.75 at 2/2, −0.66 to +0.64 at 4/4, −0.27 to +0.55 at 6/6, and −0.02 to +0.21 at 8/8. Leave-one-out changes are largest for mat-a CJ06/CJ16 and mat-A CJ09/CJ19; deleting CJ16 moves the full-group estimate to about +0.239, while deleting CJ19 moves it to about −0.122.

The four one-isolate combinations from the two documented pairs are also heterogeneous. CJ15(mat-A)/CJ16(mat-a) gives D≈−0.797 across 19,724 informative sites; the CJ03/CJ04 pair has only four informative sites and its +0.5 D is not informative. The two cross-pair combinations range from near zero to negative with different site support. Thus the combined pair result is dominated by a highly asymmetric individual combination, not replicated by two equivalent informative pairs. The available data do not identify whether this comes from isolate-specific history, geography/structure, missingness, or variant ascertainment. It contradicts any claim that sample behavior is uniform genome-wide.

The L10 individuals are described as UK in Corcoran et al.; finer geography for these isolate-level patterns is not established in the local metadata used here. No mating-type labels were inferred from isolate IDs or trees; CJ07, CJ11, CJ12, CJ18 remain excluded.

## Reproduction

The per-isolate and per-draw table is `sample_sensitivity.tsv`. Its entries are sensitivity analyses over shared linked sites, not independent demographic replicates or formal p-values.
