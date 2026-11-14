# Scientific interpretation: Anopheles 2La nucleotide Patterson D

Frozen grid: 988 windows; observed records: 15808 (8 designs × 2 filters); invalid coordinates/counts/D range: 0; undefined D rows: 0.

## Findings

Pooled D uses the ratio of summed ABBA and BABA counts. Boundary windows are excluded from inside/outside contrasts; 2La is 20,524,058–42,165,532 (1-based inclusive). The four predeclared quartet designs are not interchangeable causal tests.
- **A_merus:** inside D=-0.7349 (n=180,425 informative site contributions), outside D=-0.0134 (n=200,142); inside−outside=-0.7215, 500-kb block-jackknife 95% CI [-0.7526, -0.6905].
- **B_merus:** inside D=0.8088 (n=185,700 informative site contributions), outside D=-0.0022 (n=160,736); inside−outside=0.8110, 500-kb block-jackknife 95% CI [0.7861, 0.8359].
- **C_merus:** inside D=-0.7028 (n=182,297 informative site contributions), outside D=0.0006 (n=188,237); inside−outside=-0.7034, 500-kb block-jackknife 95% CI [-0.7434, -0.6633].
- **D_merus:** inside D=0.8187 (n=213,555 informative site contributions), outside D=0.0124 (n=203,504); inside−outside=0.8063, 500-kb block-jackknife 95% CI [0.7782, 0.8343].

## Answers to the analysis questions

1. Inside-versus-outside contrasts are consistently nonzero under the block jackknife for the primary merus designs, but magnitude is highly design dependent; see `jackknife_results.tsv` at 250/500/1,000 kb for uncertainty sensitivity.
2. Direction depends strongly on quartet design: A contrasts are negative inside, while B contrasts are positive; C/D are sensitivity designs and should be interpreted with their smaller P2 groups.
3. The direction and inversion localization persist with A. melas in place of A. merus. Magnitude shifts modestly to moderately: the inside-minus-outside difference changes from -0.7215 to -0.6452 for A and from +0.8110 to +0.8330 for B. This is outgroup sensitivity, not independent replication.
4. MAC≥2 is a filtering sensitivity, not a guarantee of unbiasedness. It left pooled D and informative-site contributions unchanged to numerical precision in all reported region/design combinations, while reducing passing sites; see paired summaries and coverage figures.
5. The effects are distributed across the inversion rather than limited to one short segment: the sign is consistent in 99.5–100% of 50-kb windows for the primary merus designs, and each of four equal window-count inversion quartiles has a large pooled effect. Magnitudes vary along 2L; see `inversion_subregion_summary.tsv` and the resolution tracks, which preserve empty bins without interpolation.
6. Nucleotide D aligns strongly in space with Stage-2R topology-shift and other genealogy statistics in the nonoverlapping 500-kb descriptive comparison; see `topology_spatial_correlations.tsv`. Their definitions differ, so this is co-localization, not equivalence.
7. The design dependent pattern is compatible with heterogeneous genealogies and recombination partitioning under MSRC, but is not evidence uniquely supporting that mechanism.
8. Introgression, ancestral population structure, geographic structure among arrangement-defined samples, ancestral inversion polymorphism, and incomplete lineage sorting remain alternatives. Nonzero D alone does not prove introgression or MSRC.

## Limitations and validation

Grid keys/coordinates/count identities and D bounds were checked for all 15808 rows. Undefined-D output rows: 0; empty spatial bins remain explicit. Four selected windows were independently re-fetched and scalar-recomputed (24 design/filter comparisons; status PASS; maximum absolute discrepancy 8.53e-14). Site contribution sums count per-window pattern-bearing observations and can count genomic positions across different designs; they are not unique variant counts. Spatial correlations use nonoverlapping 500-kb blocks but are descriptive and not a formal spatial permutation test.
