# Atlantic cod SNAPP time-tree structure analysis

This analysis summarizes the existing published 250-kb SNAPP MCC population time trees. It does not rerun SNAPP or infer a new joint genome-wide tree. A window is the spatial observation; neighboring windows are dependent. No window-matched posterior samples were available, so uncertainty below is spatial-null uncertainty, not posterior uncertainty.

## Input and tree audit

- Source rows: 426 across LG01, LG02, LG07, LG12; valid rooted ultrametric trees: 425; invalid/excluded: 1.
- All trees use the same 12 cod-population tips; no outgroup tips occur in the local MCC trees. Thus outgroup features do not drive this PCA, and cod-population-only is the primary analysis.
- Branch lengths are the published SNAPP time units (Ma); root stem, if present, is excluded from root age. Zero non-root branches are retained. Missing/negative lengths, non-ultrametric trees, or zero root age are not repaired.
- Per-window posterior files were not available in the source inventory. Supergene-wide posterior files are not matched to individual windows and are not treated as independent loci.
- Invalid tree rows: LG07:LG07_023000001_023250000 (negative_branch_length)

## Time-tree PCA

- absolute_time: PC1 explains 71.1%, PC2 22.7%; PC1/root-age correlation r=0.991.
- relative_time: PC1 explains 35.3%, PC2 28.0%; PC1/root-age correlation r=0.805.
- relative_feature_zscore_sensitivity: PC1 explains 32.4%, PC2 26.3%; PC1/root-age correlation r=0.846.
- Largest PC1 pair loadings in absolute_time: Gadmor_avc_spc–Gadmor_avo_spc, Gadmor_avc_spc–Gadmor_lfo_spc, Gadmor_avc_spc–Gadmor_ico_spc, Gadmor_avo_spc–Gadmor_low_spc, Gadmor_avo_spc–Gadmor_lfc_spc.
- Largest PC1 pair loadings in relative_time: Gadmor_avo_spc–Gadmor_two_spc, Gadmor_lfo_spc–Gadmor_two_spc, Gadmor_ico_spc–Gadmor_two_spc, Gadmor_bat_spc–Gadmor_icc_spc, Gadmor_icc_spc–Gadmor_twc_spc.
- LG01 mean root age inside minus outside = 0.7492 Ma (cyclic-shift p=0.009091, BH q=0.01449).
- LG02 mean root age inside minus outside = 0.6336 Ma (cyclic-shift p=0.01087, BH q=0.01449).
- LG07 mean root age inside minus outside = 1.091 Ma (cyclic-shift p=0.008475, BH q=0.01449).
- LG12 mean root age inside minus outside = 0.2915 Ma (cyclic-shift p=0.08571, BH q=0.08571).

Per-linkage-group PC1/PC2 contrasts and cyclic-shift nulls are in `pca_inside_outside_spatial_tests.tsv`. The cyclic-shift test rotates the frozen interval mask over the ordered available windows and leaves each observed score track intact. Relative-time PCA removes overall root-age scaling; the separate feature-z-score PCA is only a scaling sensitivity.

## Hierarchical clustering

- euclidean: best silhouette among k=2–7 is k=5 (silhouette 0.581); 80% contiguous-block subsample median ARI=1.000 (2.5–97.5% 0.884–1.000).
- cityblock: best silhouette among k=2–7 is k=3 (silhouette 0.601); 80% contiguous-block subsample median ARI=1.000 (2.5–97.5% 0.970–1.000).

Cluster count/enrichment tables use independent cyclic shifts of the interval mask within each LG, with equal-LG weighting. Do not read cluster membership as a discrete biological history; inspect the representative observed medoid trees and the MRCA features. Cluster enrichment is descriptive and window-local, with only four chromosomes.

Spatial cluster enrichment:
- euclidean cluster 2.0: inside-minus-outside proportion -0.949; cyclic-shift p=5e-05, BH q=0.000133.
- cityblock cluster 1.0: inside-minus-outside proportion -0.949; cyclic-shift p=5e-05, BH q=0.000133.
- cityblock cluster 3.0: inside-minus-outside proportion +0.704; cyclic-shift p=5e-05, BH q=0.000133.
- euclidean cluster 4.0: inside-minus-outside proportion +0.473; cyclic-shift p=0.0003, BH q=0.0006.
- euclidean cluster 1.0: inside-minus-outside proportion +0.244; cyclic-shift p=0.0191, BH q=0.0255.
- cityblock cluster 2.0: inside-minus-outside proportion +0.244; cyclic-shift p=0.0181, BH q=0.0255.
- euclidean cluster 3.0: inside-minus-outside proportion +0.236; cyclic-shift p=0.0682, BH q=0.0779.
- euclidean cluster 5.0: inside-minus-outside proportion -0.005; cyclic-shift p=0.495, BH q=0.495.

## Rooted-triplet analysis

- All 220 cod-population triplets were enumerated. Q contrasts and temporal separations are included for all three possible rooted resolutions; triplets with zero/near-zero internal separation are retained as unresolved.
- Local topology frequencies are frequencies among MCC window point trees, not posterior support for a genome-wide species tree. Temporal separation is reported separately from topology frequency, in Ma and divided by root age.
- Published-tree agreement is computed from the induced rooted topology in `baseline_population_tree_12taxa.nwk`; it compares topology only, not times. The source baseline includes an artificial root-stem value that is irrelevant to the induced topology.
- For all triplets and linkage groups, the summary includes inside/outside topology frequencies, unresolved frequencies, mean/median temporal separation, baseline agreement, effects, and circular-shift null p-values (with BH q-values by metric family).
- Track examples were selected objectively: top two triplets each for maximum LG-specific normalized separation shift, topology-frequency shift, and added baseline disagreement inside. Selected triplets: Gadmor_ico_spc,Gadmor_lfc_spc,Gadmor_lfo_spc, Gadmor_avc_spc,Gadmor_avo_spc,Gadmor_low_spc, Gadmor_bor_spc,Gadmor_twc_spc,Gadmor_two_spc, Gadmor_icc_spc,Gadmor_twc_spc,Gadmor_two_spc, Gadmor_avo_spc,Gadmor_lfc_spc,Gadmor_twc_spc, Gadmor_bat_spc,Gadmor_bor_spc,Gadmor_two_spc.

Per-linkage-group triplet summaries:
- LG01: inside/outside modal rooted topology differs for 120/220 triplets; average normalized separation shift across triplets is +0.290; agreement with the published baseline is 43.4% inside vs 64.9% outside. The largest normalized mean-separation shift is +0.771 for Gadmor_ico_spc,Gadmor_lfc_spc,Gadmor_lfo_spc (p=0.00909, BH q=0.0146).
- LG02: inside/outside modal rooted topology differs for 125/220 triplets; average normalized separation shift across triplets is +0.276; agreement with the published baseline is 42.9% inside vs 62.6% outside. The largest normalized mean-separation shift is +0.643 for Gadmor_avc_spc,Gadmor_avo_spc,Gadmor_lfc_spc (p=0.0109, BH q=0.0146).
- LG07: inside/outside modal rooted topology differs for 106/220 triplets; average normalized separation shift across triplets is +0.382; agreement with the published baseline is 49.6% inside vs 61.8% outside. The largest normalized mean-separation shift is +0.739 for Gadmor_kie_spc,Gadmor_lfo_spc,Gadmor_low_spc (p=0.00847, BH q=0.0146).
- LG12: inside/outside modal rooted topology differs for 130/220 triplets; average normalized separation shift across triplets is +0.156; agreement with the published baseline is 37.7% inside vs 66.8% outside. The largest normalized mean-separation shift is +0.656 for Gadmor_ico_spc,Gadmor_lfc_spc,Gadmor_low_spc (p=0.0857, BH q=0.0898).

## Interpretation

These analyses can establish that the observed local MCC time trees differ in absolute divergence ages, relative MRCA structure, triplet ordering frequencies, or temporal separation across the frozen inversion intervals. They cannot establish that inversions caused those changes, that the MCC trees represent posterior uncertainty, or that a local effect biases a genome-wide divergence-time estimate. Comparisons are descriptive, spatially adjusted by cyclic shifts, and not independent replication across 426 windows.

## Reproduction

`python3 empirical/atlantic_cod/analysis/time_tree_structure/scripts/run_time_tree_structure.py`

Validation tests: `python3 -m unittest discover -s empirical/atlantic_cod/analysis/time_tree_structure/tests -v`.
