# Full363 input and method comparison

| Item | Current reconstruction | Published PNAS source | Assessment |
|---|---|---|---|
| Full gene-tree source | `data/stage4d/named_recollapsed.tre` via current Stage 4E inputs | `63430.aLRT-0.95-collapsed.gene.trees.gz` | Different source collections |
| Full T0 locus count | 63,430 records | 63,430 published collapsed records | Count agrees |
| PNAS exclusion | 1,431 frozen coordinate-selected loci | Published chr4 outlier intervals | Membership/count previously verified; tree identity remains separate |
| T_PNAS retained loci | 61,999 | 61,999 expected after 1,431 exclusions | Count agrees |
| Unresolved source correspondence | 317 records; 200 source-discrepancy classifications | No complete source-side locus resolution in cached snapshot | Provenance incomplete |
| ASTRAL | v1.25.4.8, `-u 2`, `--genelength 1000`, `--sulengthtype 1` | Version/command not present beside cached `score-columbea` trees | Do not attribute RF 2 to version alone |

The current T0 matching the published full-data topology is evidence that the
current taxon labels and the baseline ASTRAL implementation are compatible at
the topology level. It does not validate the reconstructed outlier-removed
gene-tree collection. The RF 2 difference between current T_PNAS and the
published outlier-removed tree is consistent with an input/preprocessing
difference, but the existing 317-record audit does not identify a sufficient
causal subset.

The PNAS repository scripts establish the QQS definition and identify
`score-columbea` as the Figure S9 species-tree directory. They do not provide a
complete reproducible ASTRAL command for the cached species trees in the
versioned snapshot.
