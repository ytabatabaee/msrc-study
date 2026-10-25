# Empirical QSS simplex data audit

All sources are frozen postprocessing tables used by the existing `qss_dw_aw_combined_all_datasets` figure. No gene trees or ASTRAL analyses were rerun.

## Sources
### atlantic_cod
- Source: `data/atlantic_cod/processed/stage4a_window_topology_scores.tsv`
- SHA-256: `3ef09caa4edddd247ad986ff168a4fc86b79d66559293bec51d291ee8cf9328e`
- Topology mapping: q1=baseline/species; q2=arrangement; q3=third
- Region source/convention: `data/atlantic_cod/metadata/region_manifest.tsv`; frozen source inside_inversion/overlaps_inversion flags
- Raw standardized rows: 426
- Retained rows: 426
- Observation unit: one frozen window or 100-kb aggregate, never an individual quartet.

### neoaves_chr4
- Source: `empirical/neoaves_chr4/results/stage10_multiresolution_tracks/binned_quartet_support_100kb.tsv`
- SHA-256: `9cf34736401cdaaa3c19825c42e93f8ebe199b5d83141fa5c425dde3c6e33a17`
- Topology mapping: q1=species history; q2=alternative 1; q3=alternative 2
- Region source/convention: `empirical/neoaves_chr4/results/canonical_structural_regions.tsv`; 100-kb bin midpoint; boundary when a bin crosses an annotated interval without an inside midpoint
- Raw standardized rows: 868
- Retained rows: 868
- Observation unit: one frozen window or 100-kb aggregate, never an individual quartet.

### fire_ants_chr16
- Source: `data/fire_ants_chr16/processed/stage4a_window_quartet_support.tsv`
- SHA-256: `b8ff3ea95ba31ae5553f56f0aeb45a393bdb0ba2b19081fc1aca822dfdd70758`
- Topology mapping: q1=species (S); q2=haplotype/arrangement (H); q3=third
- Region source/convention: `data/fire_ants_chr16/metadata/region_manifest.tsv`; frozen Stage 4A chr16 region labels
- Raw standardized rows: 96
- Retained rows: 96
- Observation unit: one frozen window or 100-kb aggregate, never an individual quartet.

### anopheles_2la
- Source: `empirical/anopheles_2la/results/three_quartet_support/window_three_quartet_support.tsv`
- SHA-256: `ecea36f5d4307d5be6051fc87c2e210866a21d60177a1cd990bf1a74d3e9ec53`
- Topology mapping: q1=species; q2=arrangement; q3=third
- Region source/convention: `source region_class; frozen 2La boundaries are retained in source pipeline`; source region_class
- Raw standardized rows: 988
- Retained rows: 978
- Observation unit: one frozen window or 100-kb aggregate, never an individual quartet.

### house_mouse_t_complex
- Source: `empirical/house_mouse_t_complex/results/stage2_fixed_quartet_scan.tsv`
- SHA-256: `9ad62b579fae816d04b05b234f075e2b6f143198bdb5e05df0df592bddfbcbcd`
- Topology mapping: q1=species; q2=t-complex alternative; q3=other
- Region source/convention: `empirical/house_mouse_t_complex/results/t_complex_inversion_intervals.tsv`; 100-kb bin midpoint; boundary when a bin crosses an annotated interval without an inside midpoint
- Raw standardized rows: 248
- Retained rows: 248
- Observation unit: one frozen window or 100-kb aggregate, never an individual quartet.

## Individual-quartet availability

The retained QSS tables contain window-level or binned aggregates. They do not retain three support frequencies for each fixed quartet across multiple loci. A true individual-quartet simplex distribution therefore cannot be reconstructed without returning to lower-level source data; no such figure is generated.

## Region classes

Inside and outside points use the frozen source labels where present. For interval-based sources, the combined-figure convention is applied by midpoint; intervals crossed by a window whose midpoint lies outside are marked boundary. Boundary and unknown observations remain in the standardized table and are shown separately when present.
