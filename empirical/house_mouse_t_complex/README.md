# House mouse t-complex empirical analysis

This directory stages an empirical MSRC analysis of the house-mouse t-haplotype / t-complex on chromosome 17 using the published Kelemen & Vicoso local-gene-tree archive.

Stages 0 and 1 audit and freeze the published trees. Stage 2 adds exact fixed-quartet scoring, balanced-sampling ASTRAL controls, seven-population sensitivity, filtering robustness, spatial thinning, and a validated main figure. Stage 2B adds split-specific ASTRAL annotation correction, arrangement-state quartet decomposition, 1,000-replicate exact balanced resampling, and the complete spatial phase sweep.

Run Stage 0 with a local copy of the upstream archive:

```bash
python empirical/house_mouse_t_complex/scripts/00_audit_house_mouse_data.py \
  --archive /path/to/IST-2017-78-v1+1_Data.zip
```

The primary dataset is:

```text
Data/2-Coverage-and_AlleleRatio-Filtered_RAW_SNPs/
  1-Trees_for_all_5kb_windows/
    ML_trees.zip
```

The PASS-SNP and coverage-filtered RAW-SNP ML archives are inventoried only as possible robustness datasets.

Important limitations are enforced in the scripts:

- The 2017 chr17 window trees are not treated as a genome-wide or outside-inversion background.
- No tip-to-subspecies, population, or t-status mapping is inferred without source-backed evidence.
- Tentative mappings are excluded automatically from ASTRAL mapping files.
- Stage 1 stops unless the Stage 0 mapping gate passes.

Synthetic tests do not require the 144 MB source archive:

```bash
python empirical/house_mouse_t_complex/scripts/00_audit_house_mouse_data.py --run-tests
pytest empirical/house_mouse_t_complex/tests
```

Stage 2B outputs are in `results/stage2b_report.md`. The visualization and
recombination extension is generated with:

```bash
python empirical/house_mouse_t_complex/scripts/02l_import_recombination_classification.py \
  --archive /path/to/IST-2017-78-v1+1_Data.zip
python empirical/house_mouse_t_complex/scripts/02n_summarize_q_heterogeneity.py
python empirical/house_mouse_t_complex/scripts/02m_make_visualizations.py
```

It produces four-group and seven-population ASTRAL drawings, direct quartet
tracks, the source-defined phylogeny-based recombination-state track, the
aligned q/recombination figure, controls, and
`figures/house_mouse_t_complex_main_v4.{png,pdf}`. The recombination track is
an inferred local-tree classification, not a direct cM/Mb rate estimate. The
Stage 2 scripts use repository-relative output paths and retain source archive
provenance by basename, SHA256, and byte size.

The q/recombination figure also includes a published inversion-boundary panel.
Its four historical intervals use sequence-validated GRCm38/mm10 coordinates
recorded in `results/t_complex_inversion_intervals.tsv`, with supporting
evidence and the fifth overlapping inversion documented in
`results/t_complex_inversion_provenance.md`. The archive and structural study
use the same mm10/GRCm38 coordinate system. The descriptive inside/outside
summary is `results/t_complex_inversion_interval_summary.tsv`.

The independent recombination-suppression audit is in
`results/recombination_suppression_source_audit.tsv` and
`results/recombination_suppression_report.md`. The direct literature evidence
is arrangement-specific but only marker-pair or region-wide: it does not
provide a coordinate-resolved +/+ versus +/t rate track for the frozen 5-kb
windows. The separate four-panel schematic
`figures/house_mouse_t_complex_suppression.{png,pdf}` shows the published
inversion intervals, interval-level cM observations, and the frozen quartet
support and contrast profiles. The source-defined phylogeny-based
recombination-history classification remains supplementary and is not treated
as a quantitative recombination rate.
