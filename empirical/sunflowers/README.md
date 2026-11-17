# Sunflower Stage 0/1 audit

This directory records a reproducible availability and compatibility audit of the sunflower resources of Todesco et al. (2020) and Owens et al. (2023). Raw published tree and coordinate files are kept under `data/sunflowers/raw/` and excluded from Git.

The local tree sources were downloaded only after confirming their public sizes and contents. Rebuild the committed inventories with:

```sh
python3 empirical/sunflowers/scripts/build_stage0_inventories.py
```

The script reads the coordinate-bearing published 10.2-kb tree set and the transcription of the 2020 Extended Data Table 1. Haploblock endpoints in that paper are reported to the nearest Mbp; interval overlap counts therefore use the published rounded bounds and are screening counts, not exact boundary assignments.

See `results/stage0_feasibility_report.md` for the initial decision gate and `results/stage1_analysis_report.md` for the arrangement-blind Stage 1 pilot. Reproduce Stage 1 and its publication-format figures with:

```sh
python3 empirical/sunflowers/scripts/stage1_spatial_analysis.py
python3 -m unittest discover -s empirical/sunflowers/tests -p 'test_stage1_spatial_analysis.py' -v
```

Stage 1 fixes q1 as the independently reported *H. annuus* + *H. argophyllus* pair, then compares quartets with *H. debilis* and *H. paradoxus* as the fourth species. It scores all 48 induced individual combinations per set and averages within each genomic window. It excludes 0.5 Mb at each rounded interval endpoint for the primary interior, keeps genomic gaps visible, and reports branch lengths in substitutions/site. The fixed `A_species_pair_distance_q2_minus_q1` is arrangement-blind and is not an arrangement-defined MSRC A(w). No sample-to-arrangement crosswalk or matched same-chromosome controls could be recovered, so results remain descriptive and the decision is NO-GO for a strong plant MSRC claim. See the arrangement audit for the genotype limitation. Figures are in PNG, PDF, and SVG; raw inputs remain ignored under `data/sunflowers/raw/`.
