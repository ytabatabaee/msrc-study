# Sunflower Stage 0/1 audit

This directory records a reproducible availability and compatibility audit of the sunflower resources of Todesco et al. (2020) and Owens et al. (2023). Raw published tree and coordinate files are kept under `data/sunflowers/raw/` and excluded from Git.

The local tree sources were downloaded only after confirming their public sizes and contents. Rebuild the committed inventories with:

```sh
python3 empirical/sunflowers/scripts/build_stage0_inventories.py
```

The script reads the coordinate-bearing published 10.2-kb tree set and the transcription of the 2020 Extended Data Table 1. Haploblock endpoints in that paper are reported to the nearest Mbp; interval overlap counts therefore use the published rounded bounds and are screening counts, not exact boundary assignments.

See `results/stage0_feasibility_report.md` for the decision gate and limitations. No Stage 1 topology or branch-length comparison has been performed. The published records do not identify the arrangement genotypes of the exact 2023 tree tips at the candidate intervals, which is required to define an arrangement-associated quartet independently.
