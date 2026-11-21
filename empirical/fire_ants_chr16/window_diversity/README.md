# Fire ants window diversity

This species-specific bundle contains the window-based FST, dXY and within-population pi analysis scripts, results and figures. Current values use observed SNP panels; dXY and pi are per observed SNP, not per callable nucleotide.

## Reproduce

Run from the repository root:

```bash
./.venv-malaria310/bin/python empirical/fire_ants_chr16/window_diversity/scripts/analyze_sinvicta_snp_panel.py
./.venv-malaria310/bin/python empirical/fire_ants_chr16/window_diversity/scripts/finalize_sinvicta_tracks.py
```

The scripts read this species' existing genotype/tree metadata sources and write results and figures to the neighboring directories. See `results/` for the methods/results report and `data_audit.md` for sample, source, coordinate and missing-data details.

## Outputs

- `scripts/`: reproducible analysis and plotting scripts.
- `results/`: per-window tables, spatial contrasts, joined tracks, reports and audit notes.
- `figures/`: publication-quality track and inside/outside plots.
