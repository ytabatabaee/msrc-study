# House mouse window diversity

This species-specific bundle contains the window-based FST, dXY and within-population pi analysis scripts, results and figures. Current values use observed SNP panels; dXY and pi are per observed SNP, not per callable nucleotide.

## Reproduce

Run from the repository root:

```bash
./.venv-malaria310/bin/python empirical/house_mouse_t_complex/window_diversity/scripts/analyze_t_pseudohap_snp_panel.py
./.venv-malaria310/bin/python empirical/house_mouse_t_complex/window_diversity/scripts/analyze_full_chr17_snp_panel.py
./.venv-malaria310/bin/python empirical/house_mouse_t_complex/window_diversity/scripts/plot_chr17_inversion_barplots.py
```

The scripts read this species' existing genotype/tree metadata sources and write results and figures to the neighboring directories. See `results/` for the methods/results report and `data_audit.md` for sample, source, coordinate and missing-data details.

## Outputs

- `scripts/`: reproducible analysis and plotting scripts.
- `results/`: per-window tables, spatial contrasts, joined tracks, reports and audit notes.
- `figures/`: publication-quality track and inside/outside plots.
