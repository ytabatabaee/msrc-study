# Exact-population local gene-tree quartet comparison

This extension scores the 4,046 archived Kelemen & Vicoso chr17 5-kb local trees using the exact frozen sample membership in `../config/frozen_quartet_members.tsv`. For each tree and each design, vectorized four-point distances classify each induced quartet once as q1=(P1,P2)|(P3,O), q2=(P1,P3)|(P2,O), or q3=(P2,P3)|(P1,O). Tied four-point sums are unresolved. Frequencies are within-tree frequencies over resolved induced quartets and sum to one; trees are the spatial inferential unit. `D_topology=q3-q2` corresponds to nucleotide ABBA−BABA direction.

The tree window names are offset by 5 Mb in the source archive. The verified local-tree coverage is chr17:5,000,000–37,055,000 (mm10/GRCm38, half-open internal windows); no tree comparison is extended beyond it. Carrier individuals are represented by the source paper's pseudo-t-haplotype sequence tips; the VCF scan instead uses unphased diploid genotypes.

## Run

```sh
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/exact_quartet_comparison/scripts/run_analysis.py
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/exact_quartet_comparison/scripts/match_and_correlate.py
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/exact_quartet_comparison/scripts/summarize_and_plot.py
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/exact_quartet_comparison/tests/test_quartet_scoring.py
```

Run the polarization validation first. Main tables are under `results/`; six figure families in PNG, SVG, and PDF are under `figures/`. Aggregate q statistics weight each gene tree equally; induced quartet totals are diagnostics and are never treated as independent loci. The spatial-correlation table includes informative-SNP thresholds, resolved-quartet coverage thresholds, 500-kb block-bootstrap intervals, and paired circular-shift spatial-null results.
