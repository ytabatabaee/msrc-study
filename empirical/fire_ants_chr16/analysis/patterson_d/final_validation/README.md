# Fire-ant chr16 Patterson's D final validation

Additive finalization of the frozen chromosome-16 nucleotide analysis. It does
not alter original quartet membership, input data, inversion annotations, or
primary results. No cross-dataset analysis, f4, or dXY is included.

## Reproduction

Run from repository root. The sensitivity scan streams the existing verified
VCF once; it does not download data. It evaluates frozen Q1/Q2 with DP 8, 10,
12, 15 crossed with ingroup call fractions 0.50, 0.60, 0.70, 0.80. These are
diagnostics; the original DP≥10/50% estimator remains primary.

```sh
python3 empirical/fire_ants_chr16/analysis/patterson_d/final_validation/scripts/run_sensitivity_scan.py
python3 empirical/fire_ants_chr16/analysis/patterson_d/final_validation/scripts/summarize_branch_lengths.py
MPLCONFIGDIR=/private/tmp/fireant-mpl python3 empirical/fire_ants_chr16/analysis/patterson_d/final_validation/scripts/make_figures.py
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s empirical/fire_ants_chr16/analysis/patterson_d/final_validation/tests -p 'test_*.py' -v
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s empirical/fire_ants_chr16/analysis/patterson_d/tests -p 'test_*.py' -v
```

Original outputs and figure generation remain reproducible with the commands in
the parent Patterson's D README. The local VCF is
`data/fire_ants_chr16/raw/genotypes/gt.vcf.gz`; sample membership, mapping,
windows, trees, and inversion provenance are the frozen inputs referenced by
the additive script and reports.

## Package map

See `results/current_state_audit.md` for prior work and limits;
`quartet_null_audit.md` and `sample_provenance_validation.tsv` for design
validation; `filter_attrition.tsv` and `genotype_filter_sensitivity.tsv` for
filter results; and `reference_bias_limitations.md` for what GT:DP can and
cannot establish. Inversion outputs retain the source's midpoint assignment.
Genealogy results reuse exact induced quartet frequencies and block uncertainty
from the original run; branch-length availability is explicitly reported.

Figures in this directory are additive PNG/PDF/SVG outputs. Sparse tracks use
author-window coordinates, not reconstructed continuous SNP coordinates. No
interpolation is implied across unobserved regions.

The four figure sets are `figure_1_D_tracks`,
`figure_2_D_genealogy_tracks`, `figure_3_filter_sensitivity`, and
`figure_4_sample_gene_reference_diagnostics`, each saved as PNG, PDF, and SVG.
