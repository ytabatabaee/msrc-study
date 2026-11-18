# Neurospora L10 Patterson's D final validation

This is an additive validation package for the frozen analysis in the parent
directory. It retains the original quartet, estimator, coordinates and primary
results. No input data are downloaded. Scripts read the existing raw archive
and alignment from the repository's ignored `data/neurospora/raw/` directory.

## Reproduction

Run from the repository root. The full validation scan takes several minutes.

```sh
python3 empirical/neurospora/analysis/patterson_d/final_validation/scripts/run_final_validation.py
python3 empirical/neurospora/analysis/patterson_d/final_validation/scripts/summarize_final_validation.py
MPLCONFIGDIR=/private/tmp/neurospora-mpl python3 empirical/neurospora/analysis/patterson_d/final_validation/scripts/make_figures.py
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s empirical/neurospora/analysis/patterson_d/final_validation/tests -p 'test_*.py' -v
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s empirical/neurospora/analysis/patterson_d/tests -p 'test_*.py' -v
```

`run_final_validation.py --samples-only` refreshes the isolate-sensitivity
analysis without repeating the raw genotype scan. The summary and figure scripts
read the generated result tables.

## Results guide

- `results/final_scientific_interpretation.md`: decision and qualifications.
- `results/current_state_audit.md`, `quartet_null_audit.md`: frozen inputs,
  topology and null-model scope.
- `results/estimator_validation.tsv`: independent scalar estimator comparison.
- `results/sample_provenance_audit.tsv`, `sample_sensitivity_report.md`,
  `sample_sensitivity.tsv`: isolate provenance and sensitivity.
- `results/ascertainment_audit.md`, `ascertainment_diagnostics.tsv`,
  `ascertainment_common_support.tsv`: SNP ascertainment, callability and
  window-matching checks.
- `results/spatial_extent_interpretation.md`, `spatial_extent_analysis.tsv`,
  `spatial_profile_bins.tsv`, `spatial_null_audit.md`,
  `spatial_uncertainty_validation.tsv`: spatial extent and uncertainty.
- `results/genealogy_interpretation.md`, `genealogy_correspondence.tsv`:
  nucleotide D and source-tree topology-weight correspondence.
- `results/unassigned_219_audit.tsv`: source rows outside the fixed 50-SNP
  tree windows.
- `results/final_results.tsv`, `validation_checklist.tsv`,
  `manuscript_methods.md`, `figure_captions.md`: compact outputs for review
  and manuscript use.
- `figures/`: four figure sets in PNG, PDF and SVG.

The primary D estimates are conditional on the filtered variable-site matrix.
TWISST q3 is not an exact rooted D topology because the deposited trees do not
include the *N. discreta* outgroup. The circular-shift result is descriptive;
its wraparound and stationarity assumptions do not justify strong significance
language. See the interpretation reports for the full limits.
