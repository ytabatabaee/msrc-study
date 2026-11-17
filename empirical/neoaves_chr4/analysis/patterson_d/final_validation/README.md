# Neoaves chromosome-4 nucleotide D final validation

This additive package validates existing Stiller intergenic nucleotide Patterson’s D outputs, tests the D=0 rooted-quartet assumption against published species-history evidence, compares exact matched quartet genealogies, and produces manuscript figures. Frozen data, primary results, and original figures remain untouched.

## Main result

The N61 inside/outside contrast is stable across the saved alternate galliform outgroup, transitions-only and transversions-only counts, separate PNAS intervals, interval exclusions, and modest source-coordinate boundary perturbations. Columbea has a large positive D both inside and outside but no pooled regional difference; its frozen quartet is disputed by the published Stiller genome-wide species tree versus broader Jarvis Columbea evidence. N62 has a small nonzero baseline but no supported all-sites regional shift. The N61 result is site-pattern/topology-concordant on exact matched trees; this does not uniquely identify introgression or MSRC.

See `results/final_scientific_interpretation.md`, `results/species_tree_null_audit.md`, and `results/region_coordinate_validation.md` before using any estimate.

## Reproduction

From the repository root, using Python with Biopython, NumPy and Matplotlib available:

```bash
python empirical/neoaves_chr4/analysis/patterson_d/final_validation/scripts/exact_quartet_comparison.py
python empirical/neoaves_chr4/analysis/patterson_d/final_validation/scripts/augment_species_tree_audit.py
python empirical/neoaves_chr4/analysis/patterson_d/final_validation/scripts/prepare_validation.py
python empirical/neoaves_chr4/analysis/patterson_d/final_validation/scripts/spatial_block_bootstrap.py
MPLCONFIGDIR=/private/tmp/msrc-mpl python empirical/neoaves_chr4/analysis/patterson_d/final_validation/scripts/make_figures.py
python -m unittest discover -s empirical/neoaves_chr4/analysis/patterson_d/final_validation/tests -p 'test_*.py' -v
```

The exact-quartet command scans the existing compressed gene-tree archive (305 MB), not the 3.6-GB nucleotide alignment archive. The output preparation scripts reuse saved locus-level nucleotide counts. The alignment archive scan was not rerun.

## Package contents

- `results/`: current-state, null-topology, estimator, coordinate, region, spatial, exact genealogy, interpretation, methods, captions and checklist files.
- `scripts/`: additive generators for exact quartet matching, species-history annotation, result tables, block bootstrap, and figures.
- `tests/`: independent scalar ABBA/BABA cases and integration checks.
- `figures/`: Figure A–D as PNG, PDF and SVG.

Final assessment: **CONDITIONAL GO** for reporting the N61 regional site-pattern contrast and qualified descriptive results. The quartet null, representative sampling, finite-sites process, assembly crosswalk, and mechanism remain limitations.
