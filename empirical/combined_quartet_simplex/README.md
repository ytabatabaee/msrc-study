# Empirical quartet-support simplexes

This directory builds a standardized, postprocessing-only view of the five empirical QSS datasets used by the existing `qss_dw_aw_combined_all_datasets` figure. The command reads frozen repository tables, harmonizes the topology order to `q1 = baseline/species`, `q2 = arrangement or focal alternative`, and `q3 = third topology`, writes the standardized table and descriptive summaries, and renders PDF, PNG, and SVG simplex figures.

Run from the repository root:

```bash
python empirical/combined_quartet_simplex/build_quartet_simplex.py
```

Outputs are written under `empirical/combined_quartet_simplex/`, including `standardized_qss.tsv`, `inside_outside_summary.tsv`, `data_audit.md`, `provenance.json`, and the figures directory. Tests can be run with:

```bash
pytest empirical/combined_quartet_simplex/tests/test_quartet_simplex.py
```

The plotted observation is one frozen window or binned aggregate. The source tables do not retain a repeated three-frequency vector for each fixed quartet across loci, so an individual-quartet simplex distribution cannot be reconstructed honestly from these inputs. Inside/outside differences are descriptive and do not establish a mechanism or causality.
