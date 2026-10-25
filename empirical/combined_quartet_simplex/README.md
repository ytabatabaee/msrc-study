# Empirical quartet-support simplexes

This directory builds a standardized, postprocessing-only view of the five empirical QSS datasets used by the existing `qss_dw_aw_combined_all_datasets` figure. The command reads frozen repository tables, harmonizes the topology order to `q1 = baseline/species`, `q2 = arrangement or focal alternative`, and `q3 = third topology`, writes the standardized table and descriptive summaries, and renders PDF, PNG, and SVG simplex figures.

Run from the repository root:

```bash
python empirical/combined_quartet_simplex/build_quartet_simplex.py
```

Outputs are written under `empirical/combined_quartet_simplex/`, including `standardized_qss.tsv`, `inside_outside_summary.tsv`, `overlap_counts.tsv`, `data_audit.md`, `provenance.json`, and the figures directory. Tests can be run with:

```bash
pytest empirical/combined_quartet_simplex/tests/test_quartet_simplex.py
```

The plotted observation is one frozen window or binned aggregate. The source tables do not retain a repeated three-frequency vector for each fixed quartet across loci, so an individual-quartet simplex distribution cannot be reconstructed honestly from these inputs. Inside/outside differences are descriptive and do not establish a mechanism or causality.

Identical q1/q2/q3 coordinates are grouped separately by panel and region class. Each unique coordinate is drawn at its unchanged simplex position with a small, translucent marker whose area follows `min(8, 2.5 + 1.25*log2(n_overlap))` and whose alpha is 0.20; lower-multiplicity groups are drawn first. Inside and outside use blue and orange; boundary uses gray. Circle and square cores also remain distinct when inside and outside observations coincide. The compact size legend shows multiplicities of 1, 10, and 100. No halos, outlines, jitter, or smoothing are used.
