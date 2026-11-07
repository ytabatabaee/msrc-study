# Heliconius local-tree MSRC pilot

This is an exploratory, staged pilot using the Edelman et al. (2019) Dryad local-tree resources. It is deliberately separate from the validated five-dataset manuscript figures.

Run from the repository root:

```bash
python3 empirical/heliconius/scripts/heliconius_pipeline.py download
python3 empirical/heliconius/scripts/heliconius_pipeline.py all
python3 -m unittest discover -s empirical/heliconius/tests -p 'test_*.py'
```

The download command targets only the small Dryad TSV/ZIP files listed in `data/heliconius/metadata/source_manifest.tsv`; the 11.15-GB HAL alignment is not downloaded. `audit` verifies actual Newick content, coordinates, taxa, and branch-length availability. `normalize` produces `data/heliconius/processed/local_trees.tsv`; `analyze` produces per-tree and 50-kb binned quartet/distance tables; `figures` emits PDF, PNG, and SVG tracks when parsed data exist.

The primary assembly is H. erato demophoon v1. Hmel2.5 data are never merged with it. The quartet is a topology-only exploratory design; no arrangement state is inferred from local topology. 15b/2b coordinates are unresolved until an assembly-validated structural source is added, so no inside/outside claim is made in the initial report.
