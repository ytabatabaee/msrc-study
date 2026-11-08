# Heliconius local-tree MSRC pilot

This is an exploratory, staged pilot using the Edelman et al. (2019) Dryad local-tree resources. It is deliberately separate from the validated five-dataset manuscript figures.

Run from the repository root:

```bash
python3 empirical/heliconius/scripts/heliconius_pipeline.py download
python3 empirical/heliconius/scripts/heliconius_pipeline.py all
python3 -m unittest discover -s empirical/heliconius/tests -p 'test_*.py'
```

The download command targets only the small Dryad TSV/ZIP files listed in `data/heliconius/metadata/source_manifest.tsv`; the 11.15-GB HAL alignment is not downloaded. `audit` verifies actual Newick content, coordinates, taxa, and branch-length availability. `normalize` produces `data/heliconius/processed/local_trees.tsv`; `analyze` produces per-tree and 50-kb binned quartet/distance tables; `figures` emits PDF, PNG, and SVG tracks when parsed data exist.

The primary assembly is H. erato demophoon v1. Hmel2.5 data are never merged with it. The analysis enumerates all four-taxon subsets of each source-specific seven-taxon universe (35 quartets for complete trees), freezes a per-quartet baseline from a genome-wide background excluding chromosomes 2 and 15, and assigns the two alternatives deterministically. Per-quartet assignments and aggregate summaries are retained separately. No arrangement state is inferred from local topology. H. erato 2b and Hmel2.5 2b/15b are independently documented; H. erato 15b remains unresolved in the available assembly mapping.
