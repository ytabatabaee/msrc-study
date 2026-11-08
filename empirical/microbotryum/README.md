# Microbotryum Stage 0 feasibility audit

This directory contains a lightweight feasibility audit for using the Lucotte
et al. (2025) *Microbotryum* mating-type chromosome gene trees in the MSRC
empirical analyses. It is intentionally not a full coordinate-indexed
pipeline.

Run from the repository root:

```bash
python3 empirical/microbotryum/scripts/stage0_audit.py
python3 -m unittest discover -s empirical/microbotryum/tests -p 'test_*.py'
```

The audit requires the Dryad files `gene_trees`, `speciestree.tre`, and
`README.md` in `data/microbotryum/raw/`. The files are now present locally and
the audit parses the deposited trees directly. Large assembly files remain
outside Git.

Relevant provenance:

- Dryad: DOI `10.5061/dryad.0k6djhb89`
- Article: Lucotte et al. (2025), *Nature Communications*, DOI
  `10.1038/s41467-025-60222-5`
- Author code: `https://github.com/eliselucotte/MicrobotryumOnDianthus`
- Code archive inspected: Zenodo DOI `10.5281/zenodo.14918233`
