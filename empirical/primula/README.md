# Primula S-locus feasibility audit

This is a staged MSRC evaluation of Stubbs et al. (2023), not a presumption that the Primula S-locus generates a spatial genealogy signal. The Dryad landing page was verified on 2026-11-16, including the tree, species-tree, S-locus-tree, README, and 11 chromosome partition filenames. File payload retrieval from the published Dryad download endpoints returned HTTP 403 in this environment, so no tree data have been analyzed and no coordinates or statistics are reported yet.

The published methods state that the authors split the 11-chromosome *P. veris* assembly into nonoverlapping 50-kb fragments and selected every other fragment, yielding 3,623 alignments and gene trees. Adjacent sampled trees therefore represent at least 100 kb start-to-start; intervening fragments must remain unobserved. Tree order alone is not accepted as a coordinate map.

## Reproduction

After placing the verified Dryad files in `data/primula/raw/dryad_2022/`, run:

```sh
python3 empirical/primula/scripts/stage0_audit.py
python3 -m unittest discover -s empirical/primula/tests -p 'test_*.py' -v
```

The raw data directory is excluded from Git. The script inventories Newick records, taxon sets, branch-length availability, malformed records, and partition-file syntax. It intentionally refuses to infer chromosome coordinates unless an explicit, verified tree-to-partition mapping is available.

Stage 1 is not run until Stage 0 verifies tree-coordinate mapping and local sampling around an independently annotated S-locus interval. The S-locus topology's published discordance from nuclear history is not by itself evidence of a local spatial genealogy shift.

## Sources

- Stubbs et al. (2023), [New Phytologist](https://doi.org/10.1111/nph.18525).
- [Dryad dataset](https://doi.org/10.5061/dryad.b2rbnzsgz), version published 2022-10-13.
- Reference assembly cited in the paper: chromosome-scale haploid *P. veris* assembly of Potente et al. (2022).
