# Pectinate validation record

The compact independent validation was run with the authoritative local
`msrc-sim` checkout and a fresh Moran frequency history plus genealogy for
every replicate. It used 8 parameter cells:

- `M=10`, initial fraction `0.5`;
- `rho` in `{0, 0.17}`;
- sister length in `{0.3}`;
- intermediate length in `{0, 0.2}`;
- root length in `{0, 0.5}`;
- 1,000 independent histories per cell.

Summary:

| metric | value |
|---|---:|
| cells | 8 |
| histories per cell | 1,000 |
| RMSE across topology rows | 0.0127672 |
| maximum absolute error | 0.0299637 |
| fraction within 2 Monte Carlo SE | 1.000 |
| fraction within 3 Monte Carlo SE | 1.000 |

The tabular output is `unbalanced_species_tree_validation.tsv`. This is a
smoke validation, not a publication-scale uncertainty study. The prepared
larger grid is described by
`../cluster/unbalanced_validation_manifest.tsv` and is intentionally not
submitted automatically.
