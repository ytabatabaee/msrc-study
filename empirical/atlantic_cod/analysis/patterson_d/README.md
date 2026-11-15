# Atlantic cod nucleotide Patterson's D

This analysis calculates population-frequency ABBA/BABA directly from the deposited Matschiner et al. (2022) VCF genotypes. It does not infer site patterns from local trees. Primary coordinates are the repository-frozen 1-based inclusive `cod_phylogenomics` intervals; matched windows are the existing 250-kb `cod_window_trees.tsv` intervals.

## Reproduce

Raw downloaded VCFs are intentionally excluded from Git under `data/atlantic_cod/raw/zenodo/local_downloads/`. They are the Zenodo record 4560275 files and are validated by `results/download_validation.tsv`. With Python 3 and NumPy, SciPy, Biopython, and Matplotlib installed:

```sh
python3 empirical/atlantic_cod/analysis/patterson_d/scripts/run_analysis.py
python3 -m unittest discover -s empirical/atlantic_cod/analysis/patterson_d/tests -v
python3 empirical/atlantic_cod/analysis/patterson_d/scripts/make_figures.py
```

`run_analysis.py` scans the four inversion SNAPP VCFs and the outside-supergene SNAPP VCF, computes population allele frequencies for frozen 2/2/2-sample ingroup quartets, and uses only homozygous Borsai calls for fixed outgroup polarization. It pools ABBA/BABA contributions before computing regional D. The exact four-taxon tree comparison is unavailable: none of the 426 published SNAPP local trees contains Borsai. The `rooted_triple_contrast_q3_minus_q2` output is explicitly a separate three-population genealogy diagnostic.

Sensitivity runs vary outgroup, ingroup sample size, transitions, and inversion-VCF callability. They do not vary the frozen population triples; geographic/ecotype confounding is a remaining limitation and is listed in the results report.

## Main output

See `results/scientific_interpretation.md`, `results/inside_outside_summary.tsv`, and `results/window_patterson_d.tsv`. Figures A–F are each available in PNG, PDF, and SVG. Raw genotypes remain out of Git. No commit or push was made.
