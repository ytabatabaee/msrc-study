# Genealogy-based decomposition of Neoaves N61 Patterson's D

This independent analysis asks whether the observed nucleotide ABBA--BABA
contrast can be quantitatively reconstructed from the induced local quartet
topologies and their internal branch lengths.

The analysis reads the existing Stiller alignment archive and named gene-tree
archive locally. It does not download or modify the frozen Patterson's D,
QQS, distance, or integration results. Exact alignment-member identifiers are
used to recover callable sequence exposure per locus.

For each rooted quartet, q2 is `(P1,P3)|(P2,O)` and q3 is
`(P1,O)|(P2,P3)`. Thus q3 contributes to ABBA and q2 contributes to BABA.
The topology prediction uses callable exposure, and the tree prediction uses
exposure multiplied by the induced quartet internal edge.

Run from the repository root:

```bash
python3 empirical/neoaves_chr4/analysis/patterson_d_genealogy_decomposition/scripts/decompose.py
python3 empirical/neoaves_chr4/analysis/patterson_d_genealogy_decomposition/scripts/threshold_sensitivity.py
python3 -m unittest discover -s empirical/neoaves_chr4/analysis/patterson_d_genealogy_decomposition/tests -p 'test_*.py'
```

The finite-sites likelihood correction is intentionally not fit: a valid
correction would require an explicit quartet likelihood and independent model
calibration. The report records this limitation.
