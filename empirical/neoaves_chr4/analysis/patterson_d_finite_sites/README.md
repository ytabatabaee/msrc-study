# Finite-sites Neoaves Patterson's D prediction

This analysis extends the cached genealogy decomposition by evaluating ABBA
and BABA site-pattern probabilities on the complete induced quartet tree,
including all four pendant branches and the internal branch. It uses JC69,
empirical HKY85, and empirical GTR with four fixed gamma-rate categories.

The original alignment archive and named gene trees are read locally. Previous
Patterson's D, QQS, and decomposition outputs are not modified. Global
substitution parameters are estimated from callable chr4 N61 quartet bases;
they are not fit separately inside and outside the PNAS intervals.

Run from the repository root:

```bash
python3 empirical/neoaves_chr4/analysis/patterson_d_finite_sites/scripts/finite_sites.py
python3 empirical/neoaves_chr4/analysis/patterson_d_finite_sites/scripts/postprocess_finite_sites.py
python3 empirical/neoaves_chr4/analysis/patterson_d_finite_sites/scripts/make_figures_and_report.py
python3 -m unittest discover -s empirical/neoaves_chr4/analysis/patterson_d_finite_sites/tests -p 'test_*.py'
```

The model is conditional on the published gene trees. Since the trees and
site patterns originate from the same alignments, agreement is not an
independent validation of the inferred genealogies.

The postprocessing script reuses the locus-level predictions to calculate
500-kb, 1-Mb, and 2-Mb block-bootstrap intervals and residual diagnostics.
The figure/report script creates the SVG, PNG, and PDF diagnostics. Large
alignment and tree archives remain external inputs.
