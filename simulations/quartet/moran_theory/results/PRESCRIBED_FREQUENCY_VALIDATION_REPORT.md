# Prescribed-frequency quartet validation

The new implementation evaluates the balanced and pectinate quartet
recursions conditional on supplied, piecewise-constant background-frequency
histories. It does not replace the finite-M Moran-integrated solver.

## Commands and environment

```text
/Users/ytabatabaee/opt/anaconda3/bin/python -m pytest simulations/quartet/moran_theory/tests/test_prescribed_frequency_quartet.py -q
/Users/ytabatabaee/opt/anaconda3/bin/python -m pytest simulations/quartet/moran_theory/tests -q
/Users/ytabatabaee/opt/anaconda3/bin/python simulations/quartet/moran_theory/validate_prescribed_frequency.py --simulator-replicates 1000
```

The configured interpreter is Python 3.9.7 with NumPy, SciPy, and pytest.
The full theory suite passed 42 tests; the new prescribed-frequency tests
passed 5 tests.

## Numerical checks

The independent segmentwise DOP853 integration agreed with the ordered matrix
exponential propagation for one through four lineages. The largest absolute
difference was `3.75e-13`. A constant-frequency balanced-tree calculation
agreed with the existing fixed-frequency recursion to `3.89e-16`.

A complete balanced-tree calculation with different piecewise histories on
the terminal, sister, and root branches conserved probability exactly to the
reported floating-point precision. The joint terminal-pattern calculation
also passes the marginal identities in the tests, including undefined
conditional probabilities for impossible patterns at fixation.

The validation table is `prescribed_frequency_validation.tsv`; the diagnostic
figure is `prescribed_frequency_validation.pdf` and
`prescribed_frequency_validation.png`.

## Independent simulator smoke check

For a constant-frequency balanced tree, 1,000 independent genealogies from
`/Users/ytabatabaee/Desktop/msrc-sim` were compared with the prescribed
calculation. The observed frequencies were:

| topology | theory | simulation | absolute difference | standardized residual |
|---|---:|---:|---:|---:|
| Q1 | 0.358136 | 0.373 | 0.014864 | 0.98 |
| Q2 | 0.320932 | 0.323 | 0.002068 | 0.14 |
| Q3 | 0.320932 | 0.304 | 0.016932 | -1.15 |

All three deviations are within two binomial standard errors. This is a
conditional constant-history smoke check, with the same branch lengths,
population sizes, switching rates, and terminal sampling convention.

## Scope and remaining validation

The prescribed solver is exact conditional on its supplied trajectories. A
complete stochastic-history validation requires converting jointly generated
Moran histories, including daughter-count initialization and the arrangement
origin boundary, into branch trajectories and averaging the conditional result
after each history. That averaging has not been represented as an
unconditional result here. The existing finite-M count-indexed solver remains
the reference for that law. In particular, the current implementation does
not apply a second inheritance kernel to already conditioned histories.

The appropriate paper formula is therefore

```text
P(Q_i) = E_history[q_early,i(history) + F(history)^T P4(history) h_i].
```

Using one mean trajectory or averaging generators before propagation would be
a different approximation.
