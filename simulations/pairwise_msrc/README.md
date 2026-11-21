# Pairwise MSRC theory-to-data analysis

This package derives expected pairwise coalescence times under the legacy
single-locus MSRC structured-coalescent model, maps them to infinite-sites
diversity statistics, validates the numerical operator against the existing
`msrc-sim` genealogy generator, and compares the predictions with existing
callable-base Anopheles 2La tracks.

## Reproduce

From the repository root:

```sh
.venv-malaria310/bin/python -m unittest discover -s simulations/pairwise_msrc/tests -p 'test_*.py' -v
.venv-malaria310/bin/python simulations/pairwise_msrc/scripts/validate_msrc_generator.py --histories 100
.venv-malaria310/bin/python simulations/pairwise_msrc/scripts/anopheles_comparison.py
```

The validation requires the sibling, read-only `msrc-sim` checkout at the
recorded commit. No external data are downloaded. The analysis writes under
this package's `results/` and `figures/`.

## Contents

- `derivation.md`: equations, assumptions, simulator audit, finite-origin
  treatment, validation, empirical caveats, and identifiability.
- `scripts/pairwise_msrc.py`: reusable exact piecewise-constant expectation
  calculator and infinite-sites moment transforms.
- `scripts/validate_msrc_generator.py`: independent frequency-history and
  direct event-time validation using the existing genealogy API.
- `scripts/anopheles_comparison.py`: matched-callable-site tracks, block
  summaries, joins to frozen genealogy windows, parameter-sensitivity outputs,
  a labeled equal-weight moment-distance grid comparison, and publication
  figures.
- `tests/`: focused numerical and simulator consistency tests.
- `results/`, `figures/`: generated data summaries and PDF/PNG outputs.

The grid ranking is a descriptive identifiability diagnostic, not an
inferential fit: it has no likelihood, formal parameter uncertainty, or
geographic/introgression null. Its best points are not published as MSRC
parameter estimates.

Anopheles relationships to `A_distance`, mean pair distance, tree diameter,
and quartet tracks are exploratory spatial associations. The NJ distances are
substitution units, and the arrangement-stratified samples do not estimate
natural arrangement frequency. The other species' current diversity outputs
are SNP-conditioned and are not used as absolute per-base estimates.
