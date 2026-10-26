# Local validation of the pectinate Moran quartet recursion

## Environment

Both repositories were present:

- `msrc-study`: `/Users/ytabatabaee/Desktop/msrc-paper`
- `msrc-sim`: `/Users/ytabatabaee/Desktop/msrc-sim`

Validation used the existing SciPy-enabled interpreter
`/Users/ytabatabaee/opt/anaconda3/bin/python`:

- Python 3.9.7
- NumPy 1.26.4
- SciPy 1.13.1
- pytest 6.2.4

The repository's `.venv-malaria` was not changed. It contains NumPy but not
SciPy, and package installation was unavailable in the restricted environment.

## Commands

Theory tests:

```bash
/Users/ytabatabaee/opt/anaconda3/bin/python -m pytest \
  simulations/quartet/moran_theory/tests -v
```

Smoke validation:

```bash
/Users/ytabatabaee/opt/anaconda3/bin/python \
  simulations/quartet/moran_theory/validate_unbalanced_species_tree.py \
  --msrc-sim /Users/ytabatabaee/Desktop/msrc-sim \
  --smoke --replicates 1000
```

The smoke command uses eight cells and creates a fresh Moran frequency history
and genealogy for every replicate. It completed in approximately 6 seconds.

The targeted extension used 300 replicates per cell at `M=10`, `rho=0.17`,
`T_S=0.3`, `T_I=0.2`, `T_A=0.5`, initial fractions `0.25` and `0.75`, and
selection coefficients `+0.25` and `-0.25`. Each selection sign was run as a
separate two-cell command and completed in under 2 seconds.

## Unit-test result

All **20 tests passed**. The tests cover generator row sums and absorbing
states at `K=0`, interior `K`, and `K=M`; inheritance and conservation;
probability decomposition into early and root contributions; the ordinary MSC
limit; zero-duration branches; and taxon relabeling symmetry.

## Eight-cell smoke result

The eight cells used `M=10`, initial fraction `0.5`, `rho` in `{0, 0.17}`,
`T_S=0.3`, `T_I` in `{0, 0.2}`, and `T_A` in `{0, 0.5}`. The generated
topology-level table is `results/unbalanced_species_tree_validation.tsv`.

| metric | value |
|---|---:|
| cells | 8 |
| histories per cell | 1,000 |
| topology comparisons | 24 |
| RMSE | 0.0127672 |
| maximum absolute error | 0.0299637 |
| fraction within 2 Monte Carlo SE | 1.000 |
| fraction within 3 Monte Carlo SE | 1.000 |

The largest residuals were `q1=0.02996` in condition 1 and `q1=0.02924`,
`q3=0.02712` in condition 7. Their standardized residuals were +1.95,
-1.97, and +1.80, respectively. No systematic discrepancy was apparent at
nonzero switching or nonzero intermediate duration.

The figure is
`results/unbalanced_species_tree_validation.pdf` and
`results/unbalanced_species_tree_validation.png`. Panel A shows analytical
versus simulated probabilities with Monte Carlo error bars. Panels B and C
show standardized residuals and maximum absolute errors by condition index.
The exact condition parameters are in the validation TSV.

## Targeted off-center and selection result

The additional four-cell check was not included in the primary figure. It
returned the following independent summaries:

| selection | initial fractions | cells | replicates/cell | RMSE | maximum absolute error | within 2 SE |
|---:|---|---:|---:|---:|---:|---:|
| +0.25 | 0.25, 0.75 | 2 | 300 | 0.0113211 | 0.0211685 | 1.000 |
| -0.25 | 0.25, 0.75 | 2 | 300 | 0.0195921 | 0.0361603 | 1.000 |

These checks support the implementation in off-center frequency and selected
Moran regimes, while remaining too small to establish publication-scale error
bounds.

## Assessment and next steps

The implementation agrees with independent `msrc-sim` simulations for the
tested neutral and selected smoke conditions. This supports the pectinate
recursion as an exact calculation for the stated hierarchical simulator model;
it is not by itself a proof of the mathematical theorem or a validation of
every asymmetric branch-size combination.

No mathematical or implementation bug was found. The broader Bridges-2
validation manifest and array script remain available under
`cluster/`; no cluster jobs were submitted. A useful next step is the prepared
11-cell grid with 2,000 histories per cell, followed by branch-asymmetric
population-size cases if the manuscript requires that generality.
