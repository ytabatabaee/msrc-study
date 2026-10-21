# Moran-integrated quartet validation

This result directory is reserved for the standalone comparison of the sparse
finite-M operator with direct Moran simulations from the authoritative
`msrc-sim` checkout. Each tabular result records the seed, cell parameters,
replicate count, analytic value, Monte Carlo estimate, absolute error, Monte
Carlo standard error, and standardized error.

Run a compact check with:

```text
python simulations/quartet/moran_theory/validate_against_msrc_sim.py --smoke --replicates 1000
```

The publication grid is configurable through `--M`, `--fraction`, `--rho`,
and `--duration`; the default is `M={20,40,80}`, initial fractions
`{0.25,0.5,0.75}`, `rho={0,0.02,0.1,0.5}`, and `T/M={0.25,0.5,1}`.

The retained compact validation run used 36 cells with
`M={20,40,80}`, initial fractions `{0.25,0.5,0.75}`, `rho={0,0.1}`,
`T/M={0.25,0.5}`, 200 replicates per cell, and a minimum of 25 observations
before a conditional pattern row was reported. It gave RMSE `0.03237`, maximum
absolute error `0.20`, 94.8% of reported estimates within 2 MC SE, and 98.1%
within 3 MC SE. The largest residuals were in low-frequency conditional
patterns near the reporting threshold; unconditional topology and pattern
probabilities were consistent with the analytic values at the expected Monte
Carlo precision. This is a validation pilot, not a publication-scale power
run; increase `--replicates` for tighter conditional estimates.
