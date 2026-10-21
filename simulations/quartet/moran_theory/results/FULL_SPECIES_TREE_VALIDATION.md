# Full species-tree validation

`validate_full_species_tree.py` compares one genealogy from each independent
forward Moran history against the exact operator. It uses
`msrcsim.moran.simulate_frequency_history` and
`msrcsim.structured_coalescent.simulate_genealogy` directly. Loci sharing a
frequency history are not treated as independent replicates.

The script also writes `rho_sensitivity.tsv`, a numerical diagnostic of
`q(rho)` and a finite-difference derivative over a configurable grid. This is
an exploratory sensitivity analysis, not an identifiability proof.

The retained 72-cell pilot used `M={20,40}`, initial fractions
`{0.25,0.5,0.75}`, `rho={0,0.1,0.5}`, sister lengths `{0.25,0.5}`, root
lengths `{0.25,0.5}`, and 200 independent histories per cell. It gave RMSE
`0.03031`, maximum absolute error `0.09273`, 94.1% within 2 MC SE, and 99.3%
within 3 MC SE. A focused two-cell 300-history check gave RMSE `0.01411`,
maximum absolute error `0.02899`, with every reported quantity within 2 MC SE.

The sensitivity file used `M=20`, initial frequency `0.5`, sister length
`0.25`, and root length `0.25`. Across `rho=0,...,1`, `q1` changed from
approximately `0.35057` to `0.35041`; `q2=q3` changed correspondingly from
approximately `0.32471` to `0.32479`. This is a weak but nonzero numerical
signal over this setting, not evidence of injective or statistically
identifiable recombination intensity.

As a separate fast-switching exploration at `M=20`, `p0=0.5`, sister lengths
`0.25`, and root length `0.25`, the analytic values moved from
`(0.35058,0.32471,0.32471)` at `rho=0` to
`(0.34980,0.32510,0.32510)` at `rho=1000`. The matched finite-root ordinary
MSC calculation using pair rate `1/(2Ne)` is approximately
`(0.34979,0.32510,0.32510)`. This supports numerical convergence in this
setting; it is not a general limit proof.
