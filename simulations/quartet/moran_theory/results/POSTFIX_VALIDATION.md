# Post-fix Moran quartet validation

This report records validation after correcting the authoritative `msrc-sim`
genealogy boundary. The structured genealogy now stops at the explicit
rearrangement origin age and enters an ordinary panmictic ancestral process.
It no longer holds the oldest frequency record fixed indefinitely.

## Tests

- `msrc-sim/tests/test_structured_origin_boundary.py`: 4 passed.
- `msrc-sim/tests/test_events.py` and `test_piecewise_frequency.py`: passed.
- `msrc-study/simulations/quartet/moran_theory/tests`: 15 passed.
- The boundary regression uses `rho=0`, pattern `1010`, finite structured
  duration, and verifies the finite-interval formula rather than a topology
  probability forced to one.

## Shared-population validation

The run used 1,000 independent forward Moran histories per cell with
`M={20,40}`, initial frequencies `{0.25,0.5,0.75}`, effective
background-switching intensities `{0,0.1,0.5}`, and durations
`{0.25,0.5,1.0}` in the existing generation-equivalent units.

| quantity | value |
|---|---:|
| cells | 54 |
| reported rows | 2,670 |
| RMSE | 0.0535295 |
| maximum absolute error | 0.3218391 |
| within 2 MC SE | 94.20% |
| within 3 MC SE | 98.90% |

The maximum error is from a low-frequency conditional pattern row. The
unconditional topology rows have maximum absolute error about 0.0417 and show
no systematic discrepancy or forward/backward ordering bias.

## Full species-tree validation

The run used 1,000 independent forward histories per cell for
`((1,2),(3,4))`, with `M={20,40}`, initial frequencies
`{0.25,0.5,0.75}`, switching intensities `{0,0.1,0.5}`, sister lengths
`{0.25,0.5}`, and root structured lengths `{0.25,0.5}`. The literal simulator
inheritance rule was retained, so zero-duration speciation is distinct from
the collapsed shared-population kernel.

| quantity | value |
|---|---:|
| cells | 72 |
| replicates per cell | 1,000 |
| RMSE | 0.0375635 |
| maximum absolute error | 0.3008072 |
| within 2 MC SE | 94.47% |
| within 3 MC SE | 98.88% |

For the unconditional topology rows, grouped mean signed residuals were
`+0.00088` for `Q1`, `+0.00410` for `Q2`, and `-0.00498` for `Q3`; grouped
RMSE values were `0.01452`, `0.01433`, and `0.01467`, respectively. Grouped
mean residuals were effectively zero for `rho=0`, `0.1`, and `0.5`, and for
both tested root lengths. These results do not show a systematic topology
bias after the boundary correction.

## Joint law and sensitivity

`full_species_tree.py` now returns the exact `16 x 3` joint array
`P(Q_r,Z=z)`, its pattern marginal, and `P(Q_r|Z=z)` from the same operator
recursion. The validation script records marginal, joint, and sufficiently
observed conditional rows. The sensitivity diagnostic evaluates marginal
`q(rho)` and conditional quartet probabilities for patterns `1010`, `0101`,
`0011`, and `0110`. In the tested symmetric setting, marginal changes are
small, while conditional derivatives such as `d P(Q1|Z=1010)/d rho` are much
larger. In the tested `M=20` setting, marginal finite-difference derivatives
ranged from about `-7.17e-4` to `3.58e-4`, while conditional derivatives
ranged from about `-1.00e-1` to `6.21e-2` across the reported patterns and
topologies. This is a numerical sensitivity result, not an identifiability
claim.

The exact model remains the hierarchical simulator model: forward Moran
frequency histories, backward structured genealogy conditional on those
histories, the simulator's binomial speciation inheritance, and ordinary
panmictic ancestry before rearrangement origin. The effective switching
intensity `rho` is not identified with a physical crossover rate here.
