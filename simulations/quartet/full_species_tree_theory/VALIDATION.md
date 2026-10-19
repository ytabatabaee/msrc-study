# Numerical validation

`validate_against_msrc_sim.py` calls `msrc-sim`'s `simulate_genealogy` with a
fixed piecewise-constant frequency history, so the comparison shares the
simulator's species-boundary and event implementation while evaluating the
theory independently.

The checked asymmetric run used `N=10`, terminal length `0.2`, internal sister
lengths `0.1, 0.3, 0.6` when the full grid is requested, `p=0.25`, and several
switching rates. The script evaluates all 16 terminal patterns and writes one
row per pattern and topology. A focal 2,000-locus run for `Z=1010`, `p=0.25`,
`rho=0.1` gave theory `(0.16230, 0.67738, 0.16033)` and simulation
`(0.15350, 0.68700, 0.15950)`, consistent with Monte Carlo variation.

The generated files are `theory_vs_msrc_sim.csv`,
`theory_vs_msrc_sim.png`, and `convergence_to_msc.png` in the selected output
directory. High switching rates are supported by the script but are slower in
the direct event-driven simulator because the number of switching events per
locus grows with `rho`; the exact matrix evaluator remains fast.
