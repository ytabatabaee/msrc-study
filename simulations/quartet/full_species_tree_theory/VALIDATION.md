# Numerical validation

The current exact matrix evaluator and the scalar evaluator agree to below
`1e-12` across all 16 patterns for
`m={0,1e-4,1e-2,0.1,1,10,100}` and multiple values of
`lambda,T,T_A`. The pair closed form also agrees with
`scipy.linalg.expm(KT)` to machine precision.

The direct simulator validator is
`validate_against_msrc_sim.py`. It forces terminal patterns, uses the
symmetric simulator setting `p=1/2`, maps `rho=2m`, and uses
`lambda=1/N`. It produces:

- `scalar_vs_matrix_validation.png`;
- `quartet_probabilities_vs_m.png`;
- `conditional_topology_reversal_region.png`;
- `msc_limit_vs_m.png`;
- `theory_vs_msrc_sim.csv`.

Selected 5,000-locus cells with `lambda=0.1`, `T=0.3`, and long shared
root extension gave the following theory versus simulation results:

| Z | m | theory (q1,q2,q3) | simulation (q1,q2,q3) |
|---|---:|---|---|
| 0000 | 0.01 | (0.372043, 0.313978, 0.313978) | (0.3782, 0.3146, 0.3072) |
| 0001 | 0.1 | (0.352853, 0.323573, 0.323573) | (0.3558, 0.3212, 0.3230) |
| 0011 | 0.1 | (0.556667, 0.221667, 0.221667) | (0.5734, 0.2196, 0.2070) |
| 0101 | 0.1 | (0.235922, 0.529858, 0.234220) | (0.2280, 0.5278, 0.2442) |
| 1010 | 0.1 | (0.235922, 0.529858, 0.234220) | (0.2324, 0.5438, 0.2238) |
| 1010 | 10.0 | (0.351394, 0.324303, 0.324303) | (0.3450, 0.3328, 0.3222) |

The earlier 500-locus p-not-equal-to-one-half validation bundle remains
separate evidence for the prior general theory and is not used to validate
these symmetric scalar formulas.
