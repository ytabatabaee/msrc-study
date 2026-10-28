# Local image of the conditional finite-M Moran quartet map

The exact conditional map was evaluated with the original count-indexed
Moran solver. The fixed family varies root duration and the terminal duration
of taxon 1 while holding `M=10`, initial count `K=5`, `rho=0.1`, sister
durations `(0.2,0.3)`, and the other terminal durations `(0.1,0.2,0.3)`.

At `(T_root,T_terminal1)=(1,0.2)` and `Z=0101`, finite-difference Jacobians
over step sizes from `2e-3` to `5e-5` had determinant ranges:

| shape | determinant range |
|---|---:|
| balanced | `[-1.29743670e-4, -1.29743643e-4]` |
| pectinate | `[-6.49491037e-4, -6.49490358e-4]` |

The corresponding points are approximately `(q1,q2)=(0.23991861,0.52314100)`
and `(0.13807319,0.72866960)`. These stable nonzero determinants are strong
computational evidence of rank-two local reachability for `Z=0101` in both
shapes.

The inverse-function theorem gives an open neighborhood conditional on a
rigorous nonzero determinant certificate. The current certificate is a
high-precision/step-size stability check, not formal interval arithmetic.
This work therefore does not claim a fully rigorous computer-assisted proof of
the determinant sign.

The combined figure samples the two-duration family for all three paired
patterns in both rooted shapes. It is descriptive; mesh points and apparent
boundaries are not proofs of global image boundaries.

The bounded exact fits found no exact witnesses for the requested noncentral
targets in the restricted `Z=0101`, balanced, `M=4` family. Best residuals
were `0.06346` for `(0.2,0.5,0.3)`, `0.03819` for `(0.05,0.9,0.05)`, and
`0.35736` for `(0.7,0.1,0.2)`. The center is attained at zero root duration.
These failed restricted fits do not disprove coverage by a larger valid family.

Complete interior coverage, coverage by a single pattern, and coverage after
the union over patterns remain unresolved. The result does not imply
unconditional coverage or an unconditional introgression-indistinguishability
theorem.
