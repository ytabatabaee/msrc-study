# Stage 7B theorem-linked coalescent decomposition

Stage 7B expresses the existing Stage-7 coalescent-time signal in the form of the MSRC prediction `E[T_D]-E[T_S] > 0`. Because the same fixed pair sets contribute to every window, pair-specific residualization changes A(w) only by an LG-specific constant. Therefore `Delta_coal = Delta_A` up to rounding; Stage 7B is an interpretive decomposition and robustness analysis, not an independent statistical validation.

| LG | Stage-7 Delta_A | Stage-7B Delta_coal | max difference | circular p | circular BH p | physical p | physical BH p |
|---|---:|---:|---:|---:|---:|---:|---:|
| LG01 | 0.757029649 | 0.757029649 | 3.1e-10 | 0.00917431 | 0.010989 | 0.0217391 | 0.0289855 |
| LG02 | 0.653994231 | 0.653994231 | 1.97e-10 | 0.010989 | 0.010989 | 0.0133333 | 0.0266667 |
| LG07 | 0.435501004 | 0.435501004 | 4.04e-10 | 0.00854701 | 0.010989 | 0.011236 | 0.0266667 |
| LG12 | 0.324552387 | 0.324552386 | 5.06e-10 | 0.00961538 | 0.010989 | 0.0892857 | 0.0892857 |

LG01, LG02, and LG07 are supported under both the circular-shift and physical-coordinate placement nulls. LG12 is supported under the primary circular-shift null but is not significant under the physical-placement sensitivity analysis, so it is less robust to the spatial-null choice.

m_eff is only the descriptive effective long-duration approximation 1/(2 Delta_coal), in inverse published SNAPP time units. It is not a recombination rate and does not identify biological m.

## Arrangement-class pair summaries

| LG | class | n pairs | mean inside residual | median inside residual | IQR | fraction positive |
|---|---|---:|---:|---:|---:|---:|
| LG01 | opposite_arrangement | 32 | 0.777633 | 0.78357 | 0.0250635 | 1 |
| LG01 | same_arrangement | 34 | 0.0208625 | 0.00778923 | 0.0498156 | 0.764706 |
| LG02 | opposite_arrangement | 32 | 0.650137 | 0.655868 | 0.0192036 | 1 |
| LG02 | same_arrangement | 34 | -0.00373785 | -0.00726274 | 0.00916079 | 0.205882 |
| LG07 | opposite_arrangement | 35 | 0.794454 | 1.10478 | 1.1222 | 0.771429 |
| LG07 | same_arrangement | 31 | 0.358949 | 6.08793e-05 | 1.10972 | 0.516129 |
| LG12 | opposite_arrangement | 27 | 0.317486 | 0.323574 | 0.0244218 | 1 |
| LG12 | same_arrangement | 39 | -0.00740564 | -0.00625993 | 0.0231868 | 0.307692 |

LG01, LG02, and LG12 are cleaner arrangement-conditioned examples: all opposite-arrangement pairs have positive inside residuals and their mean residuals greatly exceed the same-arrangement means. LG07 is the strongest topology/time decoupling example, but its pair-level effect is heterogeneous: same-arrangement mean/median inside residuals are about 0.359/0, while opposite-arrangement mean/median values are about 0.794/1.105.

The population jackknife is a dependence and sensitivity analysis, not an iid argument. All leave-one-population-out estimates remain positive. LG01, LG02, and LG12 are tightly stable; LG07 is substantially more variable.

| LG | Delta_D | topology p | Delta_coal | coalescent p | BH p | class |
|---|---:|---:|---:|---:|---:|---|
| LG01 | 1.23277 | 0.009174 | 0.75703 | 0.00917431 | 0.010989 | topology_and_time_shift |
| LG02 | 0.892332 | 0.010989 | 0.653994 | 0.010989 | 0.010989 | topology_and_time_shift |
| LG07 | -0.005981 | 0.555556 | 0.435501 | 0.00854701 | 0.010989 | time_shift_without_topology_shift |
| LG12 | 1.04317 | 0.038462 | 0.324552 | 0.00961538 | 0.010989 | topology_and_time_shift |

The finite-duration theory predicts a positive arrangement-conditioned delay, while the infinite-duration symmetric limit is `E[T_D]-E[T_S]=1/(2m)`. Atlantic cod shows a theory-consistent empirical signature; it does not establish MSRC as the unique historical cause or separately identify m, lambda, and tau. The analysis uses published per-window SNAPP MCC point estimates because posterior samples are unavailable.
