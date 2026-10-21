# Results

The 72-cell discovery grid used 100 independent Moran histories per cell. It
identified no regime with a clearly large time effect and negligible topology
effect after accounting for Monte Carlo noise. The largest discovery values
were therefore carried forward as candidates rather than treated as findings.

Independent validation used four balanced-frequency cells at `M={20,40}`,
`rho={0,0.5}`, sister length `0.25`, and root structured length `0.5`, with
5,000 histories per cell. The strongest candidate was `M=40, rho=0.5`, where
`ΔQ=0.0186` and `ΔT=0.0396`; its approximate Monte Carlo standard errors were
`0.0141` and `0.0196`, respectively. The `M=20, rho=0.5` candidate had
`ΔQ=0.0216` and `ΔT=0.0209`. These are small and do not establish topology/time
signal separation.

A separate long-structured-phase validation used asymmetric initial
frequencies, `rho={0,0.5}`, sister length `0.5`, root length `1.0`, and 1,000
histories per cell. The largest standardized time difference was `0.0685`
(`M=40`, initial frequency `0.75`, `rho=0.5`), with approximate Monte Carlo
error `0.0454`; it was not independently decisive.

The direct finite-M theory comparison remained within the expected Monte Carlo
uncertainty for the validated topology cells. The theory comparison is limited
to `P(Q)` and does not claim to predict the time summaries.

The mutation controls behaved as expected. A common fixed rate doubled pairwise
distances relative to the time scale but left `R_T` unchanged. Window-wide
rate multipliers changed raw pairwise summaries while leaving `R_T` unchanged
to numerical precision. Branch-specific multipliers changed both raw summaries
and `R_T`; within-tree normalization is therefore conditional, not universal.
