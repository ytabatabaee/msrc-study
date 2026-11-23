# Stochastic Moran chromosome validation (2026-11-23)

The companion msrc-sim experiment now uses genuine shared finite-copy-count
Moran histories, binomial speciation inheritance, and exact pre-genealogy
terminal-pattern rejection. Its code and results are in
`validation/chromosome_moran_validation/` and
`docs/MSRC_ARG_10MB_MORAN_VALIDATION.md` in that repository.

The independent M=20 balanced control has nonzero terminal/sister/root branches
and compares 4000 independent ARG chromosomes with the original full stochastic
multispecies recursion. Maximum unconditional q error is 0.01022; all 15
unconditional q criteria meet the practical +/-0.05 interval criterion. There
are 100 naturally accepted 1010 patterns (2.5%, exact probability 2.53302%).
Conditional maximum q error is 0.05979: three criteria meet equivalence, twelve
remain unresolved by precision; none demonstrates error beyond +/-0.05.
No history or genealogy is discarded based on survival or success.

At M=2000, a fixed 64-prior-draw audit found one 1010 pattern, with exact
95% probability interval [0.0003955,0.0840103]. A separately documented
sequential rejection extension accepted its first additional draw (seed64).
That extension is not pooled into a fixed-n acceptance estimate. Conditional
ARG failures are retained, not replaced by another accepted chromosome.

A new separately profiled large-M exact attempt still times out. Stable
binomial inheritance remains fixed; the new profile is dominated by
expm_multiply/sparse count propagation before the full recursion completes.
No prescribed, single-population or diffusion prediction replaces the missing
large-M target. No new theory formula or operator change is made here.

The 10-Mb runs use root age 200, structured root 1000, initial K=1000 of M=2000,
marker 5 Mb, suppressed interval [3,7)Mb, r=1e-8, kappa=.01 and s=0. Repeated
scheduling across dense frequency histories is expensive. Finished and censored
pilot outcomes, exact seeds, resource limits and the separate T200 sensitivity
are archived in msrc-sim. Analysis has 40/200/400/1000 span-weighted windows on
the same native trees. Completed subsets of censored cohorts are labeled
examples, not unbiased unconditional/posterior expectations.

The intended 100 independent chromosomes per main condition and exact large-M
marginal validation remain outstanding. Next priorities are exact-preserving
ARG scheduling, faster finite-M propagation, and adequate independent replicate
counts. Small-M theory agreement does not certify the 10-Mb large-M joint law.
