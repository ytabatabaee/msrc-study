# Stable inheritance probabilities at large gene-copy counts

The inheritance kernel still evaluates
`choose(M_daughter,a) p^a (1-p)^(M_daughter-a)` for each parental count
p=k/M_parent. It now uses scipy.stats.binom.pmf to evaluate that same probability
without converting a huge integer combination to float first. No row
renormalization, approximation to the binomial law, or model change is applied.

New tests check the old direct formula at small equal/unequal sizes, M=2000
finite/nonnegative values, normalization, exact absorbing endpoints, and an
independent log-gamma formula. All 86 moran_theory tests pass, including the
balanced and pectinate quartet tests (23.27 seconds in msrc-sim's environment).

Full M=2000 balanced Moran calculations were retried with 30- and 120-second limits.
The inheritance preflight now passes; the full recursion times out constructing
sister-population operators, using about 306–309 MiB peak RSS at those stops. Thus
the numerical overflow is fixed, but practical full finite-M predictions at
these demographies remain unestablished. No prescribed-path prediction is
presented as a stochastic Moran prediction.

The matched 10-Mb fixed-frequency positive control uses the existing full
multispecies prescribed-frequency operator. Its results are recorded in
msrc-sim's validation/chromosome_biological_validation and
`docs/MSRC_ARG_CHROMOSOME_BIOLOGICAL_VALIDATION.md`. That experiment is conditional
on a specified path and terminal pattern; it is not a posterior Moran sample.
