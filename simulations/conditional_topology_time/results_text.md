# Results

The exact conditional quartet calculation agrees with direct simulation at
the joint-law level. Across the 4-cell validation run, the RMSE of the 48
`P(Q,Z)` absolute errors was 0.00205 and the maximum absolute error was
0.00572. Conditional rows are noisier because they use only the histories in
one pattern stratum: the RMSE of `P(Q|Z)` errors was 0.02865 and the maximum
was 0.08266. Pattern probabilities and joint probabilities are therefore the
appropriate primary theory checks; rare conditional rows are flagged in the
tables.

Conditional time summaries are more separated than the marginal summaries in
some strata, but the effect is not uniform across population sizes or
patterns. In the balanced validation cells, the mean pairwise time for 1010
minus 0011 was 1.90 time units for M=20, rho=0; 1.59 for M=20, rho=0.5;
-0.67 for M=40, rho=0; and 2.66 for M=40, rho=0.5. The corresponding 95%
history-level bootstrap intervals for the two means overlap substantially,
so these contrasts are evidence of a conditional signal worth investigating,
not an independently established large effect. Topology-conditioned rows are
still smaller strata and are treated as exploratory.

The marginal MSRC-versus-MSC total-variation topology differences in the four
validation cells were 0.0044, 0.0036, 0.0220, and 0.0098. Standardized
marginal mean-time differences were 0.0048, 0.0125, 0.0064, and 0.0403.
The exact stratum decomposition reconstructs the observed marginal quartet
probabilities and mean pairwise time to machine precision. This shows how
pattern averaging can dilute a conditional difference: each individual
informative 2:2 pattern has only about 5--8% probability in these cells, and
the marginal is a probability-weighted sum across all 16 patterns.

The result is simulation-supported, not a new time-distribution theorem.
The validated operator predicts the joint and conditional quartet law, while
the continuous-time summaries are calculated from realized genealogies. The
latent pattern `Z` is generally unavailable to an empirical detector, so the
conditional signal cannot be treated as directly observable without an
additional arrangement-inference model.
