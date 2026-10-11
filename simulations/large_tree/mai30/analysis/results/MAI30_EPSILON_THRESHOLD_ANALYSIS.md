# Mai30 epsilon threshold analysis

The implemented mixture is

\[
P_\varepsilon=(1-\varepsilon)P_{MSC}+\varepsilon P_{MSRC},
\]

with \(\varepsilon\in\{0,0.1,0.25,0.5,0.75,1\}\). Counts use
`floor(500 epsilon + 0.5)`, so the source counts are exactly
`n_MSRC=500 epsilon` for this grid. Sampling is without replacement within
each realized 500-tree source pool. Three destination mixtures were made for
each epsilon, treatment, and matched complete source pair. The manifest keeps
the source pool, source index, history ID, and destination locus.

Only the six complete MSRC pools were used: all three suppression treatments
for Mai replicates 007 and 042. The incomplete 074 and 005 pools were excluded
from this analysis.

For each sampled informative quartet, the source-pool margin is

\[
\Delta=P(Q_{species})-\max(P(Q_{alt1}),P(Q_{alt2})).
\]

The candidate threshold is computed only when the observed MSC and MSRC
source-pool margins have opposite signs. The machine-readable results are in
`analysis/results/mai30_epsilon_quartets.tsv`. The threshold bootstrap in
`mai30_epsilon_threshold_bootstrap.tsv` resamples the 25 MSRC histories with
replacement and resamples MSC loci conditionally from the realized MSC pool.

The extreme-ILS replicate 007 has many apparent endpoint sign changes in the
finite 500-tree source pools: 42.3–46.7% of the 2,000 sampled quartets, with
source-pool median candidate thresholds around 0.30–0.35. For the first 20
candidate quartets per treatment, history-level bootstrap reversal fractions
are broad (median 0.80–0.875; 5–95% ranges roughly 0.70–1.00), and threshold
intervals span approximately 0.005–0.81. This is evidence of severe finite-
sample and source-history uncertainty in this high-ILS regime. It is not
evidence of a population-level quartet reversal.

The low-ILS replicate 042 has no observed endpoint sign changes in the sampled
quartets for any suppression treatment. Its MSC and MSRC endpoints both favor
the species quartet. No epsilon threshold is reported for that replicate.

These results distinguish three claims. Finite-sample RF degradation can be
measured from the ASTRAL mixtures. A population quartet reversal requires a
sign change after the source distributions are estimated with uncertainty.
Population ASTRAL inconsistency would require showing that the expected full
ASTRAL quartet objective is maximized by an incorrect tree. The present
source-pool analysis establishes neither the second nor the third claim.
The next scientifically justified step is to expand independent MSRC history
pools for the extreme-ILS replicate and test whether candidate sign changes
survive history-level uncertainty.
