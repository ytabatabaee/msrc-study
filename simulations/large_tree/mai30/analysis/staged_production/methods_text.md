# Methods text

We selected ten original no-outgroup Mai D1 species-tree replicates from each
of the three conditions using only the original 500-gene normalized RF
summary.  Within each condition, replicates were sorted by RF and the two
endpoints of each of five equal rank strata were retained.  This fixed,
deterministic rule covers low through high original discordance without using
any MSRC outcome or structural annotation.

Each selected replicate is evaluated under MSRC switching fractions
`1.0`, `0.5`, and `0.1`, with initial frequency 0.5, `c=0.001`, 25
independent forward Wright--Fisher histories, and 20 loci per history.  The
independent history is the primary stochastic unit; the 20 genealogies within
a history are retained together for uncertainty analyses.

The ordinary MSC pilot remains an external control.  ASTRAL-IV is run with
the same settings (`-t 2 -u 1 --length CULength`) for every collection.
Quartet-frequency comparisons, topology RF, and branch-length errors are
reported separately.  Branch lengths are evaluated only on matching
bipartitions using the scaled species tree and `2*N_scaled` conversion.
