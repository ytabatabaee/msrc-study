# Methods

This experiment asks whether latent terminal arrangement patterns expose
genealogical-time differences that are diluted in marginal quartet summaries.
It uses the authoritative `msrc-sim` implementation, commit
`1658db065a3b42a2fc1d0bebda874eecc4da6d3b`, and the finite-
population quartet operator in `simulations/quartet/moran_theory/`.

Each cell uses the four-taxon tree `((1,2),(3,4))`, one genealogy per
independent forward Moran frequency history, and the same history as the
uncertainty unit. Terminal arrangements are sampled from the simulator and
recorded as `Z` in taxon order 1,2,3,4. No complementary patterns are pooled.
The ordinary MSC is simulated separately as an unconditional comparator;
because MSC has no arrangement variable, MSC quantities conditional on `Z`
are undefined.

For every MSRC genealogy we retain the true rooted Newick tree and calculate
all six pairwise MRCA times, their mean, the realized quartet topology, the
topology-specific internal branch, total tree length, `R_T`, and root
coalescence time. These are true simulator time units, not sequence-derived
branch lengths. Pattern-stratified means, quantiles, and 95% history-level
bootstrap intervals are reported. The topology-conditioned table stratifies
the same quantities by both `Z` and `Q`; rows with fewer than 20 histories are
flagged and are not treated as reliable conditional summaries.

The exact theory comparison evaluates `P(Q,Z)`, `P(Z)`, and `P(Q|Z)` using
the full-species-tree finite-M operator. It does not claim an analytic
distribution for genealogical times. Marginal identities are checked by
reconstructing `P(Q)` and `E[T]` from the natural pattern strata.

The discovery run has 12 cells, 500 histories per cell, and seed 20261121:
`M={20,40}`, initial frequency `{0.25,0.5,0.75}`, `rho={0,0.5}`, sister
length 0.25, and root structured length 0.5. The independent validation run
has four balanced-frequency cells, 5,000 histories per cell, and seed
20261131. The preceding 72-cell topology/time experiment is retained as a
separate committed result and was not silently rerun here.
