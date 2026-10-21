# MSRC topology–coalescence-time signal separation

## Finding

Under the tested finite-M Moran MSRC model and matched ordinary MSC controls,
the experiment did not find a validated regime with a substantial shift in
true genealogical coalescence summaries and little marginal quartet-topology
separation. The strongest predeclared candidates showed small effects in both
blocks. This is a negative or inconclusive result for the proposed detector
motivation, not evidence that continuous signals can never occur.

The exact quartet operator predicts topology and terminal-pattern laws, not a
general distribution of coalescence times. Direct simulation was therefore
required for every time conclusion. The observed time summaries are
simulation-supported only and are not analytic consequences of the current
finite-M theorem.

## Simulation design

The authoritative `msrc-sim` Moran implementation was imported directly. No
new simulator was introduced. Each MSRC history used the corrected explicit
structured-origin boundary and ordinary panmictic ancestry after that boundary.
The default independent unit was one forward frequency history with one
genealogy; the outputs record both counts. MSC controls used independent
genealogies in the corresponding history slots.

The discovery grid covered `M={20,40}`, initial frequencies `{0.25,0.5,0.75}`,
`rho={0,0.1,0.5}`, sister lengths `{0,0.25}`, and root structured lengths
`{0.25,0.5}`. Validation settings and seeds are recorded in the per-stage
manifests and the aggregate TSV files.

## Observables

For every true genealogy the analysis records all six pairwise MRCA times,
mean pairwise time, `R_T`, topology-specific quartet internal length, total
genealogical branch length, the simulator topology index, and the terminal
background pattern when available. These are genealogical quantities; they are
not SNAPP population-tree times, ML substitution distances, or empirical
`A(w)`/`B(w)` features.

The synthetic rate experiment applies `b_e=mu_w t_e` to the realized true tree.
Common and window-wide multipliers preserve within-tree ratios. Branch-specific
lognormal multipliers do not.

## Implications

The current evidence does not justify presenting branch-length-aware detection
as a direct consequence of the validated MSRC quartet theorem. Continuous
features may still be useful empirically, but they should enter the proposed
multi-interval detector as a separately validated observation block. The next
theory/simulation step would be to derive or validate a time-distribution
operator, or to establish a stronger simulator-supported regime before
claiming complementarity.

Raw outputs are summarized in `results/discovery/`, `results/validation/`,
`results/validation_strong/`, and `results/validation_long/`; no large raw
genealogy files are committed.
