# Extended local validation of the pectinate Moran quartet law

## Scope

The exact finite-M recursion for `(((1,2),3),4)` was tested against
independent forward Moran frequency histories and backward structured
genealogies from `/Users/ytabatabaee/Desktop/msrc-sim`. No Bridges-2 job was
submitted, and the prior smoke outputs were retained unchanged.

The extended grid contains 23 deliberately selected, non-Cartesian conditions
and 2,000 independent genealogies per condition (46,000 total). It includes
zero and long branch durations, intermediate-versus-sister branch asymmetry,
unequal chromosome-copy counts, starting frequencies 0, 0.1, 0.25, 0.5, 0.75,
0.9, and 1, all requested switching values, and neutral, positive-selection,
and negative-selection conditions. The complete condition table is
`results/unbalanced_species_tree_extended_conditions.tsv`.

The simulator tree uses an ultrametric pectinate convention. Therefore the
terminal durations supplied to the analytic calculation were `(0, 0, T_S,
T_S + T_I)` for taxa 1--4. This makes the mapping explicit while preserving
the simulator's species-boundary and binomial-inheritance rules.

## Environment and commands

The run used Python 3.9.7 from `/Users/ytabatabaee/opt/anaconda3/bin/python`,
NumPy 1.26.4, SciPy 1.13.1, and pytest 6.2.4. The configured project venv
was checked first but does not contain SciPy; no packages were installed.

```text
/Users/ytabatabaee/opt/anaconda3/bin/python -m pytest simulations/quartet/moran_theory/tests -q
/usr/bin/time -p /Users/ytabatabaee/opt/anaconda3/bin/python \
  simulations/quartet/moran_theory/validate_unbalanced_extended.py \
  --msrc-sim /Users/ytabatabaee/Desktop/msrc-sim --replicates 2000
MPLCONFIGDIR=/tmp/msrc-mpl /Users/ytabatabaee/opt/anaconda3/bin/python \
  simulations/quartet/moran_theory/plot_extended_validation.py
```

The theory test suite passed **23 tests**. The extended validation took 33.49
seconds wall time (31.21 seconds inside the driver) and used 66.57 MiB peak
resident memory. The prior eight-condition smoke validation remains recorded
in `LOCAL_VALIDATION_REPORT.md`.

## Agreement

The extended summary is in `results/unbalanced_species_tree_extended_summary.tsv`;
topology-level rows are in `results/unbalanced_species_tree_extended_validation.tsv`.

| Metric | Value |
|---|---:|
| Conditions | 23 |
| Independent genealogies | 46,000 |
| RMSE across 69 topology comparisons | 0.012183 |
| Maximum absolute error | 0.024384 |
| Topology comparisons within 2 MC SE | 89.86% |
| Topology comparisons within 3 MC SE | 100% |
| Cells with Pearson multinomial p >= 0.05 | 91.30% |
| Maximum probability normalization error | < 3e-16 |

The largest topology-level discrepancy was 0.02438 for q1 in the neutral
`p0_0.5` condition. Its standardized residual was 2.26, compatible with the
targeted Monte Carlo scale. No condition had a topology-level residual beyond
three Monte Carlo standard errors. The cell-level Pearson checks are
descriptive because multiple conditions were examined; they are not a proof
criterion.

Agreement persisted for unequal branch sizes, unequal branch durations,
nonzero switching, and both signs of selection. The extended tests support the
implemented recursion over the tested arbitrary branch-specific finite counts.
They do not prove that the operator is a microscopic ancestral dual of Moran
reproduction.

## Biological behavior in the tested grid

The species-tree topology q1 remained the most probable topology in every
condition. The two discordant probabilities were equal to numerical precision
in the unconditional conditions tested, including asymmetric population-size
cases. No alternative topology exceeded q1, and no unconditional
discordant-topology asymmetry was observed. This is consistent with
exchangeability after averaging over the deterministic ancestral count
distribution used here; it does not imply that arrangement-conditional
distributions are symmetric.

The intermediate branch changes the magnitude of q1 when sister and
intermediate durations are varied. Panel D is descriptive across the targeted
conditions, not a one-factor causal curve, because other attributes vary in
that panel.

## Theorem review and limitations

No correction to `UNBALANCED_SPECIES_TREE_DERIVATION.md` or
`unbalanced_quartet_theorem.tex` was required. The exactness claim remains
restricted to the hierarchical simulator model: finite-state Moran count
histories, binomial daughter initialization conditional on shared parental
counts, structured backward genealogies, and the simulator's finite boundary
rates. The tests now cover unequal branch counts and the simulator-consistent
ultrametric terminal-duration mapping, but not every combination of
branch-specific selection coefficients or large-M numerical conditioning.

The validation figure is
`results/unbalanced_species_tree_extended_validation.pdf` and
`results/unbalanced_species_tree_extended_validation.png`. A larger
Bridges-2 validation is not currently necessary to identify a correctness
problem, but could strengthen a broad parameter-space claim before manuscript
integration.
