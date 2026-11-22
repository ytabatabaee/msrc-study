# Exact multispecies Moran theory versus spatial MSRC simulators

The reusable workflow is maintained in `msrc-sim/validation/moran_spatial`.
It imports the existing theory modules directly, without copying or changing
the quartet theorem. The full report is
`msrc-sim/docs/MSRC_MORAN_THEORY_ARG_SMC_VALIDATION.md`; the parameter audit is
`msrc-sim/docs/MSRC_MORAN_SPATIAL_COMPATIBILITY.md`.

Theory source: this repository commit `ba260f8d4303b5d3bc884fb76a9c5a5032b7e477`.
Simulator validation commit: `5138370aa455eedc54b694c13f5aba91c9de5fa6`.
The simulator commit exists locally; its GitHub push failed because HTTPS
credentials were unavailable. The committed workflow/results can be inspected
in that local checkout.

## Completed comparison

- 13 primary settings, 4,000 independent chromosomes per method per setting:
  full balanced and pectinate recursions, fixed frequencies, fixed realized
  Moran paths and unconditional absorbing Moran histories. Nonzero terminal,
  sister/intermediate and root intervals, unequal M, selection and mixed
  root laws are included. No production simulator or theorem was changed.
- ARG: 249/252 primary practical CI passes, three unresolved, no demonstrated
  error beyond the ±0.05 tolerance. Some fixed-realized-path D intervals
  exclude zero but overlap the practical boundary; small bias is unresolved.
- SMC′: 251/252 frozen-target passes, including 51/52 reference criteria.
  Against position-dependent theory: 198 pass, 36 demonstrated mismatches,
  18 unresolved. Whole-chromosome marginal exactness is FAIL.
- Fixed-attempt terminal conditioning: 320,000 prior draws yield 10,364
  selected Z=1010/0101 chromosomes. All 16 ARG and 16 SMC′ reference checks
  pass; no selected simulator failures. The single-locus control passes
  15/16, with one unresolved D interval.
- An additional 4,000-chromosome pectinate audit passes all nine checks of
  topology mass accumulated before all four lineages reach the root.
- Full simulator tests: 335 passed, one existing warning. Adapter tests:
  19 passed. Selected original theory tests: 34 passed. Dense/sparse and
  matrix-free numerical checks agree near machine precision. Configurations,
  source hashes, probability tables, failures and PNG/SVG profiles accompany
  the full report.

## Boundary and conditioning decisions

The theorem's finite floor and the spatial strict-boundary rules differ on
the global formal state space. For the existing no-mutation absorbing Moran
process with unbiased binomial inheritance, absent backgrounds cannot
reappear in younger descendants; supported backward genealogies cannot enter
them. Every positive occupied count is at least one, so the floor is inactive.
This identifies a matched reachable domain without filtering survival or
replacing failed histories. Replayed validation histories include actual
within-branch fluctuations, losses and fixations.

This conclusion is restricted to the audited root-origin hierarchical model.
It is not a claim about arbitrary supplied paths, mutation/resurrection,
sub-one-copy frequencies, alternative founder conventions, microscopic Moran
ancestral duality, or the spatial joint law. Terminal-pattern conditioning uses
the correct posterior history given Z, not fixed tips on unconditioned prior
histories or equally weighted conditional quartets.

## Reproduction

From the simulator checkout, point the workflow at this theory checkout:

```sh
export MSRC_THEORY_SOURCE=/path/to/msrc-paper
.venv/bin/python validation/moran_spatial/run.py \
  --theory-source "$MSRC_THEORY_SOURCE" \
  --config validation/moran_spatial/cases.json \
  --output /tmp/moran-spatial --replicates 4000
```

The workflow README contains additional-setting, explicit conditional-Moran,
early-topology, numerical, collection, plotting and test commands. Root and
branch sizes are passed separately; all positions use the physical
marker-to-position integral. No constant-frequency single-population formula
substitutes for the original multispecies Moran operator.

Next steps are to refine unresolved intervals, extend explicit history/origin
adapters and boundary specifications, and assess spatial joint dependence.
A position-varying SMC′ kernel is a separate research task.
