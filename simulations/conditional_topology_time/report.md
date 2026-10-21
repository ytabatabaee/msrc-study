# Conditional topology and coalescence-time signals under MSRC

## Scope

This study tests whether weak marginal MSRC effects can conceal stronger
effects in terminal arrangement strata. It reuses the authoritative
frequency-history and genealogy simulators and does not introduce a detector
or a second stochastic process.

## Design

The discovery stage used 12 prespecified cells and 500 independent forward
Moran histories per cell. The validation stage used fresh seeds and 5,000
histories for each balanced-frequency cell with M in {20,40} and rho in
{0,0.5}. One genealogy was sampled per history, so uncertainty is not
artificially reduced by treating conditional genealogies from a shared
history as independent. The complete grid and seeds are in each stage's
`manifest.json`.

## Theory check

The finite-M full-species-tree operator was evaluated for `P(Q,Z)`, `P(Z)`,
and `P(Q|Z)`. Joint-law agreement was strong: validation RMSE 0.00205 and
maximum absolute error 0.00572. Conditional errors were larger (RMSE 0.02865,
maximum 0.08266) because the pattern strata contain only tens to hundreds of
histories. This is consistent with Monte Carlo uncertainty and is why the
tables retain pattern counts and conditional MC standard errors.

## Findings

Some arrangement strata have visibly different mean pairwise coalescence
times. For example, 1010 versus 0011 differed by 1.90 time units at M=20,
rho=0 and by 2.66 at M=40, rho=0.5. The contrast changed sign or became
small in other validation cells. The marginal topology and time effects
remained weak: marginal total-variation topology differences ranged from
0.0036 to 0.0220, and standardized mean-time differences from 0.0048 to
0.0403.

The pattern-weighted decomposition exactly reconstructs the marginal
statistics. Thus dilution by low-probability patterns is a real mechanism in
these runs, but the data do not establish that dilution is the sole cause;
opposing pattern contributions also occur and must be inspected in the
decomposition table. Topology-conditioned time rows are exploratory because
their effective counts are smaller.

## Interpretation

The experiment supports a conditional, latent-stratum signal in some finite-M
settings. It does not establish a general MSRC theorem for genealogical-time
distributions, nor does it show that a detector without arrangement labels can
recover the signal. The current evidence justifies retaining true genealogy
time summaries as a simulation feature for future work, while requiring an
observation model and additional independent validation before presenting
branch-length-aware detection as a consequence of MSRC theory.

## Reproduction

Use the configured analysis interpreter and run:

```bash
python simulations/conditional_topology_time/scripts/conditional_experiment.py \
  --output simulations/conditional_topology_time/results/discovery \
  --histories 500 --seed 20261121 --stage discovery \
  --M 20 40 --initial-fraction 0.25 0.5 0.75 --rho 0 0.5 \
  --sister-length 0.25 --root-length 0.5

python simulations/conditional_topology_time/scripts/conditional_experiment.py \
  --output simulations/conditional_topology_time/results/validation \
  --histories 5000 --seed 20261131 --stage validation \
  --M 20 40 --initial-fraction 0.5 --rho 0 0.5 \
  --sister-length 0.25 --root-length 0.5

python simulations/conditional_topology_time/scripts/make_figures.py
```

The authoritative simulator is expected at `../msrc-sim` relative to this
repository, or can be passed with `--msrc-sim`.
