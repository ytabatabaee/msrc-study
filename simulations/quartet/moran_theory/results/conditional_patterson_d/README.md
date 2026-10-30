# Conditional MSRC Patterson's D pilot

This directory contains a small feasibility pilot and isolated convergence
diagnostic. Three independent
Moran frequency histories were accepted by a predeclared terminal pattern
`0101`; each accepted history was then held fixed while 200 independent
genealogies were simulated. No history was selected using D or topology.

The branch-length statistic uses the explicit internal edge between the first
and second coalescences. ABBA is the `P2,P3` branch and BABA is the `P1,P3`
branch. The standard MSC row is an independent no-arrangement baseline.

The pilot is not a production estimate across the distribution of histories.
The `convergence_10x500/` result contains ten additional accepted histories
with 500 genealogies each; it is a local precision check, not the planned
production analysis. The planned 100-history/5,000-genealogy run is prepared
as a Bridges-2 array in
`../../cluster/conditional_patterson_d_expanded.sbatch` and has not been
submitted. Before submission, create the log directory because Slurm opens
log files before the batch shell starts:

```bash
cd simulations/quartet/moran_theory/cluster
mkdir -p logs
MSRC_REPO=/ocean/projects/cis260126p/shared/neoaves_stage4e/repo \
MSRC_SIM=/path/to/staged/msrc-sim \
PYTHON_BIN=/path/to/linux/python \
sbatch --array=0-9%2 conditional_patterson_d_expanded.sbatch
```

The parameter manifest lists four predeclared switching-rate batches. The
array script itself defaults to `rho=0.2`; use separate reviewed submissions
for the other rows rather than changing parameters after inspecting D.
Linked-site uncertainty and between-history uncertainty require the complete
history-block analysis.
