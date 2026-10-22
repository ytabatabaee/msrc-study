# Optimized simulator and quartet-reversal report

## Completed

- Preserved the original backward implementation as `simulate_genealogy_reference`.
- Added an exact grouped-hazard implementation used by `simulate_genealogy`.
- Added regression tests and a benchmark script in `msrc-sim`.
- Measured approximately 2.15x speedup on one 30-taxon Mai pilot genealogy.
- Ran a 20-history, three-scenario exploratory quartet pilot for replicate 007.

## Not completed

The 1,000-history all-tree forward validation and an independent quartet
reversal validation were not completed locally. The older-root scenario and
long branches make forward history generation expensive; high-rho backward
runs add another event-count bottleneck. No full 30-taxon origin-placement
collection or ASTRAL analysis was generated.

## Scientific conclusion

The optimization makes the backward genealogy simulation feasible for small
and moderate pilots without changing the specified model. It does not yet
establish a dominant discordant quartet. Structural exposure and forward
persistence remain separate from genealogical distortion, and the current
results must not be used to claim quartet reversal or ASTRAL bias.
