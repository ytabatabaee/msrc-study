# Large-tree runtime benchmark

## Benchmark setup

The benchmark uses Mai condition
`outgroup.0.species.1.5.genes.1.5`, replicate `042`, which has 30 taxa,
SimPhy `N_H = 68,955`, and a species-tree root at 16,022,324 generations.
The exact mapped ordinary-MSC size is 34,477.5; the finite Wright–Fisher
benchmark uses the documented integer value 34,478.

The existing Wright–Fisher implementation simulates one generation at a time.
To obtain measurable pilot points, both species-tree generations and the
population size were multiplied by the same factor `c`, preserving
`g/(2N_e)` in the continuous approximation. Each case generated one forward
history and 20 complete 30-taxon gene trees, with a 5-CU root-history
extension. The results are in
`analysis/results/mai30_runtime_benchmark.tsv`.

| c | scaled `N_e` | history seconds | seconds / gene tree | estimated seconds / 500 genes | peak RSS* |
|---:|---:|---:|---:|---:|---:|
| 0.0001 | 3 | 0.041 | 0.210 | 105 | 72 MB |
| 0.0005 | 17 | 0.200 | 1.081 | 540 | 94 MB |
| 0.0010 | 34 | 0.459 | 2.223 | 1,111 | 137 MB |
| 0.0020 | 69 | 0.917 | 4.503 | 2,252 | 209 MB |

\*RSS is the process high-water mark and is cumulative across benchmark cases;
it should be treated as an upper bound for this sequential benchmark, not as an
isolated per-case allocation.

The measured cost is approximately linear in `c` over these points. A simple
linear extrapolation from `c=0.001` gives about 1.11 million seconds, or 12.9
days, for 500 gene trees from one full Mai history at the original scale. The
estimate is uncertain because the scaled cases have very small populations
(`N_e=3–69`) and therefore do not reproduce all finite-population MSRC
dynamics. It is nevertheless sufficient to reject a naive full-scale,
generation-by-generation 500-gene pilot on a single ordinary workstation.

## What is and is not preserved

The transformation

\[
  N_e' = cN_e, \qquad g' = cg
\]

preserves branch lengths in coalescent units. It does not automatically
preserve the complete MSRC model. A faithful downscaled run would also need
to consider recombination rates, rearrangement-origin timing and copy counts,
selection coefficients, and mutation rates. In particular, the finite
Wright–Fisher chromosome count changes from `2N_e` to `2cN_e`, so very small
benchmark populations have stronger discrete drift and are computational
probes rather than biological substitutes.

`msrc-sim` already supports a continuous-time Moran process, but that is a
different evolutionary process and is not used as an MSC-equivalent
replacement in this validation. No new approximate simulator was introduced.

## Feasibility decision

The ordinary MSC control is computationally feasible: the direct coalescent
implementation generated 500 trees in roughly 1–2 seconds per selected
replicate. The full forward Wright–Fisher MSRC history at the original Mai
population sizes is currently computationally prohibitive for 500 genes per
history. Before an MSRC pilot, the study should use a small number of
coalescent-timescale-preserving scaled histories with explicit labels, or the
authoritative simulator should gain a validated acceleration for neutral
finite-population histories. Neither option should be presented as identical
to the original-scale MSRC dynamics.
