# SimPhy to `msrc-sim` population-size mapping

## Confirmed convention

The Mai D1 detailed SimPhy logs label `-sp` as **haploid effective population
size** and label species-tree times as **generations**. SimPhy's manual also
describes `-sp` as the tree-wide effective population size and its MSC model as
a Wright–Fisher coalescent. The D1 records use one sampled value per replicate;
the recovered files do not contain branch-specific population sizes.

The ordinary coalescent pair rate represented by the SimPhy haploid count
`N_H` is:

\[
  \lambda_{\mathrm{SimPhy}} = 1/N_H.
\]

`msrc-sim` uses `1/(2*N)` for ordinary MSC coalescence and uses `2*N` as the
Wright–Fisher chromosome count. Equating rates gives:

\[
  N_{\mathrm{msrc}} = N_H/2,
  \qquad
  t_{CU} = g/(2N_{\mathrm{msrc}}) = g/N_H.
\]

This is implemented by
`analysis/scripts/population_mapping.py`; the division is not duplicated in
the analysis scripts.

## Integer handling

The continuous MSC mapping is exact even when `N_H` is odd. The current
`msrc-sim` `SpeciesTree` and Wright–Fisher implementation require an integer
effective population size, so a finite-history run uses the nearest of
`floor(N_H/2)` and `ceil(N_H/2)`, selected by the smaller pair-rate error.
The mapping object records both the exact target and the integer value plus
the signed relative rate error. No odd value is silently truncated.

For ordinary MSC validation, the integer approximation is used only because of
the current simulator interface. The resulting rate error is at most about
`1/N_H` and is reported with every converted record. A simulator change that
accepts a floating-point coalescent `N_e` could use the exact value for MSC,
while finite Wright–Fisher copy counts would still require an integer.

## Branch and time treatment

The conversion is applied per branch when branch-specific sizes are supplied:

\[
  t_{CU,b}=g_b/(2N_{\mathrm{msrc},b}).
\]

The Mai D1 audit recovered a single global `N_H` per replicate, so all branches
of that replicate use the same mapping. Species-tree Newick lengths are in
generations according to the SimPhy logs; they are not rescaled by the
substitution rate. Generation time is fixed at 1.0 in every D1 command.

## Direct implementation checks

`analysis/tests/test_population_mapping.py` verifies:

- an even haploid size maps exactly and gives `g/N_H` coalescent units;
- an odd haploid size retains the exact half-integer target and reports the
  finite-history approximation error;
- conversion from coalescent units back to generations is consistent.

The simulator source checks used for this audit are
`msrcsim/structured_coalescent.py` (`1/(2*N)` MSC rate),
`msrcsim/wright_fisher.py` (`total=2*N`), and
`msrcsim/moran.py` (`total=2*N`). The local SimPhy binary's help output names
`-SP` as “tree-wide effective population size,” while the recovered D1 logs
provide the decisive haploid label.
