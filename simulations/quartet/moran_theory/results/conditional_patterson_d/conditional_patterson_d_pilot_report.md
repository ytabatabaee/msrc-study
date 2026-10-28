# Conditional MSRC Patterson's D pilot report

## Design

The pilot used the existing `msrc-sim` implementation, a pectinate tree
`(((P1,P2),P3),O)`, `M=20`, neutral Moran history generation, and effective
background-switching rate `rho=0.2`. The terminal arrangement pattern `0101`
was specified before simulation. Histories were generated and terminal
arrangements sampled until the target pattern occurred; the accepted history
and sampled arrangement were then fixed while 200 independent genealogies were
generated. This is rejection sampling from the valid history distribution
conditional on the prespecified arrangement pattern, rather than selection on
the observed genealogy or D value.

Three independent histories were retained. Their acceptance attempts were 27,
22, and 45. The standard MSC baseline used 600 independent genealogies on the
same rooted species tree.

The local run took 1.91 seconds with a recorded peak resident set size of
64.6 MiB. These measurements include history rejection, genealogy simulation,
bootstrap resampling, and mutation-count generation.

## Results

| Condition | q1 | q2 | q3 | branch D (ABBA-BABA) | topology q2-q3 | branch q2-q3 | sequence D |
|---|---:|---:|---:|---:|---:|---:|---:|
| MSRC history 1 | 0.235 | 0.500 | 0.265 | -0.461 | 0.307 | 0.461 | -0.420 |
| MSRC history 2 | 0.280 | 0.450 | 0.270 | -0.187 | 0.250 | 0.187 | -0.176 |
| MSRC history 3 | 0.280 | 0.415 | 0.305 | -0.159 | 0.153 | 0.159 | -0.239 |
| Standard MSC | 0.333 | 0.338 | 0.328 | 0.022 | 0.015 | -0.022 | 0.097 |

The branch-length values are conditional estimates from 200 genealogies per
MSRC history. Genealogy-bootstrap standard errors for branch D were 0.089,
0.114, and 0.115 for the three MSRC histories and 0.071 for the MSC baseline.
The site estimates used 100 mutation replicates per genealogy and were based
on 84--108 informative ABBA/BABA sites in the MSRC rows; they are feasibility
checks rather than independent-locus inference.

All three accepted MSRC histories showed the same negative ABBA-minus-BABA
sign, with the corresponding positive q2-minus-q3 topology diagnostic. The
MSC branch estimate was close to zero relative to its genealogy-bootstrap
uncertainty. The sequence estimates followed the branch-length predictions in
sign, as expected, but their small mutation counts are too noisy for a formal
comparison.

## Interpretation and limitations

This pilot demonstrates that a realized, biologically generated MSRC history
can produce a nonzero conditional ABBA/BABA branch-length asymmetry without an
introgression event in the simulator. It does not establish that the
expectation over all histories is nonzero, nor that the effect persists after
averaging over terminal arrangement configurations. Three histories are also
insufficient to estimate between-history variation reliably.

The statistic is a quartet branch-length asymmetry. It is not the same as the
topology diagnostic `q2-q3`, and neither statistic alone identifies
introgression. Outgroup rooting is explicit in the branch mapping: ABBA is
the `P2,P3` internal branch and BABA is the `P1,P3` internal branch, with `O`
as the fourth taxon. No mutations on terminal or non-two-two branches are
counted in the branch prediction.

The complete outputs are `pilot_results.tsv`, `pilot_manifest.json`, and the
PDF/PNG figure in this directory. A larger study should sample many
independent histories under a predeclared arrangement-pattern design and use
history blocks, rather than linked mutations, for uncertainty.
