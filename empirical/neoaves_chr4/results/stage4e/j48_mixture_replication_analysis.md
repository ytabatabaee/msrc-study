# Neoaves J48 mixture replication analysis

## Design and validation

This report postprocesses 25 completed ASTRAL-IV runs: five independent with-replacement mixtures at each epsilon in 0.02, 0.05, 0.075, 0.10, and 0.15. Every mixture contains 62,945 gene-tree records. All 25 runs passed validation with 48 taxa, identical CU/SU unrooted topologies, valid output checksums, zero exit codes, and topology_su execution metadata.

## Topology recovery

J48_T0 first appears at epsilon = 0.075, becomes the majority topology at epsilon = 0.075, and all five replicates recover it at epsilon = 0.075. These are observed grid summaries, not a fitted threshold.

| epsilon | T0 | NONCHR4 | intermediate/distinct | novel-split replicates |
|---:|---:|---:|---:|---:|
| 0.020 | 0/5 | 0/5 | 5/5 | 1/5 |
| 0.050 | 0/5 | 0/5 | 5/5 | 2/5 |
| 0.075 | 5/5 | 0/5 | 0/5 | 0/5 |
| 0.100 | 4/5 | 0/5 | 1/5 | 1/5 |
| 0.150 | 4/5 | 0/5 | 1/5 | 1/5 |

## RF summary

| epsilon | T0 RF mean/median/range/SD | non-chr4 RF mean/median/range/SD |
|---:|---|---|
| 0.020 | 4.4/4/4–6/0.894 | 2.4/2/2–4/0.894 |
| 0.050 | 4.8/4/4–6/1.1 | 2.8/2/2–4/1.1 |
| 0.075 | 0/0/0–0/0 | 6/6/6–6/0 |
| 0.100 | 0.8/0/0–4/1.79 | 6.8/6/6–10/1.79 |
| 0.150 | 0.4/0/0–2/0.894 | 6.4/6/6–8/0.894 |

## Focal relationships

N61 is q1 in all five replicates at 0.02 and 0.05, then q2 in all replicates from 0.075 onward. Columbea changes from q2 to q1 at 0.075. N62 remains q1 across the sampled replicate grid. These are frozen biological states; they are not branch-local q1/q2/q3 labels.

Three exact novel split identities occur across the 25 trees. One recurs in four replicates and the other two occur once each; no novel split is present at 0.075, where all five trees equal J48_T0. Exact split strings and frequencies are in `j48_mixture_replication_novel_splits.tsv`.

## Interpretation and limitations

The experiment establishes sensitivity of a fixed 48-taxon species-tree inference to reweighting empirical genealogy pools. It does not establish recombination suppression as the cause, equate epsilon with a physical genomic fraction, estimate an MSRC parameter, or make Patterson-style or QQS crossover claims. The QQS values N61 approximately 0.1540 and N62 approximately 0.01354 are different estimands.

Five replicates per weight provide limited binomial precision. Adjacent weights are independent resamples rather than nested additions of loci. The smallest useful follow-up is additional seeds at weights with mixed recovery, especially 0.05, 0.10, and 0.15; fixed-candidate quartet scoring on the predeclared disagreement quartets would be a complementary analysis if independently validated.
