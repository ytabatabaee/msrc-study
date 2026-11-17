# Independent species-history null audit

## Scope and sources

The frozen nucleotide designs and orderings are retained exactly. The conventional D=0 interpretation requires that P1 and P2 are sisters in the rooted species quartet, with O outside the three ingroups, plus the usual MSC and polarization assumptions. Correct site counting does not by itself validate this null.

The locally available Stiller reference is `empirical/neoaves_chr4/external/stiller2024/63K.tre`, a published 363-taxon genome-wide reference tree. Its induced quartets were calculated from the reference Newick using exact taxon labels and unrooted pairwise path sums, then rooted with *Gallus gallus*. This reference uses genome-wide Stiller data, rather than only chr4, but it is not an independent dataset from the nucleotide alignments. The separate Jarvis et al. 2014 48-taxon TENT provides independent historical species-tree evidence ([Jarvis et al. 2014](https://doi.org/10.1126/science.1253451)); the repository contains its published supplement but not a machine-readable Newick, and *Ciconia maguari* is not in that sampling. Consequently Jarvis evidence is recorded as broad clade context, not as an exact induced quartet.

The Stiller paper reports a rapid early Neoaves radiation and substantial topology heterogeneity across genomic partitions ([Stiller et al. 2024](https://doi.org/10.1038/s41586-024-07323-1)). Alternative early Neoaves backbones are also represented in the literature, so one resolved tree is not certainty about the species history.

## Design-by-design result

| Frozen design | Stiller 363-taxon induced topology | Jarvis 48-taxon evidence | D=0 null assessment |
|---|---|---|---|
| Columbea: *Columba*, *Phoenicopterus*, *Ciconia*, *Gallus* | q2, pairing *Columba* with *Ciconia*, not frozen P1/P2 | Supports broad Columbea and pigeon/flamingo clade context, but *Ciconia maguari* is absent; exact quartet unavailable | **Disputed / not established.** The historical Columbea hypothesis is a reasonable exploratory question, but the exact sampled quartet is not stable across references. A nonzero D is not diagnostic of introgression under this design. |
| N61: *Tauraco*, *Columba*, *Ciconia*, *Gallus* | q1, the frozen P1/P2 pairing | Broad placements are available, but *Ciconia maguari* is absent; exact quartet unavailable | **Conditionally supported by Stiller, not independently confirmed by exact 48-tip sampling.** The observed baseline D is not zero, so a simple symmetric-ILS null is not a good description even if q1 is accepted. |
| N62: *Columba*, *Pterocles*, *Tauraco*, *Gallus* | q1, the frozen P1/P2 pairing | Broad clade membership is consistent with Columbimorphae within Columbea and Tauraco in Passerea; exact local Newick unavailable | **Best-supported of the three, but still conditional.** Galliform rooting and the broad sister-group assignment are reasonable; the empirical genome-wide D baseline is nevertheless nonzero. |

The machine-readable taxon identities, source, induced topology, broad independent evidence, uncertainty, and interpretability are in `quartet_topology_validation.tsv`.

## Outgroup

*Gallus gallus* is outside Neoaves in both the 363-taxon reference and the Jarvis 48-taxon design and is a defensible external outgroup for these sampled neoavian ingroups. *Meleagris gallopavo* is also Galliformes and was frozen as an alternative external galliform. Substitution of *Meleagris* changes the numerical estimates modestly but does not remove the N61 regional shift or the Columbea/N62 near-null regional contrasts. Two Galliform outgroups do not independently test long-branch or ancestral-state bias shared by Galliformes.

## Consequences for interpretation

- ABBA/BABA counts can be computationally correct while the rooted quartet null is wrong or uncertain.
- A significant D outside the PNAS intervals tests `D_outside = 0`; it does not test whether the focal interval changed D.
- `D_inside - D_outside` tests a different question. It is the relevant regional contrast, but it does not establish a causal recombination-suppression mechanism.
- Rapid radiation, ILS, ancestral structure, ancient introgression, topology uncertainty, and finite-sites homoplasy are plausible model complications. The current results do not identify which caused a nonzero background D.
- No reoriented quartet estimates were generated: the frozen three designs remain unchanged, and the D=0 interpretation is qualified rather than silently optimized.
