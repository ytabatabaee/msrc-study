# Quartet and sample-provenance audit

## Frozen quartet

The nucleotide estimator uses `((L10 mat-a, L10 mat-A), N. hispaniola), N. discreta)` with 10 verified mat-a and 8 verified mat-A isolates, 8817@999 as P3, and Ndis@999 as O. The D calculation is valid as an explicitly defined four-group site-pattern contrast because Ndis is present and coordinate-validated. It does not by itself establish that the L10 mating-type groups are distinct population branches across the whole genome.

Martin & Van Belleghem (2017) explicitly note that mat-a and mat-A L10 samples are shown as separate branches for topology weighting but, except within the nonrecombining LGI region, are representatives of one recombining population. This makes the classical no-gene-flow species-tree null and its expectation D=0 conditional on a stronger demographic/exchangeability model than the rooted quartet alone. A correctly polarized quartet is not proof that the null is biologically exact.

*N. hispaniola* is a defensible P3 for the published L10 hypothesis because Corcoran et al. report the closest similarity/introgression signal in L10 mat-A and report positive autosomal D involving L10 and *N. hispaniola*. *N. discreta* is an external species used for polarization in the prior study and is present in the separately supplied alignment. These choices are biologically motivated independently of the new D estimates.

## Known competing history

Corcoran et al. (2016) report evidence for introgression from *N. hispaniola* into L10 mat-A and signatures of autosomal introgression involving all L10 strains. They also discuss admixture involving L10, L4, and *N. hispaniola*. Martin & Van Belleghem selected L10 as an example with a known strong introgression-associated topology shift. This prior work makes introgression a plausible contributor to the present asymmetry; it does not establish that every high-D flank/outside segment has the same origin or that the current ascertained-panel D uniquely measures introgression.

Geographic and population history remain relevant. Corcoran et al. identify the sampled L10 lineage as UK. The available data do not establish that the 10/8 verified individuals are exchangeable at every genomic segment, nor do they remove within-lineage structure, lineage-specific history, or sample ascertainment as alternatives.

## Assignments

`sample_provenance_audit.tsv` reports all 22 CJ isolates. Eighteen assignments are supported by Corcoran et al. Figure 1 or a directly documented opposite-mating-type heterokaryon relationship (CJ02 derived as the alternate homokaryon from mat-A CJ01). The four unresolved CJ labels remain excluded. No assignment is inferred from numbering, allele clustering, or local trees.

## References

- Corcoran et al. 2016, [Introgression maintains the genetic integrity of the mating-type determining chromosome of *Neurospora tetrasperma*](https://doi.org/10.1101/gr.197244.115).
- Martin & Van Belleghem 2017, [Exploring evolutionary relationships across the genome using topology weighting](https://doi.org/10.1534/genetics.116.194720).
- Deposited analysis data: [Dryad record 10.5061/dryad.4jm83](https://doi.org/10.5061/dryad.4jm83); source population-genomic study: [Corcoran et al. Dryad record 10.5061/dryad.162mh](https://doi.org/10.5061/dryad.162mh).
