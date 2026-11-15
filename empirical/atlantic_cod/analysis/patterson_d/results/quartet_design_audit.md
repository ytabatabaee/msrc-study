# Frozen quartet and sample audit

The four primary designs in `config/frozen_quartets.tsv` were frozen before calculating nucleotide D. Their P1/P2 sister relationships were chosen from the independently defined outside-inversion population-history baseline; arrangement states were taken from the repository's existing population arrangement table. Exact sample IDs are in `config/frozen_sample_membership.tsv` and are independently cross-checked against the downloaded VCF headers in `results/vcf_sample_inventory.tsv`.

Each quartet uses 2 individuals per ingroup population and one diploid outgroup individual. The primary outgroup is Borsai (*Boreogadus saida*), because the upstream `cod_phylogenomics/src/run_dsuite.slurm` explicitly renames Borsai to `Outgroup`. Gadoga (*Gadus ogac*) is a predeclared alternative from the three external labels in the baseline tree. Per-site Borsai or Gadoga calls are used only when homozygous; missing and heterozygous calls do not polarize.

| LG | P1, P2 baseline sister | P3 | O | Arrangement contrast | Key confounder |
|---|---|---|---|---|---|
| LG01 | avc, lfo | avo | Borsai | P2/P3 both derived | More vs Lofoten geography and migratory ecotype covary |
| LG02 | avc, lfo | avo | Borsai | P2/P3 both ancestral | Geography/ecotype and arrangement are not crossed |
| LG07 | lfo, avc | avo | Borsai | P2/P3 both ancestral | P1/P3 share migratory ecotype; P2/P3 share Møre geography |
| LG12 | bor, kie | low | Borsai | P2/P3 both derived | Three distant stationary populations; Bornholm geography is distinctive |

Population labels use the established prefix mapping from the original pipeline, not local-tree tip assignment. The SNAPP VCFs have 29 sample IDs, the same gadMor2 reference and corresponding linkage-group positions; outside-supergene VCF calls are complete. Dsuite inversion VCFs contain a different Lofoten individual (`Gadmor_lfo3`) and missing calls, which is preserved as an explicit sensitivity roster. The local trees contain 12 cod-population tips and no Borsai/Gadoga outgroup, so an exact 4-tip `D_topology=q3-q2` is not calculable from them.

The biological basis is an allele-sharing contrast conditional on a rooted population-history quartet. No quartet was selected to maximize nucleotide D. A nonzero D is not uniquely diagnostic of introgression or MSRC.
