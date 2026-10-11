# `msrc-sim` large-tree compatibility audit

Audited checkout: `/Users/ytabatabaee/Desktop/msrc-sim` at commit
`f5157d2c64068eb32e0eabfe9bec629c28a404f9` (2026-09-20). The checkout is
dirty and contains local generated outputs; the commit and status are recorded
here so a clean simulator checkout can be used for production.

## Capability matrix

| Requirement | Finding | Evidence / limitation |
|---|---|---|
| Arbitrary rooted species trees | Supported | `SpeciesTree` parses rooted ultrametric Newick; examples include 10 and 50 taxa. The Mai audit contains 30-leaf trees without an outgroup and 31-leaf trees with the SimPhy outgroup; a direct smoke run against a Mai backbone remains to be recorded. |
| Branch-specific `N_e` | Supported | `species_tree.branch_parameters.<branch>.effective_population_size`; both arbitrary MSRC and MSC paths read the branch value. |
| Forward WF histories across many species | Supported | `simulate_population_history` traverses all species-tree branches and writes split events. Runtime and memory for 30 taxa still require measurement. |
| Rearrangement on selected ancestral branch | Supported | `rearrangement.origin_branch` and `origin_time_from_branch_start`. |
| Transmission / segregation at speciation | Supported | Existing daughter initialization and optional `speciation_split` kernel; nonzero sorting requires binary splits. |
| Recombination suppression | Supported | `recombination.effective_cross_arrangement_fraction`; in structured coalescent the effective switch rate is the baseline rate times this fraction. |
| Full n-taxon genealogy | Supported | `simulate_genealogy` builds a complete `GenealogyResult` and validates all expected taxa once. |
| Full Newick export | Supported | `local_gene_trees.nwk`, one Newick record per simulated genealogy block. |
| Multiple loci sharing a rearrangement history | Partially supported | One forward history is sampled per arbitrary-tree run and reused across emitted blocks. |
| Exact requested gene count | Not yet guaranteed for linked blocks | `num_loci` controls the requested block count only when no extra breakpoint process changes the block list; the current runner emits one tree per block and truncates/returns the block list. A production wrapper must assert exactly 500 trees. |
| Ordinary MSC control | Supported | `simulate_msc_genealogy` uses the same species tree and branch-specific population sizes, but does not use the MSRC history. The Mai source's haploid `-sp` convention must be mapped to the simulator's convention before validation. |
| Reproducibility | Supported | Seeded NumPy RNG and JSON config provenance; production runs must record a clean git commit and dependency lock. |

## Current decision

No separate coalescent simulator is needed for full Mai trees. The existing
arbitrary-tree genealogy engine is the correct implementation path. Before a
biological pilot, add a hard output-contract check requiring exactly 500 full
trees, exactly the Mai taxon set, one occurrence of every taxon per tree,
finite nonnegative branch lengths, and `n_taxa - 1` coalescences. The Mai
archive and demographic records are now recovered; the remaining scientific
gate is explicit haploid/diploid mapping, root-extension treatment, and
ordinary MSC validation.
