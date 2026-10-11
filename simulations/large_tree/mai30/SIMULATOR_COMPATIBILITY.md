# `msrc-sim` large-tree compatibility audit

Audited checkout: `/Users/ytabatabaee/Desktop/msrc-sim` at commit
`f5157d2c64068eb32e0eabfe9bec629c28a404f9` (2026-09-20). The checkout is
dirty and contains local generated outputs. This validation adds a narrow
ultrametric-serialization tolerance and cached frequency-history lookups in
the working tree; record those changes in a clean simulator commit before
production.

## Capability matrix

| Requirement | Finding | Evidence / limitation |
|---|---|---|
| Arbitrary rooted species trees | Supported with serialization tolerance | `SpeciesTree` parses rooted Newick; the Mai trees are ultrametric within SimPhy's eight-decimal serialization precision. A narrow `ultrametric_tolerance` option was added and used at `1e-6` generations for Mai inputs. |
| Branch-specific `N_e` | Supported | `species_tree.branch_parameters.<branch>.effective_population_size`; both arbitrary MSRC and MSC paths read the branch value. |
| Forward WF histories across many species | Supported but generation-expensive | `simulate_population_history` traverses all species-tree branches and writes one record per generation. The measured scaled benchmark extrapolates to roughly 12.9 days for 500 trees at one selected original-scale Mai replicate. |
| Rearrangement on selected ancestral branch | Supported | `rearrangement.origin_branch` and `origin_time_from_branch_start`. |
| Transmission / segregation at speciation | Supported | Existing daughter initialization and optional `speciation_split` kernel; nonzero sorting requires binary splits. |
| Recombination suppression | Supported | `recombination.effective_cross_arrangement_fraction`; in structured coalescent the effective switch rate is the baseline rate times this fraction. |
| Full n-taxon genealogy | Supported | `simulate_genealogy` builds a complete `GenealogyResult` and validates all expected taxa once. |
| Full Newick export | Supported | `local_gene_trees.nwk`, one Newick record per simulated genealogy block. |
| Multiple loci sharing a rearrangement history | Partially supported | One forward history is sampled per arbitrary-tree run and reused across emitted blocks. |
| Exact requested gene count | Validated for direct MSC controls and the scaled pilot | The Mai MSC wrapper and scaled-WF pilot assert the requested count, exact taxa, bifurcation, finite branch lengths, and 29 coalescences per genealogy. The arbitrary-tree breakpoint builder now preserves forced rearrangement boundaries while subdividing spans to the requested block count. |
| Ordinary MSC control | Supported | `simulate_msc_genealogy` uses the same species tree and branch-specific population sizes, but does not use the MSRC history. It continues coalescence indefinitely in the root population; the configured root extension does not truncate the ordinary MSC. |
| Reproducibility | Supported | Seeded NumPy RNG and JSON config provenance; production runs must record a clean git commit and dependency lock. |

## Current decision

No separate coalescent simulator is needed for full Mai trees. The existing
arbitrary-tree genealogy engine is the correct implementation path. The Mai
MSC wrapper adds a hard output-contract check requiring exactly 500 full trees,
the exact Mai taxon set, one occurrence of every taxon per tree, finite
nonnegative branch lengths, and `n_taxa - 1` coalescences. The archive,
demographic mapping, root treatment, and ordinary MSC validation are complete;
the remaining blocker is the computational cost of original-scale forward
Wright–Fisher MSRC histories.

The Mai MSC controls generated 500 complete 30-taxon trees for four
observed-discordance representatives. Their mean normalized RF differences
from the original 500-tree collections were small relative to the observed
discordance gradients, but this is a distributional validation rather than a
claim of identical random realizations. Original SimPhy gene-tree branch
lengths are substitution lengths, so the current comparison is topological;
coalescent-time distributions require the SimPhy parameter database's
`Tree_h_cu` values or another explicitly recovered time field.

The Mai validation added a cached/bisected lookup for forward frequency-history
records in `wright_fisher.py`; this removes repeated sorting and linear scans
from the backward genealogy loop but does not change the biological process.

The scaled-WF pilot uses the same full-tree engine with `c=0.001`. Six primary
500-tree collections passed the complete output contract; six additional
conditions have validated 20-tree smoke collections pending cluster runtime.
