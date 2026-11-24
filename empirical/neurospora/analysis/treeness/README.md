# neurospora treeness analysis

## Method

Treeness is internal branch-length sum divided by total tree length, following the Systematic Biology CASTLES/CoalBL definition. For rooted trees, internal edges are non-root edges leading to nonterminal nodes; for unrooted trees, internal edges split at least two retained tips on each side. Terminal edges lead to one tip. Root stems are excluded; no missing edge length is set to zero. Zero-length edges remain zero. Missing/nonfinite/negative branches are counted and trees with any such branch have NA treeness; a zero-total tree also has NA treeness. Internal and terminal branch sums are retained separately.

## Inputs and provenance

Tree sets: L10_50SNP_local_BIONJ. Coordinate-indexed local trees processed: 45303; valid: 32472; invalid/zero-total: 12831. The per-window file records tree set, assembly chromosome/coordinates, branch-length units, tree method/type, rooting convention, taxon-set hash, branch-length counts, zero internal/terminal branches, and explicit outgroup sensitivity where taxon labels supported it. No species-wide reference tree is used. Posterior tree samples were not available as coordinate-matched samples in the included local-tree sources; these are point trees, not posterior summaries.

Suppressed-recombination masks used when validated and source-compatible: LGI_1p5_8Mb. Boundaries are not transferred across assemblies or source-specific coordinate axes. Region labels from unvalidated/approximate masks are not treated as validated comparisons.

## Comparisons

See `results/inside_outside_summary.tsv` and `results/branch_component_comparison.tsv`. Contrasts are descriptive median/mean differences. Mean-difference intervals use a 2-Mb moving-block bootstrap over genomic runs; windows are not treated as independent biological replicates. Inference is not reported where there are fewer than three observations in either class. Multiple chromosomes/source groups remain separate. Trees with nonmodal taxon sets are excluded from inside/outside estimates and counted in `n_excluded_nonmodal_taxa`; full tracks retain every valid tree.

Outgroup-pruned values are provided only where the outgroup was explicitly identifiable from source labels/metadata; otherwise no outgroup is guessed. Rooted/unrooted sensitivity can be assessed by the unrooted edge-split rule; root stem exclusion is invariant. Branch-length rescaling is not applied. See `tree_set_audit.tsv` for taxon consistency and outgroup influence.

## Reproduction

From repository root: `python3 empirical/neurospora/analysis/treeness/scripts/run_treeness.py`. Shared math is in `empirical/shared/treeness_core.py`; per-dataset results and figures remain here. Cross-dataset summary is `empirical/results/treeness_cross_dataset_summary.tsv`. Tests: `python3 -m unittest discover -s empirical/shared/tests`.
