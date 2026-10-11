# Mai30 validation results

Status: **BLOCKED** for the MSRC pilot. The ordinary MSC and output-contract
checks pass, but original-scale forward Wright–Fisher histories are currently
too expensive for the planned 500-gene MSRC datasets.

## Dataset selection

The authenticated local `D1-true.tar.gz` archive was used; it was not
redownloaded. The archive checksum matches Dryad DOI `10.6076/D1RW2G`. The
initial validation used all 300 no-outgroup collections (three conditions ×
100 replicates), each containing one 30-leaf species tree and 500 original
gene trees.

Representative replicates were selected by measured mean unrooted normalized
RF distance to their own species tree, with the normalization

\[
  RF_{norm}=RF/[2(n-3)]
\]

for fully bifurcating unrooted n-taxon trees. The selected collections span
the observed discordance range:

| Condition | Replicate | Mean normalized RF | Exact species-tree topology | Mean sampled quartet concordance |
|---|---:|---:|---:|---:|
| `outgroup.0.species.1.5.genes.1.5` | 042 | 0.00370 | 0.902 | 0.99620 |
| `outgroup.0.species.5.genes.5` | 005 | 0.30289 | 0.000 | 0.82750 |
| `outgroup.0.species.1.5.genes.1.5` | 074 | 0.57889 | 0.000 | 0.64769 |
| `outgroup.0.species.1.5.genes.1.5` | 007 | 0.97896 | 0.000 | 0.34672 |

The selection and original baseline are recorded in
`datasets/processed/mai30_original_gene_tree_selection.tsv` and
`datasets/processed/mai30_original_gene_tree_baseline.tsv`. The original
results are treated only as SimPhy reference data, never as MSRC output.

SimPhy's documentation states that species and locus branches are measured in
generations, while gene-tree branches are expected substitutions per site.
The original gene-tree files also contain branch-specific substitution-rate
heterogeneity. Consequently, this validation compares topology and
discordance distributions; it does not compare original gene-tree branch
lengths as coalescent times.

## Population-size and demographic conversion

The D1 logs label `-sp` as haploid effective population size and species-tree
times as generations. `msrc-sim` uses pairwise MSC rate `1/(2N)` and
`2N` chromosome copies in Wright–Fisher histories. The implemented mapping is:

\[
  N_{msrc}=N_H/2,\qquad t_{CU}=g/N_H.
\]

The mapping is implemented in
`analysis/scripts/population_mapping.py` and tested in
`analysis/tests/test_population_mapping.py`. Odd `N_H` values retain an exact
half-integer target for MSC calculations and separately report the nearest
integer finite-history approximation and its rate error. The recovered Mai
files use one global population size per replicate; no branch-specific size
table was found.

## Root treatment

The ordinary MSC implementation continues coalescence in the root population
without a finite endpoint. The configured root extension does not truncate or
force coalescence. A same-seed test at 5, 10, 20, and 40 coalescent units
produced identical Newick hashes for 100 trees. This behavior is documented in
`ROOT_EXTENSION_VALIDATION.md` and measured in
`analysis/results/root_extension_same_seed.tsv`.

The forward rearrangement history is separate: its Wright–Fisher records stop
at the configured root extension and the oldest frequency is held afterward.
That finite history should be treated as a modeling choice in any later MSRC
pilot.

## Ordinary MSC validation

For each selected replicate, four root settings (`R = 5, 10, 20, 40` CU) were
used, with deterministic seeds and 500 newly simulated ordinary MSC gene
trees. Every set passed the output contract. Across root settings, the
generated-minus-original mean normalized RF differences were:

| Replicate | Original mean RF | Generated mean RF range | Difference range | Maximum sorted-RF ECDF distance |
|---:|---:|---:|---:|---:|
| 042 | 0.00370 | 0.00467–0.00593 | 0.00096–0.00222 | 0.00237 |
| 005 | 0.30289 | 0.30074–0.30778 | −0.00215–0.00489 | 0.00882 |
| 074 | 0.57889 | 0.57948–0.58630 | 0.00059–0.00741 | 0.00800 |
| 007 | 0.97896 | 0.97837–0.97919 | −0.00059–0.00022 | 0.00141 |

Sampled quartet and species-split concordance showed the same broad agreement.
The new trees are independent random realizations, so this is a distributional
validation rather than a requirement for identical trees. The complete result
table is `analysis/results/mai30_msc_validation.tsv`.

## Full-tree output contract

The MSC wrapper asserts, for every generated set:

- exactly 500 Newick records;
- the exact 30-taxon set, with one occurrence per taxon;
- finite, nonnegative branch lengths;
- fully bifurcating trees;
- exactly 29 coalescence events per genealogy;
- no malformed or silently truncated tree.

All 16 selected-replicate/root-setting datasets passed. The same validation
must remain in the future MSRC linked-block wrapper, where breakpoint
partitioning can otherwise change the number of emitted blocks.

## Runtime and feasibility

The ordinary MSC implementation generated 500 trees per selected replicate in
roughly 1–2 seconds. The forward Wright–Fisher benchmark used the observed
high-discordance replicate 042 and scaled both generations and `N` by the same
factor `c` to preserve coalescent-unit branch lengths. At `c=0.001`, 20 full
gene trees took 44.5 seconds, implying about 1,111 seconds for 500 trees. A
linear extrapolation to the original scale is approximately 12.9 days for one
500-gene history, before independent-history replication. The measured points,
memory, and assumptions are in `RUNTIME_BENCHMARK.md`.

The scaling transformation is a computational approximation: it preserves
coalescent branch lengths but changes finite chromosome counts and therefore
does not automatically preserve rearrangement drift, recombination, origin
frequency, selection, or mutation dynamics. The existing Moran implementation
is also a different process and was not substituted for Wright–Fisher.

## Remaining blockers and recommended next step

The demographic mapping, representative original data, root behavior, MSC
validation, and output contract are complete. The remaining blocker is
computational: the current generation-by-generation Wright–Fisher MSRC path is
not feasible at the original Mai population sizes for 500 loci per history.
The next execution should be a small explicitly labeled scaled MSRC pilot,
only after deciding which quantities are held fixed under scaling and adding a
runtime-aware linked-block output contract. It should not be presented as an
exact original-scale MSRC simulation.
