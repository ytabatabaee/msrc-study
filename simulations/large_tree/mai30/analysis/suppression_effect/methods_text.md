# Methods

We reanalyzed the completed four-replicate Mai30 scaled-Wright–Fisher pilot.
Each replicate has 500 gene trees for ordinary MSC, MSRC κ=1, κ=0.5, and
κ=0.1. The primary suppression reference is MSRC κ=1 within the same Mai
backbone and treatment design. Ordinary MSC is retained as an external
control.

Gene-tree summaries were read from the frozen pilot collections and the
existing recalculated tables. ASTRAL outputs were produced with
`-t 2 -u 1 --length CULength`. Normalized RF uses the fully bifurcating
30-taxon denominator `2(n-3)`. Branch-length summaries use only common true
and inferred internal splits; branch lengths on incorrect splits are excluded.
True lengths are read from the scaled pilot species tree and divided by
`2*N_scaled`; the raw SimPhy tree remains in generation units and is not used.
Per-split tables include signed, absolute, and relative errors.

For uncertainty, the analysis resamples complete 20-locus MSRC histories. Five
bootstrap collections were generated per replicate and treatment as an initial
runtime-limited analysis, and ASTRAL was rerun for every bootstrap collection.
The ordinary MSC controls were partitioned into sequential 20-tree
pseudo-blocks. These blocks provide a comparable gene-tree resampling check,
but they are not independent forward frequency histories and are labeled as
such in `astral_bootstrap.tsv`.

The split analysis matches clades by descendant taxon sets. It reports every
true split as recovered or missing and every inferred split absent from the
true species tree as an incorrect split. The implementation does not pair
MSRC and MSC histories or infer statistical significance from the four pilot
replicates.
