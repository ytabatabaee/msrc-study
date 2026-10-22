# Mai30 MSRC mechanism audit

The retained pilot uses the hierarchical model in `msrc-sim`: a forward
Wright--Fisher frequency history is generated once per history, terminal
arrangements are sampled from that history, and independent backward
genealogies are generated conditional on it.  The 20 loci sharing one history
are therefore clustered observations.

For all pilot collections the configured baseline intensity is

```text
rho0 = 0.02
rho  = rho0 * kappa = 0.02, 0.01, 0.002.
```

The `c=0.001` scaling is applied to the Mai population sizes and branch times;
the simulation operates on the scaled Wright--Fisher generation-equivalent
tree.  The backward switch rates are `rho*p` for 0→1 and `rho*(1-p)` for
1→0.  Same-background coalescence uses the finite-population rates in
`structured_coalescent.py`, with the ordinary panmictic fallback after the
rearrangement-origin boundary.

The pilot does not place the rearrangement on the root population.  The
origin branches are representative internal branches: N8 for replicate 007,
N23 for 042, N20 for 074, and N28 for 005.  The rearrangement therefore
affects only the descendant clade of that branch.  Lineages outside that clade
never experience the structured background process.  This is a major reason
that a 30-taxon ASTRAL topology can be insensitive even when affected-clade
genealogies show switching.

The root extension is a simulated ordinary ancestral extension after the
species-tree root; it is not evidence that the rearrangement is present across
the entire 30-taxon tree.  The biological structured interval is bounded by
the configured origin branch/time.  The current simulator then resolves older
unstructured ancestry under the ordinary MSC process.

Relative rates are not determined by rho alone.  At p≈0.5, a single-lineage
switch rate is about rho/2, while a same-background pair coalescence rate is
approximately `1/(2*Ne*p)`.  The integrated switch hazard depends on branch
duration and the realized frequency path.  The retained diagnostics show that
the three pilot rates do produce materially different switch counts, so the
weak ASTRAL contrast is not explained by identical switching histories.

The retained outputs include `arrangement_history.tsv`, terminal states, and
per-genealogy total switch counts.  They do not include the event time and
population branch for each backward switch or the ordering of each
coalescence relative to switches.  Those quantities are therefore marked
unavailable in the diagnostic tables rather than inferred from final trees.
