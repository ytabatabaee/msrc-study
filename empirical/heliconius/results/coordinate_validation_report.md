# Heliconius inversion-coordinate validation

Validated against the Edelman et al. (2019) local-tree publication/data record and Thawornwattana et al. (2023), Supplementary file 1, Table 2.

## 2b

Table 2 contains two entries labeled 2b under different datasets:

- `etales-8spp`: `Hmel202001o:2147755–2230606` and `Herato0206:880283–1059137`.
- `hmelv25-res`: `Hmel202001o:4398947–6357046` and `n/a` for H. erato.

These are not interchangeable coordinates. They are retained as distinct dataset-specific definitions. The first pair is the only cross-assembly pair in the table and is the one used for the available erato-clade Hmel-coordinate trees. The H. erato interval is converted to H. erato demophoon v1 chromosome 2 coordinates using `Herato.transitions.tsv`: scaffold 6 begins at chromosome coordinate 8,233,830, giving 9,114,112–9,292,966 (1-based inclusive).

The available `erato_50kb_hmel` local trees cover chromosome 2 through at least 9.025 Mb and therefore physically cover both Hmel2.5 intervals. Coverage alone does not make the trees an `hmelv25-res` dataset, so the 1.96-Mb interval is not used to classify them.

## 15b

Table 2 gives `Hmel215003o:1356511–1757545` only for `hmelv25-res` and reports `n/a` for H. erato in the relevant coordinate table. Edelman et al. describe the erato-clade chromosome-15 inversion using Hmel2.5-based contig alignments, but the available local-tree files do not provide an independently validated H. erato chromosome-15 projection. No H. erato 15b mask is applied, and the Hmel2.5 15b interval is not applied to the current erato-clade Hmel-coordinate trees because their dataset identity is not `hmelv25-res`.

## Biological interpretation

The ~179-kb H. erato 2b interval and ~83-kb Hmel2.5 `etales-8spp` interval are the paired cross-assembly definition in Table 2. The ~1.96-Mb Hmel2.5 interval is a separate `hmelv25-res` definition; this analysis does not claim that it is the same physical boundary set. The coordinate manifest records scaffold names, assemblies, source applicability, and unresolved mappings explicitly.

The Hera seven-taxon trees contain `Hdem`, `HeraRef`, `Hhim`, `Hhsa`, `HmelRef`, `Hsar`, and `Htel`. The Hmel-coordinate erato-clade trees contain `Etal`, `Hdem`, `HeraRef`, `Hhim`, `Hhsa`, `Hsar`, and `Htel`. Thus the erato-sara taxa discussed for 2b/15b are represented, but the coordinate table’s `hmelv25-res` dataset identity is not established for the downloaded Hmel-coordinate tree file. The structural arrangement history and the sampled local-tree taxon set should therefore not be treated as identical evidence.

Sources: [Edelman et al. 2019](https://doi.org/10.1126/science.aaw2090), [Dryad record](https://doi.org/10.5061/dryad.b7bj832), [Thawornwattana et al. 2023](https://doi.org/10.7554/eLife.90656).
