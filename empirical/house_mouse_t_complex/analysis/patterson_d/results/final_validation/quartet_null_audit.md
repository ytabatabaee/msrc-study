# Frozen quartet null audit

## Decision

All three designs are correctly rooted *four-population contrasts*: the eight *Mus spretus* samples provide the external branch and P1/P2/P3 define ABBA/BABA orientation. At broad subspecies scale, each has P1 and P2 from the same subspecies and P3 from another, so the assumed `((P1,P2),P3),O` is a plausible species-history null for each. This does not establish a panmictic D=0 null: t-carrier selection, geography, within-subspecies structure, introgression, and ancestry can differentiate P1 from P2. All remain exploratory arrangement-associated population contrasts, and none is t-chromosome-specific.

## Sample/provenance audit

The frozen member table and source-backed mapping give D1 P1/P2 as 4/12 *M. m. domesticus* from France/Germany, P3 as 13 *M. m. musculus* from Afghanistan/Czechia/Kazakhstan, and O as eight *M. spretus* from Madrid. D2 is 8/13 *M. m. musculus* (t carriers sampled in Afghanistan, Czechia, Kazakhstan; standards from the same broad origins), 12 *M. m. domesticus* P3 from the France/Germany cohort, and the same eight outgroups. D3 is 3/7 *M. m. castaneus* from Himachal Pradesh, 13 *M. m. musculus* P3, and the same eight outgroups. The exact VCF/tree membership audit reports all 55 frozen IDs present in both sources; outgroup status is explicit in the frozen design. Verification used `config/frozen_quartets.tsv`, `config/frozen_quartet_members.tsv`, `exact_quartet_comparison/results/sample_mapping_audit.tsv`, and the indexed VCF header; arrangement labels were read from the frozen member manifest, not inferred from D. These are labels/provenance, not evidence that a population branching order is exchangeable.

## Biological interpretation of null

Population phylogenies support rapid divergence and gene flow among the three major house-mouse subspecies; several genome-wide analyses place *M. m. musculus* and *M. m. castaneus* closer to one another than either is to *M. m. domesticus*. That relationship does not contradict this quartet null: in every design P1 and P2 are two groups within one subspecies (DOM in D1, MUS in D2, CAS in D3), and P3 is another subspecies. The species tree alone therefore cannot reject P1/P2 as sisters. It also cannot validate exchangeability of selected t carriers and standards within those subspecies. D1 uses France/Germany domesticus carriers and standards and a geographically broad MUS P3; D2 uses AFG/CZE/KAZ musculus carrier and standard cohorts and France/Germany domesticus P3; D3 uses Himachal Pradesh castaneus carriers and standards and a geographically broad MUS P3. This gives some within-subspecies geographic matching but leaves population structure, ancestry, introgression, and selection-induced sampling as confounders. D3 is especially sensitive to its three carrier individuals, but its coarse species-history null is not invalidated by the CAS+MUS sister relationship. Recent Eurasian whole-genome work finds complex admixture, especially CAS–MUS, and topology sensitivity to hybrid-sample exclusion.

The samples are selected t carriers versus standard arrangement individuals within P1/P2. That is an arrangement-associated population contrast, but VCF genotypes are unphased diploid calls, so it does not isolate t-bearing chromosomes. Exact IDs do not make pseudo-t tree tips the same chromosome copies as the diploid VCF calls.

## Classification

| Design | Rooted quartet | Conventional species null | Structure/admixture concern | t-haplotype-specific? | Use |
|---|---|---|---|---|---|
| D1 Mmd t vs standard | Yes, with *M. spretus* outgroup | Plausible at subspecies level; not a validated panmictic null | Moderate: regional ancestry and cross-subspecies P3 | No | Exploratory arrangement-associated |
| D2 Mmm t vs standard | Yes, with *M. spretus* outgroup | Plausible at subspecies level; not a validated panmictic null | Moderate: geographic strata, ancestry, and cross-subspecies P3 | No | Exploratory arrangement-associated |
| D3 Mmc t vs standard | Yes, with *M. spretus* outgroup | Plausible at subspecies level; not a validated panmictic null | High: only 3 carriers; local ancestry and admixture remain possible | No | Exploratory arrangement-associated |

## Sources and scope

This classification is independent of chr17 D and local-tree results. It uses the source-backed sample IDs and geographic metadata plus published population history: [White et al. 2009, genome-wide house-mouse speciation history](https://pmc.ncbi.nlm.nih.gov/articles/PMC2915842/); [Phifer-Rixey et al. 2022, wild-mouse population structure](https://pmc.ncbi.nlm.nih.gov/articles/PMC9411160/); [Fujiwara et al. 2022, Eurasian whole-genome structure and admixture](https://academic.oup.com/gbe/article/14/5/evac068/6582302); and [Duvaux et al. 2011, isolation and gene flow in European house mice](https://pubmed.ncbi.nlm.nih.gov/22066696/). *M. spretus* is a separate species and established close outgroup for *M. musculus* in the source design; using it roots the quartet but does not establish the ingroup null. No quartets were changed post hoc.
