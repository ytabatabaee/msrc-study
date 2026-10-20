# Full-taxon data and estimator audit

The chromosome-4 local-tree universe contains 363 distinct observed taxon labels across 6140 unique loci. This is the empirical expected universe used here; no absent 363-species labels were invented.
The source file contains up to 363 observed taxa per local tree. The complete expected-universe count is 232 loci.
Stage-6 exact mapping has 47 exact entries; 47 occur in the local-tree universe and 0 are absent. Stage-6 complete-case regression used the observed fixed set.
The fixed broad set is defined before inspecting signal as taxa present in >=95% of chr4 local trees: 108 taxa.
B_observed is the mean patristic distance over all distinct observed local-tree taxa. B_fixed requires complete representation of the fixed broad set. B_residual is a pair-standardized residual centered by pair-specific strict-outside-region medians; it is not an absolute distance.
All distances are ML substitution units. They are not calibrated divergence times or coalescence times. Missing taxa are retained as missing observations and are never filled by interpolation or label guessing.
Stage-6 B regression maximum absolute difference over matched complete rows: 0.

## Estimator comparison

Across windows where both estimators are defined, B47 and B_observed have
Pearson correlation 0.9329 (2,707 windows); B_observed and B_fixed have
correlation 0.8491 (2,943 windows). Their means are 0.2098, 0.2537, and
0.1908 substitution units, respectively, over their available windows.
B_residual is a centered pair-standardized quantity with mean 0.0044 and is
not interchangeable with any absolute B track.

The canonical-region summaries are descriptive because windows are spatially
dependent and taxon composition varies. For example, B_observed minus its
outside-region mean is -0.0085 in CR_01, +0.0800 in CR_10, and +0.0167 in
CR_11. These values do not establish a structural causal effect and should
not be interpreted as independent biological replication.
