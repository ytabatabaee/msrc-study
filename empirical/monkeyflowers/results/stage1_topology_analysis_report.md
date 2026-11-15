# Stage 1 topology-only analysis report

## Recommendation

**CONDITIONAL GO for a topology-only empirical example.** The 100-kb TWISST weights are correctly mapped to 1,948 physical windows and show both a genome-wide mosaic of alternative genealogies and a strong, spatially sustained recent-introgression-associated signal at the independently reported chromosome 5 candidate interval. The pattern is shorter and less spatially coherent than the Neurospora example, with several monkeyflower runs reaching 0.6–1.4 Mb but none approaching Neurospora's longest 4.5-Mb run. This is an exploratory comparison, not a causal test of introgression versus recombination suppression.

This recommendation is limited to topology weights. The original local Newick trees are unavailable, so A(w), B(w), internal SU branches, branch-length comparisons, CU estimates and full gene-tree-based species-tree inference remain unsupported. Include the system only if the paper can clearly label these constraints and treat topology weights as TWISST support across windows, not independent one-quartet-per-window observations.

## Input and topology validation

There are 1,948 windows across LG1–LG10, each a 100-kb nonoverlapping interval on the *M. aurantiacus* reference assembly used by the study (Stankowski et al. 2019). All windows have finite topology weights and positive SNP/site counts. Coordinate keys match between the metadata TSV, the recombination-coordinate grid and the combined weight export. The separately deposited topology-count file lacks coordinates; its normalized 15-class rows match the combined export row-by-row, and the combined export then matches the metadata by genomic coordinate. This count-file alignment cannot be keyed independently because the file has no window IDs.

The paper's 15 topology labels were mapped by their Newick definitions, not column order guesses. `q1_species` is topo3. `q2_recent` is the combined weight of topo9, topo12 and topo13 (all have Island *longiflorus* and *parviflorus* as sisters). `q3_ancient` is topo10. The other 11 classes are preserved as `other_weight`. All 15 classes sum to 1 in each window. Mean normalized 15-class entropy is 0.081.

The file has 553 windows with complete species-history topology weight, 61 windows with a single recent-introgression topology at weight 1 (11 topo9, 38 topo12, 12 topo13), and 419 windows with topo10 weight 1. A further 915 windows have no one of these individual topologies at weight 1, matching the study's unresolved category. The set-level recent weight is exactly 1 in 68 windows: seven of those have their weight split among the three recent topologies rather than one individual topology at 1. This distinction reconciles the 61 exact-tree count with the 68 windows fully contained in the recent sister-pair topology set.

## Genome-wide support and sensitivity

Across all windows the unweighted mean weights are: species-history q1 = 0.3712, recent set q2 = 0.1364, ancient topology q3 = 0.3860, and other 11 topologies = 0.1064. `D(w) = q2_recent − q1_species` averages −0.2348. The ancient-introgression topology has the highest average individual topology weight, narrowly above the species-history topology. By per-window dominant class, 800 windows are ancient, 743 species-history, 217 recent-set member, and 188 other.

Excluding the two reported 100-kb chromosome 5 windows changes average weights by only +0.00038 for q1, −0.00085 for q2 and +0.00036 for q3; the dominant genome-wide class remains ancient. Excluding the approximate 300-kb sensitivity interval (three windows) changes q1 by +0.00057, q2 by −0.00130 and q3 by +0.00056. These are small genome-wide shifts. Topology-defined subsets are descriptive decompositions and do not test species-tree bias.

## Spatial runs and autocorrelation

A run requires consecutive observed windows with exact physical adjacency; runs are split at chromosome boundaries, coordinate gaps and missing windows. A run records topology-weight support only; its span is not an inferred ancestry-tract length. Threshold comparisons use strict `weight > threshold`.

| Topology-weight track | Threshold | Runs | Median run | Longest run | Transitions/Mb |
|---|---:|---:|---:|---:|---:|
| Recent set q2 | >0.5 | 131 | 0.1 Mb | 0.6 Mb | 1.34 |
| Recent set q2 | >0.7 | 106 | 0.1 Mb | 0.5 Mb | 1.08 |
| Recent set q2 | >0.9 | 73 | 0.1 Mb | 0.5 Mb | 0.74 |
| Ancient topology q3 | >0.5 | 401 | 0.1 Mb | 1.4 Mb | 4.11 |
| Ancient topology q3 | >0.7 | 394 | 0.1 Mb | 1.0 Mb | 4.04 |
| Ancient topology q3 | >0.9 | 351 | 0.1 Mb | 1.0 Mb | 3.60 |

The longest q2 runs are 0.6 Mb on LG3 (17.6–18.2 Mb) and LG10 (8.9–9.5 Mb). The longest q3 run is 1.4 Mb on LG5 (21.5–22.9 Mb). Thus most runs are a single 100-kb window, but a minority form longer coherent clusters, especially for q3. Threshold choice changes run counts substantially while leaving medians at one window and the ranking of longest runs unchanged.

At 100-kb midpoint separation, pooled Pearson autocorrelation is 0.412 for q2, 0.139 for q3 and 0.444 for D. The first relative 1/e crossing from the 100-kb value occurs at about 1.3 Mb for q2, 0.7 Mb for q3 and 1.2 Mb for D. These are descriptive crossing distances, not fitted decay parameters. Distance-bin pairs overlap heavily and are not independent; no confidence intervals or significance tests are assigned. The reported n-pair count is retained in the autocorrelation table.

## Chromosome 5 candidate region

The article identifies two adjacent 100-kb windows at 24.7–24.9 Mb with elevated Q95 ratios and f_d of 0.73/0.74, and describes the broader candidate as approximately 300 kb. Coordinates use the study's *M. aurantiacus* reference assembly. As the textual report does not define nucleotide-precise broader boundaries, the comparison reports the two published windows as primary and 24.7–25.0 Mb as an approximate 300-kb sensitivity interval.

Across the primary two windows, mean recent-set weight is 0.964, ancient weight 0.036, and species-history weight 0; both windows are recent-set dominant. Across the approximate 300-kb interval, recent-set weight averages 0.976. The recent-set >0.5 run continues through 25.1 Mb, spanning four adjacent windows/0.4 Mb. Immediate flanks of the primary interval have mean recent-set weight 0.594 and mixed recent/ancient dominance; other LG5 windows average 0.077 recent, 0.441 ancient and 0.363 species-history weight. This is a strong, sustained local topology-weight shift, with some support extending beyond the approximate candidate boundary. It remains evidence of genealogy support, not proof that the full 0.4 Mb is introgressed ancestry.

## Comparison with Neurospora

At common 100-kb physical resolution, monkeyflower q2 recent support has 131 >0.5 runs (median 0.1 Mb, maximum 0.6 Mb); q3 ancient has 401 runs (median 0.1 Mb, maximum 1.4 Mb). The existing Neurospora LGI 100-kb block series has two q3 >0.5 runs totaling 6.9 Mb, median 3.45 Mb and maximum 4.5 Mb. Its 100-kb autocorrelation is 0.950 and transition rate is 0.43/Mb, compared with monkeyflower q3's 0.139 and 4.11/Mb. Monkeyflower's recent set has stronger short-distance autocorrelation (0.412) but still much shorter runs than Neurospora. These results support a more fragmented topology-weight mosaic in monkeyflowers and longer, more coherent q3 support in Neurospora.

The Neurospora series consists of TWISST weights summarized from 50-variable-site trees on one linkage group, while monkeyflower uses 100-kb TWISST weights over 10 chromosomes and a different set of taxa. Pairwise autocorrelation estimates are descriptive; 100-kb bins harmonize physical resolution, not taxon sampling, inference pipeline, or biological meaning. The comparison does not establish that introgression rather than recombination suppression causes these differences.

## Answers and limitations

1. **Windows:** 1,948 informative windows across 10 scaffolds.
2. **Coordinate mapping:** validated by unique genomic keys across the coordinate-bearing files; count-only weights match the joined export row-by-row but have no keys themselves.
3. **Fragmented or sustained:** a mixture; many single-window runs plus sustained q2 runs up to 0.6 Mb and q3 up to 1.4 Mb.
4. **Longest alternative runs:** q2 0.6 Mb (LG3 and LG10); q3 1.4 Mb (LG5).
5. **Chromosome 5:** clear recent-set topology replacement over the reported windows, with >0.5 support extending to 0.4 Mb at this resolution.
6. **Neurospora:** its q3 support is markedly longer and more spatially coherent (4.5-Mb longest run, 0.95 lag-100-kb autocorrelation) than monkeyflower's separate introgression-associated topology tracks.
7. **Genome-wide support:** removing 200–300 kb alters mean weights by at most about 0.13 percentage points; the ranking does not change.
8. **Inclusion:** conditional GO as a topology-only spatial example. Do not present it as branch-length evidence or claim that alternative topologies uniquely identify introgressed tracts.

No local Newick trees were found in the monkeyflower analysis/data directories. Therefore branch-length statistics A(w), B(w), quartet internal branches, CU estimates and local-tree species-tree inference are not available from these inputs.
