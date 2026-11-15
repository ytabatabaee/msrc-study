# Stage 1 — Quantifying the Neurospora genealogy shift

## Findings

**The topology shift is large, spatially coherent, and linked to the published L10 mat-A introgression. It changes genome-wide quartet frequencies, but it does not change the dominant genome-wide quartet topology. It also does not provide evidence that the species-tree topology estimate is wrong.** The dataset supports an introgressed alternative genealogy preserved in a recombination-suppressed region; it does not uniquely identify MSRC or establish that inclusion of the region biases a species-tree estimator.

The archive contains 45,303 coordinate records; 32,472 have both a Newick tree and valid TWISST counts. Every actual tree has the same 26-tip set, and no non-root branch lengths are missing. The local trees are unrooted BIONJ estimates with substitution-distance branch lengths. The analysis uses four-lineage quartet weights, not a 26-tip species tree: the CJ tips are within-species strain/mating-type samples, and collapsing them as 22 independent species would be biologically invalid.

## Published topology shift

The deposited topology definitions make the resolutions explicit: q1 groups the two L10 mating-type backgrounds (putative species-associated topology); q3 groups L10 mat A with *N. hispaniola*; q2 is the remaining quartet resolution. Counts were normalized to frequencies because the source table stores subtree counts (complete rows sum to 120).

| Region | Valid windows | q1 | q2 | q3 | D = q3 − q1 | Dominant |
|---|---:|---:|---:|---:|---:|---|
| All linkage groups | 32,472 | 0.621 | 0.123 | 0.256 | −0.365 | q1 |
| LGI, published ~1.5–8.0 Mb interval | 5,767 | 0.107 | 0.102 | 0.791 | +0.684 | q3 |
| Outside interval, all other windows | 26,705 | 0.732 | 0.128 | 0.140 | −0.592 | q1 |
| Other LGI outside interval | 2,353 | 0.589 | 0.153 | 0.259 | −0.330 | q1 |
| LGII–LGVII | 24,352 | 0.746 | 0.126 | 0.129 | −0.617 | q1 |
| 250-kb LGI flanks combined | 479 | 0.180 | 0.157 | 0.663 | +0.483 | q3 |

Within the published interval q3 has 79.1% mean weight; outside it q1 has 73.2%. Genomic-span weighting gives nearly identical values (inside q3 0.793; outside q1 0.735). The q3 weight rises sharply near 1.2 Mb and falls near 8.2 Mb in the spatial series. On 100-kb block means, q3 is dominant in two long runs, 1.2–5.7 Mb and 5.8–8.2 Mb, separated by one non-q3-dominant block (69 blocks, 6.9 Mb combined). The signal therefore occupies approximately 7 Mb, with fuzzy edges extending about 0.3 Mb beyond the nominal 1.5–8.0 Mb reference interval. The average q3 weight in the two 250-kb flanks is also elevated; the reported boundary should not be interpreted as an exact topology breakpoint.

Adjacent 100-kb block D values have autocorrelation r = 0.954 (92 adjacent valid block pairs). The difference in mean D between the published interval and the rest of LGI is +1.014; a 100-kb block bootstrap 95% interval is 0.883–1.202. One randomly selected informative window per 100-kb bin preserves the direction and magnitude: inside q3 0.726 versus outside-LGI q1 0.689 (64 versus 30 bins). Equal-count sampling is also directionally stable. The raw 50-site windows are physically nonoverlapping, but these high spatial correlations mean they are not independent evolutionary replicates.

## Branch lengths

The focal quartet branch analysis uses two independently documented matched L10 pairs (CJ03 mat A / CJ04 mat a; CJ15 mat A / CJ16 mat a), plus *N. crassa* and *N. hispaniola*. q1/q2/q3 are assigned from the minimum of the three four-tip pair-sum scores. This provides a direct local-tree check of the deposited group-weight result.

For CJ03/CJ04, q3 is the resolved local quartet in 79.2% of available trees inside the interval, compared with 24.5% outside the interval on LGI and 7.4% on LGII–VII. CJ15/CJ16 shows the same pattern in the full per-window table. The cross-lineage contrast A = d(mat A, *N. hispaniola*) − d(mat a, *N. hispaniola*) is −0.379 substitutions/site inside versus −0.006 outside for CJ03/CJ04; the independent pair gives −0.378 inside versus +0.029 outside. This is consistent with the introgressed mat-A lineage being closer to *N. hispaniola*.

Mean six-pair patristic distance B is 0.824 inside versus 0.762 outside for CJ03/CJ04 (0.823 versus 0.772 for CJ15/CJ16). The mean focal quartet internal edge is 0.248 versus 0.212 substitutions/site (0.243 versus 0.204 in the second pair). These local branch-length summaries are descriptive: they are estimated from 50-variable-site trees, include rate and estimation effects, and are not coalescent units. The directional A contrast replicates strongly across two matched pairs; B and internal-edge changes are smaller and distributions overlap broadly.

## Whole-genome inference sensitivity

Including the region changes weighted quartet concordance substantially: compared with outside-region windows, pooled q1 falls from 0.732 to 0.621 and q3 rises from 0.140 to 0.256. But q1 remains the dominant resolution across all available windows. Thus this dataset does **not** show a flip in the overall genome-wide quartet topology when the introgressed region is included. The interval is 17.8% of informative windows, so it noticeably shifts concordance frequencies without overcoming the genome-wide q1 majority.

The outside-region quartet frequencies have nearly equal discordant resolutions (q2 0.128, q3 0.140) and yield an MSC-equivalent internal length of 0.911 CU under the standard four-taxon MSC formula. The other chromosomes alone give 0.964 CU. For all windows, q2 and q3 are strongly asymmetric (0.123 versus 0.256), and inside the interval they are even more asymmetric (0.102 versus 0.791); a single-tree MSC quartet formula is therefore not valid for those mixtures. Consequently, topology remains q1 overall, while a coherent alternative genealogy changes concordance and makes a pooled CU branch length non-identifiable under the single-tree MSC model. We do not claim a species-tree CU estimate changed from one meaningful value to another.

No ASTRAL estimate is reported. The dataset is designed for TWISST on four lineage groups, has two singleton heterothallic representatives, and its dense windows are highly linked; a many-tip ASTRAL tree would treat sampled within-species haplotypes as species or require an unsupported collapse/independence model. The four-lineage quartet comparison is the compatible estimator supported by these data. It is a lineage-history comparison, not a replacement for the independent autosomal species history.

## Interpretation and manuscript decision

Corcoran et al. independently inferred introgression into L10 mat-A chromosomes using sequence-divergence and population-genomic evidence; Martin & Van Belleghem then showed the same *N. hispaniola* + mat-A topology dominating central LGI. The windowed trees reproduce that strong shift and show altered cross-lineage patristic distances. Independent work documents multiple inversions and large-scale recombination suppression on the mating chromosome; it also concludes structural heterozygosity alone does not explain the full extent of suppression. The most defensible interpretation is **introgression preserved by recombination suppression**, with lineage-specific processes potentially also contributing. The local trees cannot separate the effects of introgression, ancestral mating-type genealogies, recombination suppression, rearrangements, or selection.

The data support more than ordinary stochastic gene-tree variation because the alternative genealogy is independently supported as introgression and is concentrated in a long, spatially coherent region. However, whether that is a *bias* depends on the target: it biases a naïve single-history summary of genomic genealogies, but it is also a real part of the genome's evolutionary history. It does not change the dominant four-lineage topology, and this dataset cannot establish an objectively wrong species-tree estimate. The evidence resembles the MSRC prediction of a large, spatially coherent alternative-genealogy regime in its quartet and branch-distance pattern, but known introgression plus recombination suppression is a direct explanation and is not an MSRC-specific result. The Neoaves case is not directly comparable mechanistically; both demonstrate that a localized regime can alter genome-wide quartet summaries without necessarily changing the most frequent topology.

### Recommendation: CONDITIONAL GO

Include Neurospora as an empirical example of a large introgressed, recombination-suppressed region that changes genome-wide quartet frequencies and makes a one-tree CU interpretation inappropriate. State explicitly that q1 remains the pooled dominant topology, that branch-length summaries are substitution units rather than CU, and that introgression and recombination suppression jointly explain the signal. Do not claim this dataset demonstrates species-tree topology bias, a changed identifiable genome-wide CU branch length, or unique support for MSRC.

## Reproducibility

Run `python3 empirical/neurospora/scripts/run_stage1.py` and `python3 -m unittest empirical.neurospora.tests.test_stage1 -v`. Output tables and four figure sets (PNG, PDF, SVG) are in this directory's sibling `results/` and `figures/`. The raw archive is checksummed in `source_manifest.tsv` and excluded from Git.

Sources: [Martin & Van Belleghem 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5419486/), [Corcoran et al. 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4817772/), [large-scale recombination-suppression and rearrangement study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5658415/), [Zenodo data record](https://zenodo.org/records/4946041).
