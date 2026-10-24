# J48 mixture analysis

All 13 mixture treatments completed with exit code 0 and passed taxon-set, annotation, CU/SU topology, provenance, and output-checksum validation. Each tree has 48 taxa and 45 internal splits.

The natural chr4 fraction in the stage-4d locus manifest is 6262/63430 = 0.098723. This is a locus-manifest fraction and is separate from the with-replacement mixture weight epsilon.

The inferred topology is distinct from both references at epsilon = 0.02 and 0.05. The complete J48_T0 topology first appears at epsilon = 0.10 and persists at all larger tested weights. This is an observed coarse-grid interval, not an estimated critical threshold.

The six reference-specific splits are tracked separately in `j48_mixture_split_presence.tsv`; novel splits are retained in `j48_mixture_topology_summary.tsv`. ASTRAL q1/q2/q3 values are reported only for their inferred branches. They are not connected across incompatible splits and do not provide a fixed-topology Delta-Q statistic.

The QQS crossover estimates (N61 0.1540, N62 0.01354, no Columbea dominant-topology crossover) and this ASTRAL experiment have different estimands. QQS thresholds mix published focal quartet summaries; the mixture experiment reweights resampled gene-tree records before 48-taxon species-tree inference. Their qualitative agreement is limited to showing that different focal relationships can respond at different weights. The observed whole-tree transition cannot validate either QQS epsilon* value.

## Manuscript-ready Results paragraph

To test the sensitivity of species-tree inference to a localized genealogy regime, we analyzed 13 completed 48-taxon ASTRAL-IV runs in which 62,945 gene-tree records were sampled with replacement from chromosome-4 and non-chromosome-4 strata at chr4 weights from 0 to 1. The inferred CU and SU trees had identical unrooted topologies for every mixture. The non-chromosome-4 topology was recovered at epsilon = 0, mixtures at 0.02 and 0.05 were intermediate or distinct from both references, and the full-data J48_T0 topology first appeared at epsilon = 0.10 and persisted through epsilon = 1.00. Thus a sampled chr4 contribution of roughly 10% was sufficient to recover the baseline topology on this coarse, single-resample grid, while individual focal relationships changed at different weights. Because each weight used an independent resample, these results describe observed sensitivity rather than a statistically estimated transition threshold; the mixture weight is a sampling fraction and not a measured fraction of the genome under suppressed recombination.

## Published PNAS comparison

The published PNAS taxon-removal experiment changes taxon sampling and asks how the inferred tree responds to missing taxa. The present experiment keeps the 48-taxon set fixed and changes genomic sampling weights. The reconstructed J48 inputs are not treated as exact copies of the published collapsed-tree inputs. A fixed-topology quartet-score analysis analogous to Delta-Q would require scoring each mixture gene-tree collection against the same predeclared four-taxon alternatives; the stored inferred-tree annotations do not supply that fixed-topology score, so no Delta-Q value is reported.
