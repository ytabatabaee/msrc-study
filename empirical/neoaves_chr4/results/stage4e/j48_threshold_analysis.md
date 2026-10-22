# Recovered genome-wide QQS threshold analysis

The published `clade-rec.stat.xz` source covers all 63,430 locus IDs and all chromosomes. The recovered values are full-source QQS statistics, not quartet frequencies recomputed from the 48-taxon J48 gene-tree files.

Chromosome-4 q values reproduce the existing Stage 1 table exactly: 18,109 focal rows compared and zero q-value mismatches.

The affected-locus weights below are valid-observation fractions for each published-statistics comparison. They are distinct from physical genome fractions and from the 48-taxon ASTRAL mixture parameter.

N61 changes from q1 in the non-chr4 background to q2 on chr4; its q1=q2 crossover is epsilon=0.1540, above the observed valid-locus weight 0.0985. The PNAS, primary, stringent, and inclusive region comparisons also have interior q1=q2 crossovers, but their weights and estimands are source-statistics quantities.

N62 changes from q2 in the non-chr4 background to q1 on chr4; its q1=q2 crossover is epsilon=0.01354, below the observed valid-locus weight 0.0978. The PNAS comparison has q1=q2 crossover 0.00303 and q2=q3 crossover 0.8453; the third topology must be checked at every pairwise equality.

Columbea remains q1-dominant at both non-chr4 and chr4 endpoints, so no dominant-topology crossover is supported for that comparison.

These thresholds are empirical mixtures of published full-source QQS statistics. They do not predict the 48-taxon ASTRAL transition exactly because taxon sampling, estimand, and locus weighting differ.
