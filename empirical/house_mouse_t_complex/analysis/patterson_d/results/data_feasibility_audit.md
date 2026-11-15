# House mouse nucleotide Patterson D feasibility audit

**Status: READY_TO_RUN.** Source archive inspection found no nucleotide genotypes: nested `ML_trees.zip` contains 4,046 `.contree` trees, and supplementary variation bundles contain position-only coordinates/annotations, not alleles.

The public Harr et al. AllMouse PASS VCF is bgzip/tabix indexed and its HTTP byte-range access was verified. Header inspection returned 67 samples and chr17 length 94,987,271; indexed regional fetch succeeded. The 55 sample IDs used by the frozen designs are present. No whole VCF download was made.

Sample counts recovered from the tree tip union and `stage2_t_only_sample_coverage.tsv`: 15 t-designated tree samples (Mmd 4, Mmm 8, Mmc 3), 32 standard mice, and 8 *M. spretus*. The VCF includes additional samples; the analysis uses the exact 55 tree tips. CZE6/CR29 appears in the public cohort but is absent from the tree-tip set and was excluded; the source paper flags it as a possible partial t carrier. Sample arrangement labels come from archived `_tHaplSubset` records, not topology. Geographic population codes are recorded in the sample table.

The archived source README documents the tree coordinate offset of 5 Mb. The VCF header references mm10, POS is one-based, and the validated inversion intervals are GRCm38/mm10. The nucleotide scan covers the full VCF chr17 contig (POS 1-94,987,271). Inversion/collinear/flank contrasts focus on chr17:5,000,000-42,050,000, covering all four validated intervals and a 2 Mb distal flank; exact tree overlap ends at 37,055,000.

The three designs were frozen in `../config/frozen_quartets.tsv`, and exact sample-to-role assignments were fixed in the runner before the first D scan. `../config/frozen_quartet_members.tsv` is the persisted listing of those pre-scan assignments. Each is rooted with *M. spretus*, a separate species used as a close outgroup in the source study. They are exploratory carrier-versus-standard population contrasts; they are not conventional species-tree null quartets, and unphased diploid genotypes do not isolate the t-bearing homolog. Geographic structure, admixture, ancestral polymorphism, and SNP ascertainment remain material alternatives.

Estimator: per-site group ALT allele frequencies; retain sites where all 8 outgroup samples are called and have ALT frequency zero; ABBA=(1-p1)p2p3, BABA=p1(1-p2)p3. SNP-only site counts do not estimate invariant callable bases or mutation rates. Missing groups are retained as missingness counts. Multi-allelic records are excluded.
