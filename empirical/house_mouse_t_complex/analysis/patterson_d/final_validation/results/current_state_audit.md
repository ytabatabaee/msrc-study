# Current-state audit

## Completed and verified

- The original pipeline uses indexed regional access to the Harr et al. AllMouse PASS VCF (GRCm38/mm10 chr17), fixed frozen sample IDs, genotype allele frequencies, pooled ABBA/BABA, and 5–500-kb aggregation. Existing input checksums, sample mapping, coordinate checks, and synthetic tests are present.
- The additive polarization extension correctly includes both fixed-REF and fixed-ALT outgroups (fixed-ALT sites complement ALT frequencies), while excluding polymorphic/missing outgroups. Heterozygotes contribute 0.5 ALT frequency; incomplete individuals are omitted; biallelic one-base A/C/G/T sites only. Independent scalar tests now cover pure/balanced patterns, recoding, fixed ALT, missing calls, polymorphic outgroup, and multiple individuals.
- Exact sample-ID topology scoring classifies each induced quartet once, resolves ties as unresolved, reports q1/q2/q3 and D_topology=q3−q2, and weights local trees equally. Induced quartet totals are diagnostics, not independent loci. Exact matching is restricted to shared 5-kb starts.
- The optional transversion and 80/90/100%-ingroup-call-rate scan now has 65/65 complete chunk checkpoints and regenerated pooled tables. The prior 59 valid chunks were reused; only six missing chunks were fetched. A failed remote BGZF read was recoverable on bounded retry; no full VCF download or imputation occurred.
- Corrected full-chromosome and exact-tree-matched summaries, original/corrected comparisons, block uncertainty, and all source figures remain preserved. The final package provides a new four-panel figure and exact reproduction commands.

## Superseded estimates

The original fixed-REF-only D estimates (D1 −0.13038, D2 +0.02569, D3 +0.05725) are historical. Corrected fixed-REF+ALT full-contig estimates are +0.19955, +0.13186, and +0.08359. The earlier D1 correlation near −0.097 against q_t_alt is a different topology estimand; it is superseded for exact quartet comparison by D_topology=q3−q2 and corrected polarization. Neither result should be conflated with other regions or topology definitions.

## Remaining mathematical/computational concerns

- The estimator is allele-frequency based and unphased. No error was found in allele orientation, fixed-ALT complementation, ABBA/BABA pooling, or zero-denominator handling. The test and coordinate audit cover the critical paths; output tables remain distinct by genomic support.
- The retry endpoint can truncate BGZF at high virtual offsets. The optional filter scan completed via bounded retries. A separate leave-one-out/balanced-subsampling pass did not finish: the endpoint truncated at byte offset 48,567,409,620. No sample-level estimates are included.
- At physical-span level, robustness tables include all nucleotide windows from 5–37.055 Mb, including some without trees. `robustness_exact_tree_matched.tsv` separately restricts to exact shared starts. These are not interchangeable.
- Block uncertainty is assessed at 250/500/1,000 kb, but matched-span contrasts are imprecise. The circular-shift implementation compresses to observed blocks and rotates their order; it does not preserve missing physical intervals or irregular inversion/tree coverage. Treat its p-values as descriptive.

## Biological confounding and data availability

Population-history review supports the broad *musculus*+*castaneus* relationship but also recent divergence, introgression, and population structure. All three are rooted and broadly compatible with major-subspecies relationships because P1/P2 share a subspecies while P3 is another subspecies. They remain exploratory because arrangement-selected carriers and standards are not guaranteed exchangeable; D3 has only three carriers. The fixed *M. spretus* outgroup roots the quartet but does not prove the ingroup null. These judgments use published population studies and sample provenance, not chr17 D values.

The local sources include exact sample IDs/geographies, indexed unphased diploid genotypes, published pseudo-t carrier tips/local trees, and published inversion intervals. No verified sample-matched haplotype phase or chromosome-specific wild-carrier genotypes were found. Therefore no defensible t-specific nucleotide D can be obtained from the current data; no carrier ALT alleles were assigned to the t chromosome.

See `../results/final_validation/quartet_null_audit.md`, `../results/final_validation/haplotype_feasibility.md`, and the final result tables. Current readiness: **CONDITIONAL GO** for a carefully qualified population-level allele-sharing analysis, not direct phased t-haplotype evidence.
