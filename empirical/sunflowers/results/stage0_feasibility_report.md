# Sunflower Stage 0 feasibility report

**Decision: CONDITIONAL GO for a bounded Stage 1 pilot; NO-GO for claiming sunflowers as a demonstrated MSRC empirical success at this stage.** Public local trees exist, use the same assembly label as the independent haploblock study, and have enough spatial indexing to investigate a small number of large regions. A crucial independent arrangement-to-tip link remains unverified, and reported block boundaries are rounded to Mbp. Stage 1A/1B should not be interpreted as an arrangement-associated test until these metadata are resolved.

## Answers to the Stage 0 questions

1. **Are public local Newicks available?** Yes: 3,242 coordinate-indexed 10.2-kb physical-window Newicks plus an archive of the matching `.treefile` files in the Owens et al. GitHub repository.
2. **Do they retain branch lengths and coordinates?** Yes: every record has numeric branch lengths; each tree key provides chromosome, start, and end in bp. The full 11,797-window manifest includes coordinates and called-base counts. Trees have 22 tips each in this set.
3. **Can suppression-associated regions be identified in the same assembly?** The paper states that haploblock positions are relative to HA412-HOv2, and the 2023 paper maps reads to HA412-HOv2. No liftover is needed. But the 37 published intervals are reported to whole-Mbp precision (sometimes as multiple segments); exact boundaries for broad overlay are still needed.
4. **Is there enough inside/outside coverage?** Screening on rounded published boundaries shows retained trees inside multiple large blocks: ann13.01 97, ann05.01 57, ann11.01 38, ann15.01 25. These counts are based on coarse rounded regions and are not sufficient for an exact edge/flank comparison. Tree density is highly nonuniform and must be checked with exact intervals and quality data. The primary study reports 3,267 retained windows; the public tree files contain 3,242. The public manifest has 11,797 names, of which 3,273 have ≥2,000 called bases; 31 qualifying manifest windows have no Newick. These count differences need provenance reconciliation before detailed coverage claims.
5. **Do suppressed regions shift topology/branch lengths?** Not tested. The presence of inversion/haploblock evidence is not itself evidence of a local genealogical shift.
6. **Are shifts exceptional vs genomic controls?** Not tested.
7. **Does region exclusion alter species-tree topology/CU lengths?** Not tested.
8. **How do introgression and alternatives affect interpretation?** Owens et al. report pervasive admixture/introgression and recombination-associated changes in concordance, so these are major alternatives/confounders. Any future effect requires explicit differentiation of recombination suppression, introgression, ILS, linked selection, and tree uncertainty.
9. **Should sunflowers be included as the plant empirical example now?** Not as a successful MSRC example. They are a plausible conditional pilot. No suppression-associated spatial genealogical signal or inference sensitivity has been demonstrated.

## Why the gate is conditional

The study-level assembly and coordinate schemes are compatible. Study A prose reports genetic-map suppression at ten haploblocks; the Extended Data Table 1 transcription flags 11 rows with `g`, an unresolved internal count discrepancy. Structural evidence is also reported at many blocks. Study B uses the same HA412-HOv2 reference, has coordinate-bearing local ML trees with branch lengths, and samples overlapping annual taxa. These facts establish a feasible data route.

However, the core contrast needs to be frozen from arrangement labels independent of the observed tree topology. Extended Data Table 1 gives rounded Mbp intervals and allele frequencies/selected Hi-C genotypes, but the examined Study B repository does not cross-reference its 22 tree tip IDs to their genotypes at those blocks. The paper's 2-MB supplement endpoint returned CAPTCHA HTML in this environment; a separately listed 362-byte H. p. fallax coordinate file was also inaccessible from the file endpoint (HTTP 401). The paper prose/table also differ on the number of map-suppressed blocks (10 vs 11 `g` rows). Until a sample-to-arrangement crosswalk and exact interval boundaries are obtained from a reliable source, defining q2 as an arrangement-associated quartet would be unsupported.

The README and `haploblock_inventory.tsv` preserve the published block data, including suppression/structural evidence, allele frequencies, multi-part intervals, and coarse tree overlap counts. `tree_coordinate_inventory.tsv` contains all 3,242 local tree coordinates and format checks. `study_compatibility.tsv` records the reference, taxon, ploidy, sampling, tree-method, and branch-unit comparison.

## Next steps to unlock Stage 1

1. Obtain exact HA412-HOv2 interval coordinates and any multi-segment mappings for the focal candidate blocks (prefer ann05.01 and ann13.01 because both are large and have independent inversion/recombination evidence).
2. Resolve the 2023 tree sample IDs against the original sample metadata and genotype those exact individuals at candidate haploblocks using an already deposited manageable genotype resource; do not start from raw reads.
3. Recompute exact window counts inside, in adjacent flanks, and in matched same-chromosome controls using called-base counts and tree quality.
4. Then freeze multiple biologically justified quartets, select q1/q2/q3 independently of tree signal, and run the spatial topology/branch-length and block-resampling sensitivity analyses.

No numerical Stage 1 result files or figures were created because those analyses have not been performed.
