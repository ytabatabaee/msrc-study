# Independent raw-alignment validation

Deterministic stratified sample seed: 20261118. Unique loci sampled: 46 (inside PNAS_1/2/3, outside high-information, outside low-information, and outside with ambiguity/missing calls). N61-primary informative-count lower/upper quintile cutoffs were 12/23.

The checksum-verified local archive was streamed once; 6,247 chr4 alignments matched the manifest. Alternative quartet per-locus all-site and transversion counts were computed in this pass. Raw validation covered the three original primary quartets and their three predeclared Meleagris-outgroup designs; each independently decoded row was compared to saved all-site and transversion ABBA/BABA integers.

Validation rows PASS: 544; discrepancy/format failures: 0; unavailable quartet/locus rows: 4. Unique loci: 46. Duplicate sequence IDs: 0 loci. Unequal full-alignment lengths: 0 loci.

All sampled loci had equal full-alignment lengths and zero duplicate sequence IDs. Exact taxon identifiers were retained. Four unavailable rows are the N62 *Pterocles gutturalis* absences at two loci, once for each outgroup. The independent reference implementation uses its own FASTA parser and direct four-base pattern comparisons; it does not import or invoke the production estimator. It records duplicate IDs, full-alignment length equality, four-taxon alignment length, ambiguous/gap site count, multistate non-pattern sites, and exact ABBA/BABA counts. Outgroup polarization is the observed O base; only patterns with one derived P3/O state and the exact ABBA or BABA configuration count. Transversions-only excludes A↔G and C↔T differences between P3 and O.

## Interpretation

Exact integer agreement is required for both ABBA and BABA, separately for all informative patterns and transversions-only. The per-row comparison is in `raw_alignment_validation.tsv`; the deterministic sample list and strata are in `raw_validation_sample_manifest.tsv`. No production filters were modified to force agreement.
