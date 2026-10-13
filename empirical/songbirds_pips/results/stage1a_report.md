# Songbirds PIPs Stage 1A: shared_CM027530.1 reproducibility audit

Stage 1A asked whether the published shared-PIP analysis for `CM027530.1` can be reproduced from public intermediate data or author code without raw-read reprocessing. The answer is no: public sources provide the representative sample set, coordinates, and workflow code, but not the BAMs, SAFs, pairwise SFS files, dxy matrices, or Newick/tree objects required to regenerate the published neighbor-joining analysis directly.

## Frozen representative set

- Samples: 22
- Species: 11
- Genotypes: AA=11, BB=11
- Design: one AA and one BB representative per species: `True`
- Selection source: Dataset S5, not a new Stage-1A random draw. The author code contains a random picker, but the published representatives are fixed in Dataset S5.

## Shared interval

- Reference: `GCA_001746935.2` (*Setophaga coronata*)
- Chromosome/scaffold: `CM027530.1`
- Derived shared overlap: `CM027530.1:574805-874807`
- Length: `300002` bp under the source end-minus-start convention used by Dataset S1.
- Rule: `max(PIP.start)` and `min(PIP.end)` across the 11 species-specific `CM027530.1` PIPs in Dataset S1, matching the article Methods statement that shared PIP regions are the largest possible overlap among species-level PIPs.
- All 11 species-specific PIPs are mapped to `GCA_001746935.2` (*Setophaga coronata*).
- Exact author shared-COMP interval was not present in Dataset S1, Dataset S5, or packaged B3 code; B3a expects it in a user-provided `shared_pips_meta` file.

## Reusable intermediates

No reusable public CM027530.1 VCF/BCF, BEAGLE genotype-likelihood file, SAF, pairwise SFS, dxy matrix, BAM/CRAM, or Newick tree was found. Figshare 27284553 is a preprocessing/code repository and contains no CM027530.1 shared-PIP intermediate data. NCBI BioProjects list raw sequence reads.

## Cheapest Stage 1B route

First try to obtain author intermediate BAM/SAF/dxy files or the `shared_pips_meta` table for `CM027530.1`. If those are unavailable, a limited selected-sample regional reprocessing workflow is required: retrieve the 22 samples' raw reads or BAMs if authors can provide them, align/index to `GCA_001746935.2`, run ANGSD SAF over the shared PIP and matched COMP region, compute pairwise winsfs/dxy, and only then reconstruct NJ/local-tree outputs. Do not run this automatically from Stage 1A.
