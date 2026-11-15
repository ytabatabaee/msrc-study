# Nucleotide data feasibility audit: Atlantic cod inversion regions

**Audit date:** 2026-11-15  
**Checkout:** `/Users/ytabatabaee/Desktop/msrc-paper`, branch `main` tracking `origin/main`. Existing unrelated modified and untracked files were present at audit start and were preserved. No commit or push was made.

## Decision

**CONDITIONAL GO for analysis after genotype access and outgroup validation; NO-GO for calculating/reporting Patterson's D in this run.** The published dataset inventory documents nucleotide VCFs for all four supergenes and outside-supergene cod regions. The current checkout contains no VCF/BCF or SNAPP XML, and direct retrieval from Zenodo failed (`curl`: `Could not resolve host: zenodo.org`). No values were imputed from local population trees.

A rooted D test cannot yet be certified: the source inventory lists candidate external-taxon Dsuite input and the upstream population-order files list Gadmac, Gadoga, and Gadcha, but the actual genotype headers/calls were unavailable for checking their presence, missingness, fixed state, assembly compatibility, or callable overlap with the cod population SNP data. Candidate presence in source labels or trees alone is insufficient. Exact quartet members and sample counts therefore remain unfrozen.

## What is documented

The local Zenodo API record (`data/atlantic_cod/raw/zenodo/zenodo_4560275_record.json`, record 4560275; DOI 10.5281/zenodo.4560275) is an open-access inventory for Matschiner et al. (2022). It lists: (i) SNAPP VCFs for LG01/LG02/LG07/LG12 (0.62–2.11 MB), (ii) Dsuite VCFs for those supergenes (35.3–103.1 MB), (iii) outside-supergene SNAPP and Dsuite VCFs (71.6 MB and 3.36 GB), and (iv) a 1.87-GB `gadus_arctogadus_boreogadus_dsuite.vcf.gz`. The exact names, record MD5 values, sizes, and content URLs are in `nucleotide_source_inventory.tsv`. File existence and checksums in that table describe Zenodo metadata only, not local files. The record metadata does not list `.tbi` or `.csi` indexes.

The paper states that read data for the cod specimens in Supplementary Table 4 are in ENA project PRJEB43149 and that SNP datasets in PED/VCF format plus phylogenetic inputs/outputs are deposited at Zenodo record 4560275. The source repository is `https://github.com/mmatschiner/supergenes`, `cod_phylogenomics`. The local `split_vcf.sh` identifies the original genotype source path as `../res/gatk/{gadMor2|melAeg_in_gadMor2_coords}/{original_gadoga|alternative_gadoga}/concatenated.masked.vcf.gz`; it extracts biallelic SNPs on the frozen LG01/LG02/LG07/LG12 intervals. Its full variants retain missing genotypes, while the default tree-input extraction excludes missing genotypes. The window generator uses 250,000-bp windows, at least 500 and at most 1,000 sites, and a 50-bp minimum site spacing. It labels sample IDs by their first 10 characters; LG12 explicitly excludes `Gadmor_lfc1` and `Gadmor_lfc2`. These are documented source rules, not a locally verified VCF header.

## Local availability and access checks

- The cod raw directory has source scripts/tables, source-data text, the supplementary PDF, and Zenodo record JSON; it has no local `.vcf`, `.vcf.gz`, `.bcf`, `.tbi`, `.csi`, or SNAPP XML input. The relevant Zenodo files are not already downloaded.
- Attempted the smallest inversion-restricted deposited SNAPP VCF first (LG01, 2,107,536 bytes; Zenodo MD5 `ea55dae3ce58ba43f9c438480744c98a`). DNS resolution failed before an HTTP response. Thus HTTP Range support, file transfer, gzip integrity, checksum, tabix index availability, chromosome names, and coordinate values could not be empirically verified. No large genotype archive or whole-genome reads were downloaded.
- No exact sample headers or individual genotypes are available locally. The repository maps the 12 population labels to location/ecotype and independently records arrangement states; it does not provide a verified per-VCF sample membership roster. The audit TSV keeps these population labels and metadata separate from the unverified genotype roster.
- The frozen intervals remain the exact 1-based inclusive cod_phylogenomics coordinates: LG01 9,114,741–26,192,386; LG02 18,489,307–24,048,607; LG07 13,610,591–23,019,113; LG12 638,100–14,327,837. No alternate coordinate definition was used.

## Nucleotide and outgroup assessment

- **VCF/BCF genotypes:** VCFs are publicly documented by the Zenodo record; none is locally accessible in this checkout. The record lists no index files.
- **Population frequencies and nucleotide matrices:** not separately documented as the format; SNAPP VCFs are genotype-bearing VCF inputs by filename, while population-frequency matrices would need to be calculated from individual calls. SNAPP XMLs may encode biallelic SNAPP patterns, but none is locally accessible; whether original nucleotide identity, REF/ALT orientation, individual IDs, and window mapping are retained cannot be inspected.
- **Coordinates/reference:** upstream source uses `gadMor2` and the alternate `melAeg_in_gadMor2_coords`; the frozen window intervals are in linkage-group coordinates. Actual VCF `##reference`, contig names, positions, and compatibility across the cod and external-taxon files were not verifiable.
- **External taxa:** `Gadcha`, `Gadmac`, and `Gadoga` occur in the baseline tree and upstream Dsuite order files. This demonstrates candidate taxon labels in the documented workflow, not that the relevant deposited VCFs contain usable, fixed, callable outgroup genotypes at shared sites. The putative species expansions (G. chalcogrammus, G. macrocephalus, G. ogac) and sample identities require confirmation against VCF headers and source supplementary metadata. No external candidate is approved as a polarization outgroup yet.
- **Reference alleles:** will not be treated as ancestral absent independent polarization. A fourth cod population would require the independent outside-inversion baseline tree and an explicitly justified root; none has been frozen.

## Quartet decision and blockers

No quartet designs or sample membership were frozen. The baseline tree and Stage-2 arrangement metadata can guide predeclared contrasts, but no design can satisfy the user's requirement that all four groups have verified genotype coverage and a biologically justified polarization source until the actual VCF headers/calls and shared-site overlap are inspected. In particular, all three candidate external taxa in the Stage-3 tree are not established as sequenced in the cod population VCFs merely by occurring in the tree or Dsuite order list. A four-cod-population unrooted contrast is not silently relabeled Patterson's D.

The numerical D, inside/outside contrasts, block uncertainty, spatial tests, topology matching, and figures were not run. Existing SNAPP tree topology statistics are not nucleotide ABBA/BABA data and were not substituted.

## Next executable step

Restore DNS/network access, then download only the inversion-restricted VCFs and required shared-site outgroup data from the exact Zenodo content URLs in `nucleotide_source_inventory.tsv`; prefer a demonstrably indexed source or verify Range support and download only the required files. Validate record MD5s and gzip integrity, inspect VCF headers/sample IDs/contigs/reference, reconcile each sample to Supplementary Table 4, evaluate outgroup callability/fixation, and compute overlap across files. Only then freeze quartets and run the Stage-4+ analysis. Do not fetch raw whole-genome reads.

## Sources

- Matschiner et al. (2022), “Supergene origin and maintenance in Atlantic cod,” *Nature Ecology & Evolution*, DOI [10.1038/s41559-022-01661-x](https://doi.org/10.1038/s41559-022-01661-x); [full text and Data Availability statement](https://pmc.ncbi.nlm.nih.gov/articles/PMC8986531/).
- Zenodo record [10.5281/zenodo.4560275](https://doi.org/10.5281/zenodo.4560275), API metadata cached locally as stated above.
- Source code [mmatschiner/supergenes](https://github.com/mmatschiner/supergenes), local cod_phylogenomics source snapshot and checksums in the existing download manifest.
- ENA project [PRJEB43149](https://www.ebi.ac.uk/ena/browser/view/PRJEB43149) (read archive is outside the scope of this feasibility check; no reads downloaded).

## Follow-up: local access restored and Stage 1 completed

The initial DNS failure was specific to the restricted execution sandbox. With repository-scoped network permission, the Zenodo API returned HTTP 200 with normal TLS verification; no resolver or certificate settings were changed. The files listed above were downloaded only from the record's exact content URLs and were validated against its MD5 values and with `gzip -t`. Downloads are kept in `data/atlantic_cod/raw/zenodo/local_downloads/`, excluded locally via `.git/info/exclude`, and are not Git changes. No tabix/CSI index was listed by the record, and no remote indexed-range access was assumed.

The four inversion SNAPP VCFs contain 29 individual IDs, biallelic complete-call SNPs, LG-prefixed coordinates, and `gadMor2` reference metadata; the 390,090-SNP outside-supergene SNAPP VCF has the same 29 IDs and complete calls. The inversion Dsuite VCFs include missing genotypes and one roster difference (`Gadmor_lfo3` replaces SNAPP `Gadmor_lfo2`). Their complete-call position sets are exactly the SNAPP VCF position sets. Header sample prefixes independently map to population IDs. The source `run_dsuite.slurm` explicitly maps Borsai to `Outgroup`; Borsai (*Boreogadus saida*) is callable in the same VCFs, providing fixed homozygous genotypes at a subset of SNPs. Gadoga was tested as a single-individual alternative outgroup. No reference base was treated as ancestral.

The quartet designs were frozen in `config/` before site-pattern D was inspected. All four primary population groups have two individual samples each; the outgroup has one diploid sample. `results/scientific_interpretation.md` and downstream tables now document the empirical results. A precise limitation remains: the 426 published cod-population local trees omit Borsai and all other outgroup tips, so exact four-taxon induced q1/q2/q3 and exact `D_topology=q3-q2` cannot be matched to nucleotide D. The reported rooted three-population contrast is labeled separately.
