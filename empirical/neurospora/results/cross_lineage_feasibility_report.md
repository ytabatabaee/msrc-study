# Neurospora cross-lineage genealogy comparison: feasibility audit

## Decision

**Do not start a cross-lineage local-tree analysis yet.** In this checkout there is no coordinate-indexed local Newick or Relate set for L1, L5, L6, or L8. The only identified genome-wide window-tree set is the Martin & Van Belleghem L10 dataset already represented by the Stage 0/1 outputs. The requested L1/L5/L6/L8 topology and branch-length tracks cannot be made from the current tree files.

A reconstruction path exists without raw-read processing: Corcoran et al.'s Dryad package contains lineage VCFs and masked consensus sequences mapped to the 2509a coordinate system, plus an outgroup alignment. But the relevant lineage-specific VCF files are hundreds of MB each and the shared consensus archive is about 993 MB. Building new alignments, inferring local trees, and validating outgroup/sample coverage would be a substantial new analysis. No such files were downloaded in this audit, and no trees were reconstructed.

## What is on disk and what was checked

`data/neurospora/raw/` currently contains only `.gitignore`; the ignored `data/neurospora/raw/neurospora.tar.gz` and extracted `data/neurospora/raw/neurospora/` tree files are absent in this checkout. The tracked `source_manifest.tsv` and Stage 0 report preserve the prior archive identity and contents:

- Zenodo `neurospora.tar.gz`, 15,577,951 bytes, MD5 `571324702709eb1c2732a3fd993edc93`.
- Expected local Newick file: `data/neurospora/raw/neurospora/trees/cras_sit_hist_tet_L10.min20MAC1.phyml_GTR_bionj_w50sites.trees.gz`.
- Matching window-coordinate file: `data/neurospora/raw/neurospora/trees/cras_sit_hist_tet_L10.min20MAC1.phyml_GTR_bionj_w50sites.data.tsv.gz`.
- Published L10 topology weights: `data/neurospora/raw/neurospora/cras_his_L10.min20MAC1.phyml_GTR_bionj_w50sites.weights.tsv.gz`; the archive also includes `cra_his_L10.topos` and L10 genotype data.

The Stage 0 audit records 45,303 window records, of which 32,472 are Newick trees and 12,831 are `NA`; each actual tree has 26 tips and branch lengths. Windows contain 50 variable sites, with variable physical widths, across LGI–LGVII. LGI records cover 752–9,464,227 bp. TWISST uses 10 L10 mat-A and 12 L10 mat-a samples plus *N. crassa* and *N. hispaniola* group representatives; additional tree tips include *N. sitophila* and an L6 reference sample. This is a single L10 tree set, not independent L1/L5/L6/L8 data. The detailed inventory also preserves that individual CJ mating-type labels are not embedded in the tree tips.

The Zenodo record describes one `neurospora.tar.gz` file containing genotypes, trees, and weights. Martin & Van Belleghem state that their fungal local-tree analysis used the 22 L10 samples (10 mat A, 12 mat a), chosen because of known introgression from *N. hispaniola*. The article's tree analyses use 50-SNP windows. The public Zenodo file endpoint returned HTTP 403 during this audit, so the ignored archive could not be restored to the current checkout. The committed Stage 0 outputs and source manifest therefore serve as the archived local inspection record; no claim is made that new L1/L5/L6/L8 files were hidden in the unavailable tarball.

## Lineage histories relevant to a natural comparison

The companion machine-readable `lineage_mechanism_inventory.tsv` records exact interval endpoints and source coordinates. Corcoran et al. Table 1 puts the suppressed-region intervals on the *N. tetrasperma* 2509a reference: L1 1.15–8.29 Mb (7.14 Mb; 75%), L5 0.74–7.98 Mb (7.24 Mb; 76%), L6 0.68–8.60 Mb (7.92 Mb; 84%), L8 1.54–7.54 Mb (6.0 Mb; 63%), and L10 1.38–7.58 Mb (6.2 Mb; 66%). Sun et al. (2017) independently report suppression in all studied *N. tetrasperma* lineages, about 8 Mb / 86% overall, but stress lineage-specific extents.

- **L1:** L1 has extensive suppression and Sun et al. (2017) describes both mating-type chromosomes as collinear. It is not a clean no-introgression control: Sun et al. (2012) inferred an approximately 4.1-Mb mat-a introgression from an unidentified heterothallic taxon, while Corcoran et al. (2016) did not confirm mat introgression in L1. Treat introgression status as conflicting across studies, not excluded. The 2012 work includes only a small number of locus genealogies (for example, `eth-1` and `cys-9` in the central region), not a dense genome-wide local-tree track.
- **L5/L6:** Corcoran et al. (2016) found introgression associated with mat A in both lineages but could not confidently identify the donor. Sun et al. (2017) report a 6.73-Mb inversion in L5 mat-A chromosomes (732,329–7,462,114 bp) and the corresponding inversion-breakpoint pattern in the L6 reference genome. Sun et al. interpret the lineage-specific inversions as derived after the shared recombination-suppression state was established; the exact timing relationship of this L5/L6 inversion pattern to each lineage's suppression expansion remains unresolved. Corcoran's six-gene concatenated SR tree is not a spatial genealogy series.
- **L8:** Suppression spans about 6 Mb. Sun et al. (2017) report a 255,806-bp inversion in sampled L8 mat-a chromosomes (4,805,724–5,061,530 bp). Corcoran et al. (2016) did not confirm mat-chromosome introgression in L8; that is a non-detection, not proof of lifetime absence. No coordinate-indexed local trees were found.
- **L10:** Corcoran et al. support *N. hispaniola* introgression into mat A. Sun et al. report a 3,632,717-bp inversion only in strain CJ01 mat A (4,332,431–7,965,148 bp), not fixed across L10. The existing local tree dataset contains 32,472 usable trees. It remains the only line with the necessary dense window genealogies found here.

A prior study (Sun et al. 2012) also proposed mat-a introgression in L4 and L9, and Corcoran et al. (2016) later found other lineage/donor-specific patterns. This history cautions against treating labels such as “not confirmed” as proof of absence or defining a single introgression topology as the universal alternative across lineages.

## Are cross-lineage spatial comparisons technically possible?

The main data source is Corcoran et al. (2016), Dryad DOI `10.5061/dryad.162mh`, separate from the Martin & Van Belleghem Zenodo mirror DOI `10.5061/dryad.4jm83` (the latter record is hosted on Zenodo). The Corcoran package has these relevant sequence inputs:

- `L1_vcf.tar.gz` (595.36 MB).
- `L5_vcf_scaffold_1_to_2.tar.gz` (603.74 MB) and `L5_vcf_scaffold_3_to_7.tar.gz` (844.79 MB).
- `L8_vcf_scaffold_1_to_2.tar.gz` (554.89 MB) and `L8_vcf_scaffold_3_to_7.tar.gz` (781.29 MB).
- `Ntet_consensus_seqeunces.zip` (993.21 MB), containing masked consensus sequences for 92 strains on the seven largest scaffolds, with coordinates following the 2509 reference.
- `outgroup_alignment.zip` (63.38 MB), chromosome-wise alignments for heterothallic species with coordinates on the 2509 reference.
- `Fig5_alignments.zip` (28.11 KB), the six SR and two PAR gene alignments used for the published concatenated phylogenies.

The Dryad package therefore offers a path to reconstruct window trees from already-derived variant/consensus sequences without reprocessing raw reads. For L1, L5, and L8, a careful lineage sample list, mat-A/mat-a assignment, outgroup alignment, missing-data, and local-coordinate coverage audit would be required before tree inference. The page does not advertise any L1/L5/L8 per-window Newick or Relate genealogies. The listed VCFs are per lineage, but the described consensus and outgroup alignments are the clearest route to extracting same-coordinate sequence windows. L6 has no separately named VCF tarball on the Dryad listing; availability of suitable L6 sample sequences and matching outgroups in the shared consensus archive must be verified before treating L6 as reconstructable at population scale.

The relevant Corcoran suppression, inversion, and L10 local-tree axes are all expressed against the 2509a reference in the cited source metadata, so coordinate comparisons can in principle share that reference. However, any newly built alignment must retain those reference coordinates and confirm strand/order with the source sequence and scaffold metadata. Do not join positions just because both are labeled LGI if a future source uses another assembly.

## Biological answer at this stage

The published natural contrasts are promising but confounded: L1 is collinear yet has a prior, disputed mat-a introgression signal; L5/L6 combine introgression of uncertain donor with a large mat-A inversion; L8 has suppression and a small inversion but no introgression confirmed by Corcoran et al.; L10 combines known *N. hispaniola* introgression, a strain-specific inversion, and strong local genealogical shift. This is sufficient to motivate a cross-lineage test, but it cannot currently determine whether suppression alone produces MSRC-like topology or branch-length patterns. There are no suitable local trees for the key contrasts, and lineage history, donor, rearrangement state, and sampling differ together.

Thus the requested topology frequencies, explicitly defined q1/q2/q3 and D, manuscript-consistent A(w), B(w), internal branches, spatial autocorrelation, inside/outside contrasts, and tree-removal sensitivity estimates are **not calculated** for additional lineages. No cross-lineage spatial figure is made. In particular, no introgression-associated topology is imposed on L1/L5/L6/L8.

## Recommended next step

If the scientific priority justifies a new analysis, retrieve only the needed lineage-consensus sequences and the outgroup alignment from Dryad, first confirming exact tips and assembly coordinates. Build a narrowly scoped, pre-registered local-tree dataset for one contrast at a time (L8 is a reasonable initial candidate for “no introgression confirmed” but not “no introgression”; L1 should be treated as disputed). Only after successful coordinate/tip validation should we implement comparable quartet and branch-length calculations. The missing local Newick files are not recoverable from existing committed summary tables; they require tree inference from sequence inputs.

## Primary sources

- Corcoran et al. 2016, *Genome Research*: https://pmc.ncbi.nlm.nih.gov/articles/PMC4817772/; Dryad data: https://datadryad.org/dataset/doi%3A10.5061/dryad.162mh
- Martin & Van Belleghem 2017, *Genetics*: https://pmc.ncbi.nlm.nih.gov/articles/PMC5419486/; data record: https://zenodo.org/records/4946041
- Sun et al. 2017, *Nature Communications*: https://pmc.ncbi.nlm.nih.gov/articles/PMC5658415/
- Sun et al. 2012, *PLOS Genetics*: https://pmc.ncbi.nlm.nih.gov/articles/PMC3406010/
- Corcoran et al. 2014 lineage/species background: cited within Corcoran et al. 2016; details are not needed to determine tree-file availability.
