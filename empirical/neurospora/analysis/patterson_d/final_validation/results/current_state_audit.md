# Neurospora L10 Patterson's D — current-state audit

## Preserved primary analysis

The original analysis is under `empirical/neurospora/analysis/patterson_d/`. Its primary design is frozen as P1=L10 mat-a (10 verified isolates), P2=L10 mat-A (8 verified isolates), P3=*N. hispaniola* 8817@999, O=*N. discreta* Ndis@999. CJ07, CJ11, CJ12, and CJ18 remain excluded because no independent mating-type assignment was found. The primary estimates and original Figures A–F were not overwritten.

The checksum-verified local inputs remain the Martin & Van Belleghem/Dryad genotype archive (MD5 `571324702709eb1c2732a3fd993edc93`) and the user-supplied Corcoran et al./Dryad Ndis alignment (MD5 `95a3eb66abc7d0c44cda5cb7953b2a4f`). The original coordinate audit matches the 2509@6 alignment row to the genotype matrix at all 2,265,369 SNP positions. No coordinate remapping was repeated.

The production D code computes haploid group allele frequencies, requiring at least 50% calls in each multi-isolate L10 group and a called P3 singleton. It retains only biallelic ingroup sites for which the unambiguous outgroup matches one of the two alleles, then pools frequency-weighted ABBA/BABA contributions before calculating D. Fixed-ALT outgroup states are correctly polarized. The implementation has synthetic tests for orientation, unequal group sizes, missing calls, biallelic/multiallelic filtering, and gaps/ambiguities.

## Existing validated computations

- Primary LGI estimates: D_inside=+0.87123, D_outside=+0.44595, ΔD=+0.42527; transversion estimates are +0.87087 and +0.47193.
- The 250-kb immediate flanks are elevated (+0.71914 left; +0.78686 right).
- The frozen physical-block jackknife has 26/13/7 inside blocks and 12/6/4 outside blocks at 250 kb/500 kb/1 Mb. Intervals are wide for outside LGI, especially at 1 Mb.
- Existing LGI q3-versus-D correlation is +0.949 across 100-kb bins, +0.454 inside, and +0.843 outside. q3 is a published TWISST weight, not a rooted quartet topology D.
- The published local trees omit Ndis; no exact four-taxon `D_topology=q3-q2` is available from them.
- The original circular-shift output reports 0.010526. Audit shows 95 total rotations including the observed zero shift (94 nonidentity shifts), not 95 nonidentity rotations.

## New audit findings and discrepancies

The matrix has 2,265,369 rows. The deposited local-tree coordinate table has 45,303 records, each spanning exactly 50 genotype rows (2,265,150 rows assigned); 219 genotype rows fall outside all recorded tree intervals. The original D pipeline scans the full genotype matrix, so its summaries and physical bins include these 219 rows, while `window_patterson_d.tsv` covers only rows assigned to source windows. The 219 rows contribute no ABBA/BABA information on LGI. Across all linkage groups, only 8 informative patterns occur among them (6 on LGIV and 2 on LGV), and removing them leaves each chromosome-level D unchanged to the displayed precision. This is a small accounting distinction, not an estimator error.

The variant matrix is ascertained and has no invariant/callable-base denominator. The local archive and code do not document the complete upstream SNP discovery, genotyping, or missingness filters. The source label includes `min20MAC1`; its full criteria cannot be reliably decoded from the filename alone. The published topology analysis uses nonoverlapping 50-variable-SNP windows and requires each tree sample to be called at ≥40/50 sites/window. Those tree filters do not turn the deposited matrix into a random or genome-wide callable-site panel.

The ordinary species-tree null is biologically qualified: the published methods describe L10 mat-A and mat-a as one recombining population outside the nonrecombining LGI tract. Thus the four-tip rooted pattern is calculable, but a mat-A/mat-a “sister populations” model with expected D=0 is not automatically a valid genome-wide demographic null. Prior published introgression between L10 and *N. hispaniola* is a competing/contributing explanation that must be acknowledged.

## New analyses

See `../README.md` for reproduction and `final_scientific_interpretation.md` for the decision. New validation outputs include raw-site scalar checks, full source-window ascertainment/callability metrics, common-support filters, spatial flank estimates, isolate leave-one-out and balanced subsets, refreshed block/circular-shift audit, and a TWISST genealogy comparison with physical-block uncertainty. Results remain conditional on the deposited variable-site matrix.
