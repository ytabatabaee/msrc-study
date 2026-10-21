## Fire ants — analysis complete

Biological question: how the chromosome-16 social-supergene region affects local quartet support, summary-tree topology, and species-tree branch lengths in the Stolle et al. fire-ant dataset.

Dataset size: 213 four-BUSCO windows, 267 individual tree tips, seven recovered TWISST/ASTRAL groups, and 267/267 exact metadata matches from Supplementary Data 1.

Primary Stage-4/5 result: local chromosome-16 genealogy support shifts from the background species-history quartet outside the supergene to the cross-species SB/Sb haplotype quartet inside the independently defined supergene region.

Summary-tree result: background and all-window individual-level ASTRAL4 analyses retain the species focal topology, whereas the supergene-only analysis recovers the haplotype focal topology.

CU/SU result: fixed-topology Stage 6C shows strong focal CULength decreases when supergene windows are added, while SULength responses differ between focal branches.

Focal branch changes:
- focal_richteri_SB_Sb_pair: SU 0.00144347 -> 0.00142827 (-1.05%), CU 2.84094 -> 1.14248 (-59.79%)
- focal_invicta_macdonaghi_SB_Sb_pair: SU 0.000603369 -> 0.00108534 (79.88%), CU 1.68846 -> 0.958651 (-43.22%)

Interpretation: topology can remain stable while branch-length estimates, especially coalescent-unit estimates, are affected by a localized alternative genealogy regime.

Caveat: the source study attributes the shared supergene history to recurrent adaptive introgression. This analysis demonstrates phylogenetic consequences of the localized genealogy regime and does not establish MSRC without gene flow.

Manuscript-ready outputs: `empirical/fire_ants_chr16/results/fire_ants_final_summary.tsv`, `empirical/fire_ants_chr16/results/fire_ants_final_methods.md`, `empirical/fire_ants_chr16/results/fire_ants_final_results.md`, and `empirical/fire_ants_chr16/figures/fire_ants_final_summary.pdf`.

Status: complete / frozen unless manuscript review requires changes.

## Fire-ant chromosome 16 recombination recovery audit

Completed the expanded recombination-source audit in
`empirical/fire_ants_chr16/scripts/07_recombination_validation.py`. Wang et al.
(2013) Supplementary Tables 8–14 yielded 27,934 marker rows across seven
families, with `Si_gnF` scaffold coordinates, `lg1`–`lg16` linkage-group
labels, and cM positions. The prior zero-marker chr16 result was an assembly/
correspondence failure: the tables have no marker sequences, the linkage-group
labels cannot be equated to chromosome numbers, and the upstream
`linkage_map_supergene.txt` correspondence to the Stolle `Si_gnGA`/`chr16`
axis is absent from both the local snapshot and public upstream repository.

Ross and Shoemaker (2018) were added as a direct pedigree route. Their article
and supplementary files document pairwise recombination fractions among
`C294`, `Gp-9`, and `i_126` from 101 polygyne queen progenies, including rare
SB/Sb exchange within the largest inversion. These are retained separately as
pairwise measurements because no validated physical coordinates or machine-
readable numeric fraction table was recovered for projection to the frozen
axis; they do not support cM/Mb intervals. The Dryad raw workbook is publicly
listed but its direct stream returned HTTP 403 in this environment.

Yan et al. (2020) supplementary files and the PRJNA421367 route were audited;
the published LD/inversion evidence does not provide a numeric LD matrix or
cM/Mb track. The status remains
`QUANTITATIVE_RECOMBINATION_TRACK_NOT_RECOVERED`: zero Wang markers map to
target chr16, zero existing MSRC windows have a supported quantitative rate,
and no recombination figures or association tests were generated. New
marker-level, pairwise, window-integration, source-manifest, and provenance
outputs preserve the missing values explicitly. The frozen Stage 4A–6C
results and manifests are unchanged. Validation: 15 tests passed.

## Fire-ant final targeted recombination recovery

The final recovery attempt followed only the Wang (2013) -> Pracana (2017) ->
Stolle assembly-correction -> Ross (2018) lineage. The Stolle 2019 primary
data index yielded `S.inv.AGP.txt`, which orders and orients 31 accessioned
`GCF_000188075.1`/Si_gnG scaffold components on an optical `chr16`
pseudomolecule. This is a validated historical scaffold-to-optical-axis
correspondence and is preserved in
`data/fire_ants_chr16/processed/final_historical_chr16_scaffold_mapping.tsv`.

It does not resolve the required marker transfer: Wang's 27,934 rows use
`Si_gnF.scaffold..._nt...` identifiers and contain no marker sequences, while
the recovered AGP uses `NW_...` GCF accessioned scaffolds. No validated
Si_gnF-to-Si_gnG/Si_gnGA correspondence or sequence alignment was available.
Pracana's paper documents ALLMAPS and seven Wang maps but the cited public
supplement/data record did not expose the ALLMAPS input/order files. Ross and
Shoemaker's 2018 source supplies a 15-locus pedigree panel and pairwise
supergene measurements, not a 2,796-marker RAD map; those pairwise values lack
validated physical intervals on the frozen axis.

The final stage produces zero newly mapped Wang markers, zero cM/Mb intervals,
zero supported chromosome-16 windows out of 96, no recombination associations,
and no recombination figures. Stage 4A–6C outputs and checksums remain frozen.
The classification is
`QUANTITATIVE_RECOMBINATION_TRACK_NOT_RECOVERED`. See
`empirical/fire_ants_chr16/results/final_recombination_recovery_report.md`,
the source manifest, and `scripts/08_final_recombination_recovery.py`.

## Fire-ant chromosome 16 Stage 6C SU branch-length sensitivity

Stage 6C complete — mapped individual-level ASTRAL4/CASTLES-II SU branch-length sensitivity. The exact 267-tip mapping was recovered from Stolle et al. 2022 Supplementary Data 1 without guessing. Original RAxML local-tree branch lengths were preserved, fixed-topology background versus all-window CULength and SULength values were estimated from the same individual-level inputs, and CU/SU changes were compared. Stage 6C supersedes Stage 6B for direct CU-vs-SU comparison because CU and SU are estimated from identical individual-level inputs and the identical fixed background topology. Stage 6B remains an independent grouped/modal-topology CU sensitivity. Stage 5 remains frozen. See `results/stage6c_su_summary.txt`, `results/stage6c_su_branch_length_comparison.tsv`, and `results/stage6c_cu_vs_su_comparison.tsv`.

## Fire-ant chromosome 16 Stage 6 CU branch-length extension

The existing grouped TWISST / ASTRAL4 Stage 6 trees were reused to compare ASTRAL coalescent-unit branch lengths. The primary topology, `T_background` versus `T_all`, remains unchanged, but focal CU lengths were compared branch by branch using canonical unrooted bipartitions. `SULength` was not analyzed because the modal gene-tree inputs do not contain meaningful substitution branch lengths. See `empirical/fire_ants_chr16/results/stage6_cu_branch_length_summary.txt` and `empirical/fire_ants_chr16/results/stage6_cu_branch_length_comparison.tsv`.

# MSRC Project Status

## Controlled quartet Experiment 1: conditional grid using msrc-sim

Completed the first controlled four-taxon quartet experiment under
`simulations/quartet/`. All raw simulations used `ytabatabaee/msrc-sim` at
commit `f5157d2c64068eb32e0eabfe9bec629c28a404f9`. The study uses the exact
finite-state `H_m` calculation and Monte Carlo validation for configuration
`1010`, with durations 0.02, 0.05, 0.10, 0.25, 0.50, 1.00, and 2.00;
symmetric switching rates 0, 0.0025, 0.005, 0.010, 0.025, 0.050, 0.100,
0.250, 0.500, and 1.000; and 100,000 loci per cell.

`Q_ALT = 13|24` dominates `Q_SPECIES = 12|34` in every exact cell on this
grid. The strongest conditional bias is at duration 2.0 and `m=0`, with
exact vector `(0.00610521, 0.98778957, 0.00610521)` and
`Q_ALT-Q_SPECIES = 0.98168436`. The maximum empirical-versus-exact error is
0.0040061; 97.619% of topology comparisons are within 2 SE and all are within
3 SE. The next step is the mechanistic Wright–Fisher/Moran quartet experiment;
this conditional result is not yet a full species-tree inconsistency claim.

## Mechanistic quartet Experiment 2: Wright–Fisher pilot

Completed a validated pilot under `simulations/quartet/mechanistic/` using
`msrc-sim-replicates` at commit `f5157d2c64068eb32e0eabfe9bec629c28a404f9`.
The unconditional pilot contains 72 cells, 12 independent Wright–Fisher
histories per cell, and 20 loci per history; matched ordinary MSC controls use
20,000 loci per cell. The six controls agree with analytic quartet
expectations (maximum absolute z score 1.20). No pilot cell produced a
supported unconditional `Q_ALT > Q_SPECIES` result; all 72 point estimates of
`P(Q_A)-P(Q_S)` were negative. The terminal `1010` pattern occurred in six
cells, with maximum frequency 0.0833, so the strong conditional Experiment 1
bias is rare after mechanistic averaging in this pilot. Figures, processed
tables, audit, provenance, and follow-up plan are saved under the mechanistic
experiment directory. Next step: increase independent histories in targeted
persistence and intermediate-frequency settings while retaining unconditional
sampling and matched MSC controls.

## Mechanistic quartet Experiment 2B: persistence and bias decomposition

Completed the targeted persistence refinement under
`simulations/quartet/mechanistic/`. All 32 planned primary settings and four
short-branch controls finished with 100 independent Wright–Fisher histories per
cell and 20 loci per history. Daughter-branch persistence was measured by
requiring both ancestral branches A and B to end segregating; terminal 2:2
patterns were analyzed separately. Across 3,600 histories, `1010` occurred 24
times and `0101` 17 times. The unconditional decomposition reconstructed direct
`P(Q_A)-P(Q_S)` to within `2.22e-16` in every cell. All 36 unconditional
estimates remained negative; the least negative was `-0.159` in a short-branch
setting. Rare 2:2 categories retained strong conditional alternative support,
but their weights were insufficient for a topology reversal. Suppression
changed species-topology support modestly, and cross fraction 1.0 was shown to
be a mechanistic no-suppression condition rather than an exact ordinary MSC.
History-level bootstrap tables, five diagnostic figures, results, tests, and
the next refinement plan are saved under the mechanistic Experiment 2B
directories. No 500-history adaptive rerun was launched because no cell met
the pre-specified anomaly refinement criterion.

## House mouse t-complex Stage 0/1 scaffold

Implemented a new empirical scaffold for the Kelemen & Vicoso house-mouse
t-haplotype / t-complex local-gene-tree dataset under
`empirical/house_mouse_t_complex/` with the corresponding data root under
`data/house_mouse_t_complex/`.

Stage 0 status: passed against the local
`IST-2017-78-v1+1_Data.zip` archive. The audit script reads the upstream ZIP
directly, inventories the three nested ML tree archives, parses the primary
nested ML-tree ZIP, writes deterministic tree-tip and window inventories, and
maps tips only where supported by the archive README plus Kelemen & Vicoso
2018 paper metadata.

Source paper/data archive: Kelemen & Vicoso public archive
`IST-2017-78-v1+1_Data.zip`; primary target path
`Data/2-Coverage-and_AlleleRatio-Filtered_RAW_SNPs/1-Trees_for_all_5kb_windows/ML_trees.zip`.

Availability of published 5-kb ML gene trees: available. The primary nested
archive contains 4,046 `.fa.contree` ML trees under
`ML_trees/`.

Number of usable trees after audit: 4,046/4,046 valid Newicks, each with 55
tips and branch lengths. Internal labels/support values are present; no
malformed trees, duplicate tip labels, or negative branch lengths were detected.

Mapping status: source-backed strong mapping for all 55 canonical tips. The
mapping uses paper-supported population labels: CAST = M. m. castaneus,
FRA/GER = M. m. domesticus, AFG/CZE/KAZ = M. m. musculus, and SPRE =
M. spretus. Fifteen `_tHaplSubset.fa` tips are treated as pseudo-t/t-haplotype
sequences; the remaining M. musculus tips are standard/noncarrier sequences,
and M. spretus tips are outgroup non-t tips.

ASTRAL-ready input status: produced. Stage 1 wrote the canonical 4,046-tree
file, metadata table, readable mapping TSVs, ASTRAL4 whitespace-safe `.map`
files, and treatment inputs for `T0_STANDARD`, `T1_ALL_WINDOWS`,
`T2_SPATIALLY_THINNED_50kb`, and `T2_SPATIALLY_THINNED_100kb`. Exploratory
local ASTRAL4 diagnostics completed for those four treatments, but these are
Stage-1 diagnostics only and are not a final empirical result.

Current blocker: biological interpretation is gated. The Stage-1 ASTRAL runs
need review of mapping semantics, treatment design, and local-tree dependence
before any result can be stated. These data remain chr17 t-complex windows, not
a genome-wide or outside-inversion background.

Next step: inspect `empirical/house_mouse_t_complex/results/stage0_source_audit.md`,
`empirical/house_mouse_t_complex/results/stage1_report.md`, and
`empirical/house_mouse_t_complex/results/stage1_aster/stage1_astral_summary.tsv`
before deciding whether to promote these exploratory diagnostics into a formal
Stage-2 analysis.

## House mouse t-complex Stage 2

Stage 2B is complete for the local exact controls. The corrected focal
ASTRAL annotations are T0_STANDARD Q_SPECIES with q=(0.369910, 0.343656,
0.286434) and T1_ALL_WINDOWS Q_T_ALT with q=(0.365364, 0.323323, 0.311313).
The arrangement scan shows that SSS and TTT favor Q_SPECIES, while pooled
one-T arrangements favor Q_T_ALT; TTS and STS are the strongest
alternative-supporting patterns. The 1,000-replicate balanced exact analysis
returns Q_T_ALT in 65.8% of mixed replicates, while matched T-only returns
Q_SPECIES in 100%; exact quartet winners agree with ASTRAL in 60/60 validation
replicates.

All three filtering datasets are now compatible after exact deterministic
normalization of 55/55 labels and reproduce STANDARD_ONLY = Q_SPECIES and
ALL_TIPS = Q_T_ALT. All feasible 5-kb phase offsets were evaluated through
500 kb; topology frequencies become more variable and ASTRAL localPP falls
from 1.0 at 5 kb to a median 0.497 in the 500-kb phase sweep. The population
P0/P1 comparison now reports two changed canonical splits and split-specific CU
changes. No source-defined genomic intervals for the archive's concatenated
non-recombined regions were recovered, so `source_defined_blocks_available =
false`.

The remaining limitations are the single chr17 t-complex region, lack of an
independent genome-wide control, and the fact that these are published trees.
The result is an arrangement-state-associated topology sensitivity consistent
with an MSRC-like mechanism, not proof of MSRC or a causal result. See
`empirical/house_mouse_t_complex/results/stage2b_report.md` and
`stage2b_manifest.json`.

The house-mouse spatial figure now includes a published inversion-boundary
panel above the quartet-support track. The four historical intervals use
sequence-validated GRCm38/mm10 coordinates recorded in
`empirical/house_mouse_t_complex/results/t_complex_inversion_intervals.tsv`,
with literature and coordinate provenance in
`t_complex_inversion_provenance.md`. The 2017 archive and structural study
use the same mm10/GRCm38 coordinate system. The resulting descriptive
inside/outside summary is
`t_complex_inversion_interval_summary.tsv`; windows intersecting boundaries
are excluded and linked windows are not treated as independent replicates.

Stage 2B visualization/recombination cleanup is complete. Four-group and
seven-population ASTRAL cladograms, direct q-support tracks, the aligned
q/recombination figure, controls figure, and manuscript `main_v4` figure were
generated and passed programmatic nonempty-figure checks. The primary result
is arrangement-state-conditioned quartet distortion: spatial q means are
heterogeneous and usually near one third, while mixed structural states alter
the aggregate quartet distribution. The six mixed-state classes comprise
64.6% of induced quartets, with STS and TTS driving much of the Q_T_ALT excess.
Standard-only ASTRAL returns Q_SPECIES, all-tips ASTRAL returns Q_T_ALT, and
balanced mixed sampling shifts replicate outcomes toward Q_T_ALT. The source-defined
phylogeny-based recombination-state track directly uses the archived
`Tree_results_dom_ML`, `Tree_results_mus_ML`, and `Tree_results_cas_ML` files;
4,046/4,046 frozen windows match a source classification row for each
subspecies, and all observed source codes have an explicit mapping. This is
not a direct cM/Mb recombination-rate estimate. No coordinate list for the
concatenated non-recombined-region tree was recovered. The next step is
manuscript integration, with no new empirical question started.

## Neoaves chromosome 4 Stage 4E ASTER/ASTRAL4 preparation

## Neoaves chromosome 4 Stage 4F ASTRAL-IV CU/SU extension

Stage 4F code and local validation are complete; full empirical inference is
still pending Bridges-2 execution. The Stage 4E treatment definitions,
predeclared random-matched manifests, and input checksums were preserved.
The runner now supports `topology_cu`, `topology_su`, and `both` modes, writes
separate `tree.cu.nwk` and `tree.su.nwk` outputs, records mode-specific
provenance, and refuses large execution outside SLURM. The collector validates
CU and SU outputs independently, and the analyzer matches internal branches by
normalized unrooted bipartition and terminal branches by exact taxon label.

The official ASTER/ASTRAL-IV interface was verified from its source and
tutorial: `CULength` is coalescent units, while `SULength` is CASTLES
substitution-per-site length; SU estimation also requires an average gene
length and `sulengthtype`. The current configured default is `--genelength
1000 --sulengthtype 1`, which must be checked against source alignment metadata
before interpreting returned SU values. The local four-taxon smoke tests pass
for both modes, including finite nonnegative branch lengths and explicit
annotation parsing. Fixed-topology commands use the documented `-C -c`
interface but remain pending a validated T0 topology.

No Bridges-2 jobs were submitted by this workflow. At the current freeze,
`stage4e_collect_status.tsv` reports all treatment/unit outputs as missing, so
there are no valid 363-species topology, CU, SU, fixed-topology, or
matched-random scientific results to interpret. Pending outputs are represented
as pending/missing rather than zero. Required next commands are the generated
preflight job, then the five primary treatments (`T0`, `T_PNAS`,
`T_STRUCT_PRIMARY`, `T_BLOCK_500K`, `NONCHR4`), followed by collection and
analysis before supplementary treatments or random controls.

## House mouse t-complex recombination-suppression audit

Direct quantitative suppression evidence was recovered from primary literature,
but only as a T–tf marker-pair comparison (7–12 cM conventional versus
0.1–0.5 cM in naturally occurring t/+ heterozygotes) and a region-wide
20-cM / 100–200-fold literature summary. These observations are arrangement-
specific but not coordinate-resolved on the frozen 5-kb window grid, and no
spatial association test with quartet contrast was justified.
The source audit, report, and four-panel schematic are in
empirical/house_mouse_t_complex/results/recombination_suppression_source_audit.tsv,
recombination_suppression_report.md, and
figures/house_mouse_t_complex_suppression.{png,pdf}. The historical
phylogeny-based recombination classification remains supplementary and is not
treated as a direct rate track.

Stage 4C is complete. Historical Stage 4D ASTRAL-MP reproduction attempts are
preserved as provenance but discontinued as a prerequisite for new full-tree
analysis. Stage 4E prepares ASTER/ASTRAL4 treatment inputs and Bridges-2 SLURM
scripts; large inference must run only under SLURM, and Stage 4E is not complete
until cluster outputs are returned, collected, and validated.

## Atlantic cod Stage 7 post-freeze divergence-time sensitivity

Stage 7 is complete as a post-freeze extension to the frozen Atlantic cod analysis. It reuses the published 250-kb SNAPP MCC window trees, frozen inversion coordinates, and frozen population arrangement states to compare the same population pairs inside versus outside each inversion. SNAPP was not rerun, and Stage 5/6 results were not invalidated.

Tree-time audit: 425/426 published MCC trees passed ultrametric and nonnegative-branch validation. The excluded tree is the boundary-overlap window `LG07_023000001_023250000`, which has a materially negative branch length around -0.0032, larger than the documented numerical clamp tolerance; primary fully inside/outside window counts are unchanged. The workflow applies the same root-age prior constraint across window XMLs, but the source files do not explicitly label the Newick branch lengths as years, Ma, or substitutions, so Stage 7 interprets results as within-workflow calibrated SNAPP time units.

Primary Stage-7 result: all four LGs show positive opposite-minus-same arrangement MRCA-time shifts inside inversions under the circular-shift test. The corrected physical-coordinate sensitivity p-values are LG01 0.021739, LG02 0.013333, LG07 0.011236, and LG12 0.089286, with physical-null BH values 0.028986, 0.026667, 0.026667, and 0.089286. LG01, LG02, and LG07 remain robust under this sensitivity check; LG12 is supportive but less robust and sensitive to the spatial null definition. LG07 is topology-time discordant because it has weak topology enrichment but a strong divergence-time shift. Per-window posterior SNAPP samples were not available in the frozen Zenodo/local inventory, so posterior uncertainty was not propagated.

Outputs: `empirical/atlantic_cod/scripts/07_divergence_time_sensitivity.py`, `data/atlantic_cod/processed/stage7_pairwise_mrca_times.tsv`, `data/atlantic_cod/processed/stage7_pairwise_time_shifts.tsv`, `empirical/atlantic_cod/results/stage7_*`, and `empirical/atlantic_cod/figures/atlantic_cod_stage7_*`.

## Atlantic cod Stage 8 post-freeze summary-tree time projection

Stage 8 is complete as a post-freeze extension. It projects aggregate Stage-7 SNAPP pairwise MRCA times onto fixed Stage-6 ASTRAL topologies. ASTRAL4 was not rerun, and ASTRAL branch lengths, CU/SU lengths, support values, and annotations were ignored for the projection.

Primary T0_all result: the outside-only MRCA matrix selected the root edge separating `Gadmor_bat_spc+Gadmor_twc_spc+Gadmor_two_spc` from the other nine populations. Reusing that root, the projected root age increased from 0.0590363 outside-only to 0.0656562 with all non-boundary focal windows, a shift of 0.0066199 or 11.2% in the Stage-7 workflow time units. Internal-node shifts were mixed: 6 nodes older and 5 younger with all windows; median absolute shift was 0.0041335 and maximum absolute shift was 0.012448. Projection fit was good for outside-only inputs and poorer for all-window inputs, consistent with inversion-linked heterogeneity in the aggregate pairwise MRCA matrix. The secondary T1_outside topology selected the same root edge and the same root-age shift, with median absolute internal-node shift 0.0024188.

Outputs: `empirical/atlantic_cod/scripts/08_summary_tree_time_projection.py`, `data/atlantic_cod/processed/stage8_pairwise_mrca_matrices.tsv`, `empirical/atlantic_cod/results/stage8_*`, and `empirical/atlantic_cod/figures/atlantic_cod_stage8_*`.

## Songbirds/PIPs Stage 1B minimal raw-read pilot

Stage 1B was initialized for the frozen shared-PIP candidate
`shared_CM027530.1` on `CM027530.1:574805-874807` using reference
`GCA_001746935.2` and the 22 frozen Stage-1A author representatives.

Design validation passed:

- 22 samples.
- 11 species.
- AA = 11 and BB = 11.
- One AA and one BB representative per species.
- 231 expected pairwise comparisons.
- No sample substitutions, candidate-set edits, interval edits, tree-result
  inspection, or background/control analysis.

Execution is blocked before raw-read download because the local environment
lacks the required executables: `prefetch`, `fasterq-dump`, `bwa`, `samtools`,
`angsd`, and `winsfs`. `Rscript` is present, but the upstream raw-read/SFS
pipeline cannot run. The frozen Stage-1A sample table stores SAMN BioSample
accessions in the `sra_accession` column; an executable rerun must first resolve
those SAMN identifiers to run-level SRR accessions and record the mapping.

No FASTQ/BAM/SAF/SFS/dXY matrix/NJ tree/NJ figure was produced, and the
published AA/BB-across-species signal has not yet been reproduced. Stage 1C is
not justified until Stage 1B succeeds with usable focal-region coverage, stable
pairwise SFS estimation, and a complete or nearly complete dXY matrix.

Records are in `empirical/songbirds_pips/results/stage1b_*`,
`data/songbirds_pips/processed/stage1b_shared_pip_dxy*.tsv`, and
`data/songbirds_pips/intermediate/stage1b_shared_pip/logs/`.

## Stage 4D resumed inference

The full-363 checkpoint was attempted with historical ASTRAL-MP 5.15.1,
four CPU threads, seed 692, annotation level 3, and an 8 GB heap on the
16 GB local host. It was stopped after 429.57 seconds during memory pressure
before input loading completed. Process diagnostics show full-heap garbage
collection; no tree or scientific reproduction result was produced. This is
a computational limitation, not a topology mismatch. Attempt provenance and
diagnostics are in `results/stage4d/attempts/full363_resume_01/`.

The authors' archived final tree matches the supplied `63K.tre` exactly
(RF=0); this is archive consistency, not fresh-inference reproduction. Their
archived full run took 112,915.644 seconds on 28 threads and four P100 GPUs.
The new Jarvis-48 attempt loaded all 62,945 usable trees and entered search;
its completion status is recorded under
`results/stage4d/attempts/J48_ALL_20260919T121120529946Z/`.

Added a provenance-recording checkpoint runner and conservative extraction of
frozen focal splits. The complete Neoaves suite passes 34 tests. No novel
treatments, analysis freeze, or full-tree correction claims are released.
Details and remaining gates are in `results/stage4d_report.md`.

## Neoaves chromosome 4 Stage 4D continuation audit

Continuation reused the existing downloads, 63,430-row locus manifest, and
36-condition Figure-4 summary reproduction; none were regenerated. Existing
FULL363 and J48 ASTRAL attempts have empty output trees and no live processes,
so neither is a completed reproduction checkpoint. Figure-4 source-summary
aggregation exists, but its audit marks de novo scores unverified.

The complete Neoaves test suite returned **29 passed, 1 failed**. The failure
is a pre-existing modification of `results/stage4b_report.md`: `selected` was
changed to `select[..](../..)ed`. The HEAD version exactly matches the frozen
Stage-4B manifest checksum; the working copy does not. No frozen file was
modified during this continuation. Exact expected/current/HEAD checksums are
saved in `empirical/neoaves_chr4/results/stage4d/continuation_audit.json`.

The Stage-4D prior-inventory test passes, meaning that inventory already
contains the modified report; it does not supersede the older Stage-4B freeze.
The user subsequently authorized restoring the single typo. The report now
matches the original Stage-4B frozen SHA256
`0ebdd6697c7754ed53cbaeb9609f4d0c0065e03228a4217e696267f5e04ea03f`.
The Stage-4D prior inventory entry was corrected to that original hash; the
continuation audit retains both hashes as the pre-restoration record.
After restoration and inventory correction, the complete Neoaves suite passes:
**30 passed**. The checksum conflict is resolved.
Novel treatments and the final analysis freeze remain
unexecuted. No new full-tree, random-thinning, structural-correction, or MSRC
manuscript claim is supported by this continuation audit.

## Neoaves chromosome 4 Stage 4A

Stage 4A tests whether frozen rearrangement-associated chr4 neighborhoods
disproportionately contribute support for competing quartet relationships. No
Stage 1-3 outputs were modified, and the Stage 3B strict 2:2 eligibility result
remains zero.

Frozen structural mask rule:

- Structural-associated is the union of +/-250 kb neighborhoods around the
  frozen breakpoint set, clipped to chr4:1-91,310,470.
- The 250 kb radius is the pre-existing Stage 3B primary flank window.
- Stringent, primary, and inclusive breakpoint sets are sensitivity analyses;
  structural-background is each mask's complement.
- The masks are built before and independently of topology data.

Main results:

- N61 q2 is enriched in structural neighborhoods for all masks: fold 1.506,
  1.456, and 1.365 for stringent, primary, and inclusive, with spatial
  circular-shift p-values 0.0031, 0.0166, and 0.0292.
- Columbea and N62 show depletion of q2 in all three masks, so the correlated
  tracks do not support a general chr4-wide enrichment claim.
- N61 q2 falls from 0.454 window-weighted to 0.332 block-normalized, consistent
  with amplification by dense consecutive windows. Columbea and N62 q2 rise
  after block normalization.
- Dominance changes under block normalization for N61 (q2 to q1) and N62 (q1
  to q2), but remains q1 for Columbea.
- Negative q1-minus-q2 event margins occur in 13/16 N61 events, 3/16 Columbea
  events, and 4/16 N62 events. N61 is broadly distributed rather than driven by
  one event, though all event summaries remain association-only.

Outputs are in `empirical/neoaves_chr4/results/stage4a_*`, figures in
`empirical/neoaves_chr4/figures/stage4a_*`, and the deterministic analysis is
implemented by `scripts/04a_structural_topology_support.py`. Six regression
tests cover genealogy-independent masks, run collapsing, count conservation,
deterministic circular shifts, event assignment, and preservation of the Stage
3B zero-eligibility result.

## Neoaves chromosome 4

The existing Neoaves chromosome-4 empirical analysis is treated as frozen for
this Stage-0 Anopheles scaffold. No Neoaves files or conclusions were modified.

## Neoaves chromosome 4 recombination-map audit

The published-source audit is complete in
`empirical/neoaves_chr4/results/recombination_source_audit.tsv` and
`recombination_source_audit.md`. Two pedigree-based chicken sources were
identified: Groenen et al. (2009), with a 9,268-marker map and approximately
1-Mb rate bins, and Elferink et al. (2010), with a 13,340-SNP map and 1,819
approximately 500-kb bins carrying sex-average and sex-specific cM/Mb rates.
Rabie et al. (2004) is retained as a GGA4 radiation-hybrid/BAC correspondence
source, not as a recombination map. Stapley et al. (2010) is retained as a
pedigree-based zebra-finch comparison.

The Elferink supplementary XLS files are publicly listed but the publisher
endpoint returned HTML protection pages during recovery, so no external numeric
file was accepted as downloaded or preserved. The map is on WASHUC2, and no
validated WASHUC2-to-GalGal6 marker conversion was available in the local
snapshot. Therefore no coordinate-resolved cM/Mb table, recombination figure,
inside/outside rate comparison, or significance test was generated. The PNAS
Neoaves data provide a GalGal6 focal coordinate reference and inferred
historical suppression, not direct rates. Existing structural and quartet
outputs remain frozen, and no ASTRAL analysis was rerun. Validation is in
`empirical/neoaves_chr4/tests/test_recombination_audit.sh`.

## Anopheles 2La Stage 0

Implemented a parallel empirical scaffold under `empirical/anopheles_2la/` and
an external data tree under `data/anopheles_2la/`.

Stage-0 outputs now include:

- `config/anopheles_2la.yaml`
- `metadata/region_manifest.tsv`
- `metadata/source_manifest.tsv`
- `metadata/sample_manifest.tsv`
- `processed/candidate_strict_quartets.tsv`
- `results/stage0_data_audit.md`

Datasets identified:

- Dryad Fontaine et al. 2015, DOI `10.5061/dryad.f4114`.
- MalariaGEN Ag3.10 `fontaine-2015-rebuild`, 72 samples.
- MalariaGEN Ag3 metadata and `ag3.karyotype("2La")` API for future karyotype retrieval.

Reference/build status:

- MalariaGEN Ag3 data are recorded as AgamP4-aligned, but the exact 2La
  breakpoint coordinates trace back to PEST/AgamP3/anoGam3-style reporting.
- 2La interval recorded as `2L:20524058-42165532`.
- Exact AgamP4 coordinate equivalence should be confirmed from Fontaine et al.
  table S11, White et al. breakpoint provenance, or an equivalent
  machine-readable source before Stage 1.

Usable sample/karyotype counts from the local Stage-0 sample manifest:

- Sample-level rows: 0.
- `A0_homozygous`: 0.
- `A1_homozygous`: 0.
- `heterokaryotype`: 0.
- `unknown`: 0.
- Candidate strict 2:2 quartets: 0.

Aggregate public availability from MalariaGEN Ag3.10 `fontaine-2015-rebuild`
has also been recorded separately and is not used to emit sample-level quartets:

- Total aggregate samples: 72.
- Species counts: arabiensis 12, coluzzii 11, gambiae 26, melas 4, merus 9,
  quadriannulatus 10.
- Species-level fixed-state counts from non-polymorphic taxa: 14 candidate
  `A0_homozygous` lineage slots and 21 candidate `A1_homozygous` lineage slots;
  37 gambiae/coluzzii samples remain unknown until per-sample karyotypes are
  retrieved.

Blockers:

- Sample-level metadata and 2La karyotype calls were not retrieved locally.
- The attempted inferred Ag3.10 GCS path returned 403, so the exact metadata URL
  layout or `malariagen_data` retrieval path needs resolution.
- No topology or genealogy data should be inspected until sample/karyotype
  predictions are populated and frozen.

Recommended Stage 1:

1. Retrieve small MalariaGEN Ag3 sample metadata and `2La` karyotype tables via
   the official `malariagen_data` API or confirmed GCS object paths.
2. Populate `data/anopheles_2la/metadata/sample_manifest.tsv` with one row per
   biological sample and documented provenance.
3. Rerun Stage 0 to freeze candidate strict 2:2 quartets.
4. Only after freezing those candidates, evaluate local sequence/haplotype/tree
   resources for the 2La interval and matched flanks.

## Anopheles 2La Stage 1A

Stage 1A preserved the Stage-0 zero-row strict-quartet result and added a
metadata-freeze pipeline without accessing haplotypes, SNP calls, local trees,
distance trees, q1/q2/q3, QQS/BQS, or observed topology data.

Coordinate provenance:

- `2L:20524058-42165532` is now recorded as the working AgamP4 interval because
  Ag1000G/AgamP4-based 2La analyses explicitly use these coordinates.
- The earlier Stage-0 concern is retained: the original breakpoint structure was
  reported on PEST/AgamP3 scaffolds, so later coordinate conversions must state
  whether they are 1-based inclusive or 0-based half-open.

API retrieval:

- Installed and used official `malariagen-data==12.0.1` package code.
- Created `empirical/anopheles_2la/scripts/01_fetch_sample_karyotypes.py`.
- Attempted `malariagen_data.Ag3()` with anonymous public access, intending to
  call `sample_sets(release='3.10')`, `sample_metadata(...)`, and
  `karyotype('2La', ...)`.
- Retrieval was blocked before sample-set discovery by a 401 on
  `gs://vo_agam_release/v3-config.json`.

Sample/karyotype counts from authoritative sample-level manifest after the
blocked API attempt:

- Total samples retrieved: 0 of expected 72.
- `A0_homozygous`: 0.
- `A1_homozygous`: 0.
- `heterokaryotype`: 0.
- `unknown`: 0.
- Direct sample-level states: 0.

## Anopheles Stage 0R — MalariaGEN re-audit and analysis redesign

Reason for reopening: Ag3.10 now explicitly documents a `fontaine-2015-rebuild`
sample set, making a stronger sample-level 2La analysis plausible if sample
metadata, karyotypes, and regional variation can be accessed.

MalariaGEN access: Stage 0R used the versioned Ag3.10 public documentation and
Ag3 API documentation as provenance. Local live API access did not run because
`malariagen-data` could not be installed in the configured Python 3.12
environment: both latest `15.10.0` and repo-recorded `12.0.1` attempts failed
building native dependencies (`llvmlite`, `numba`, and `biotite`) against local
SDK/LLVM tooling. The prior Stage 1A record still documents an additional GCS
anonymous-access blocker for `malariagen-data==12.0.1`.

Sample counts from the Ag3.10 `fontaine-2015-rebuild` public summary:

- arabiensis: 12
- coluzzii: 11
- gambiae: 26
- melas: 4
- merus: 9
- quadriannulatus: 10
- christyi: 0 / not listed

2La-state availability: authoritative sample-level 2La states were not
retrieved. Stage 0R therefore counts all 72 focal Fontaine-associated aggregate
samples as `unknown` for sample-level karyotype purposes and does not substitute
literature-fixed species states for direct calls.

Phased-data availability: Ag3 documentation describes phased haplotype panels
for `gamb_colu_arab`, `gamb_colu`, and `arab`, with haplotypes representing the
two homologues. Fontaine rebuild membership and rare species support were not
verified locally.

Feasibility classification: provisional B/C. A partial modern MalariaGEN
analysis is likely feasible for `gambiae`, `coluzzii`, and `arabiensis`, but a
full six-species Fontaine-compatible sample-level analysis is not yet proven.
The Fontaine-only historical comparison must remain primary until Stage 1R
recovers sample IDs, karyotypes, and regional variation.

Proposed primary design: use only homozygous `2La/2La` and `2L+a/2L+a`
individuals, exclude heterokaryotypes from the primary unphased local-tree
analysis, and compare arrangement versus species support inside 2La against
eligible collinear 2L windows outside the inversion. The working 2La interval
remains `2L:20524058-42165532`; a 1-based 50-kb 2L grid gives 432 fully inside
windows and 434 intersecting windows.

Blockers: exact Fontaine sample IDs not recovered; per-sample 2La karyotypes
not retrieved; Fontaine rebuild haplotype membership not verified; exact
Fontaine 4,063-window coordinates not recovered; local API/package environment
needs repair or Colab execution.

Next step: Stage 1R should run in a working MalariaGEN environment, retrieve
`fontaine-2015-rebuild` sample metadata and `karyotype("2La")`, confirm phased
2L haplotype membership, perform tiny regional SNP/haplotype smoke tests, and
freeze homozygous sample eligibility before any local tree inference. Do not
mark Anopheles complete.
- Species-fixed states: 0.

Strict quartet freeze:

- Candidate strict 2:2 quartets: 0.
- Four-distinct-species strict quartets: 0.
- Frozen file: `data/anopheles_2la/processed/frozen_strict_quartets_stage1a.tsv`.
- SHA256: `b19c6ac1987b0f141b5f8c89d51a12cb8f47b125015d105b39d0054cb8850cbd`.

Strongest designs:

- Design A: none frozen from authoritative sample-level metadata.
- Design B: none frozen; requires gambiae/coluzzii A0 and A1 sample-level calls.
- Design C: none frozen; requires sample-level geography plus 2La calls.

Remaining blocker:

- Resolve authenticated/public access to the official MalariaGEN API config and
  metadata objects, or use an environment with approved MalariaGEN access.

Portability update:

- Stage 1A now includes a metadata-only authentication check at
  `empirical/anopheles_2la/scripts/01_check_malariagen_access.py`.
- Local authenticated reruns should first use
  `gcloud auth application-default login`, then run the check script, then run
  the Stage-1A fetch/freeze pipeline only after sample-set metadata are visible.
- Timestamped Stage-1A attempt summaries are retained under
  `empirical/anopheles_2la/results/provenance/`; failed attempts do not replace
  previously frozen canonical outputs.
- `malariagen-data==12.0.1` remains the package version to test first after
  authentication succeeds.

## Anopheles 2La Stage 1B / Stage 2 Synthetic Freeze

Stage 1A remains blocked on authenticated MalariaGEN metadata access, and the
real frozen strict-quartet table remains empty. No real Anopheles local trees,
SNP calls, haplotypes, sequence-derived topologies, q1/q2/q3, QQS/BQS, NJ
trees, inside/outside genealogy patterns, or published per-window genealogies
were accessed for Stage 1B.

The downstream Stage 1B/Stage 2 machinery was implemented prospectively using
synthetic fixtures only:

- Prospective plan:
  `empirical/anopheles_2la/config/stage1b_analysis_plan.yaml`
- Plan report:
  `empirical/anopheles_2la/results/stage1b_analysis_plan.md`
- Interpretation framework:
  `empirical/anopheles_2la/results/competing_model_predictions.md`
- Success criteria:
  `empirical/anopheles_2la/results/stage2_success_criteria.md`
- Synthetic fixtures:
  `empirical/anopheles_2la/tests/fixtures/`
- Local-tree ingestion:
  `empirical/anopheles_2la/scripts/10_ingest_local_trees.py`
- Tree-inference wrapper:
  `empirical/anopheles_2la/scripts/11_infer_local_trees.py`
- Synthetic validation runner:
  `empirical/anopheles_2la/scripts/run_stage1b_synthetic_validation.py`
- Power analysis:
  `empirical/anopheles_2la/scripts/20_power_analysis.py`

Frozen prospective choices:

- 2La interval: `2L:20,524,058-42,165,532`, AgamP4, 1-based inclusive.
- Primary regions: left flank `15,524,058-20,524,057`, inversion
  `20,524,058-42,165,532`, right flank `42,165,533-47,165,532`.
- Primary flank width: 5 Mb.
- Sensitivity flank widths: 2 Mb and 10 Mb.
- Primary windows: 100 kb non-overlapping physical windows.
- Sensitivity windows: 50 kb and 200 kb.
- Primary statistic: `Delta_arr = f_in - f_out`, where `f_in` and `f_out` are
  frequencies of the predeclared arrangement-predicted topology inside 2La and
  in fixed flanks.
- Primary test: block-aware circular-shift permutation, 10,000 permutations for
  final analysis.
- Diagnostic test: naive per-window permutation, explicitly not primary because
  it can be anti-conservative under linkage.
- Polytomy policy: unresolved quartets are recorded and excluded from topology
  denominators; no random resolution.
- Missing-tip policy: missing/duplicate/malformed cases are retained as
  unusable records with exclusion reasons.
- Block normalization: consecutive windows with the same quartet topology state
  for a frozen comparison collapse into one run of weight one.
- Breakpoint analysis: fixed breakpoints at 20,524,058 and 42,165,532; do not
  move breakpoints based on observed transitions.
- Primary tree estimator wrapper: IQ-TREE 2 pinned as
  `bioconda::iqtree=2.4.0`, executable `iqtree2`, deterministic seed 1729,
  command `iqtree2 -s WINDOW.fasta -seed 1729 -nt AUTO -pre OUTPUT_PREFIX -m
  MFP`; no bootstrap/support threshold is used for the primary classification.
  Structural predictions are not available to tree inference.

Synthetic validation results:

- Null/MSC-like scenario: not positive, `Delta_arr = -0.014`, block-aware
  `p = 0.760`.
- MSRC-positive scenario: positive, `Delta_arr = 0.422`, block-aware `p = 0.005`.
- Broad alternative scenario: not positive for 2La-specific enrichment,
  `Delta_arr = 0.014`, block-aware `p = 0.240`, flagged as broad-alternative-like.
- Design-B positive contrast: changed in the structurally predicted direction.
- Design-B null contrast: `null_like`.

Freeze manifest:

- `empirical/anopheles_2la/results/stage1b_synthetic/stage1b_analysis_freeze.json`
- Previous SHA256:
  `a4c366a9921ac11d3794bf29d7dbb102a278dda543a9650cd4ede1c985aa8fe2`
- Current SHA256: see
  `64c6b70f9046ea5b3b7f493150dfd68cd6ae454beb513c3147599e4d204a4274`

Real Stage 1B remains gated. Stage 1B eligibility is design-agnostic. A
scientifically complete Stage 1A may enable Design A, Design B, Design C, any
combination, or no analysis if no prospectively eligible designs exist. Zero
Design-A quartets is not itself a reason to alter the analysis plan. Real-mode
scripts refuse to run until
`data/anopheles_2la/processed/stage1a_freeze_complete.json` exists,
Stage-1A freeze/provenance validation succeeds, and at least one prospectively
frozen eligible design is nonempty. If A/B/C are all empty, Stage 1B stops with
`NO_ELIGIBLE_PROSPECTIVE_DESIGNS` before topology inspection; this is
inconclusive / underpowered, not evidence against MSRC.

The prior checksum is superseded only for pre-data hardening: design-agnostic
execution gate, tree-inference environment pinning, and test-only network
isolation. It is not a response to observed biological results.

## Songbirds/PIPs Stage 0

Started a new empirical feasibility/provenance audit under
`empirical/songbirds_pips/` with data products under `data/songbirds_pips/`.

Project goal: test, in later stages only, whether independently identified
putative inversion polymorphisms define local genealogy regimes that differ
from genomic background, using Pegan and Winger (2025), *Large Inversion
Polymorphisms are Widespread in North American Songbirds*,
DOI `10.1093/gbe/evaf205`.

Data availability: the OUP supplementary workbook and code are frozen locally.
Dataset S1 contains 174 PIPs across 28 PIP-containing species; 139 are
LD-supported in addition to PCA support and 35 have ambiguous LD. Dataset S4
contains per-individual PIP genotype assignments for 1,411 unique samples in
the 28 PIP-containing species. Article Table 1 reports 35 species and 1,660
analyzed individuals overall; the seven no-PIP species are absent from Dataset
S4. Public sequence reads are available through NCBI SRA BioProjects
`PRJNA1043688` and `PRJNA1130443`, but no public packaged analysis-ready
VCF/BCF, BAM/CRAM, phased haplotypes, sequence alignments, Newick local trees,
or whole-genome trees were identified.

Candidate count: all 174 PIPs were inventoried and scored with
topology-independent criteria. The selected initial set is:

- `Empidonax_flaviventris|CM020533.1`
- `Empidonax_minimus|CM020533.1`
- `Geothlypis_philadelphia|CM019906.1`
- `Zonotrichia_albicollis|CM042600.1`
- `shared_CM027530.1`

Blockers: the cheapest rigorous Stage-1 route depends on obtaining or
regenerating BAM/ANGSD genotype-likelihood/SAF or regional matrix inputs.
Using raw SRA reads alone is feasible but terabyte-scale for the full dataset.
Dataset S4 lacks sex and explicit population labels, and exact shared-PIP
overlap/control coordinates must be reconstructed before cross-species tree
analysis.

Downstream strategy: do not infer trees yet. Stage 1 should first try to obtain
author intermediate BAM/SAF/genotype-likelihood inputs or build a small
selected-sample SRA pilot, then freeze exact focal/control intervals and sample
sets before any local-tree inference.

Status: Stage 0 started and documented; this dataset is not complete.

## Songbirds/PIPs Stage 1A

Stage 1A reconstructed the published shared-PIP analysis design for
`shared_CM027530.1` without read downloads, BAM processing, or tree inference.

Focal object: `CM027530.1`, shared across 11 warbler species mapped to
`GCA_001746935.2` (*Setophaga coronata*). Dataset S5 fixes the author
representative set at 22 samples, with one AA and one BB representative per
species. The representative set joins cleanly to Dataset S4 metadata.

Shared-PIP coordinates were reconstructed from Dataset S1 as the common overlap
among the 11 species-specific PIPs: `CM027530.1:574805-874807`, length 300,002
bp under the source end-minus-start convention. All 11 species-specific PIPs use
the same reference assembly. The exact author shared-COMP interval was not
present in Dataset S1, Dataset S5, or packaged B3 code; the code expects that
interval in a user-provided `shared_pips_meta` file.

Published workflow: select representatives, generate ANGSD SAFs from selected
BAMs for PIP and COMP intervals, compute pairwise 2D-SFS with `winsfs`, convert
to pairwise `dxy`, then plot neighbor-joining genetic-distance trees with
`ape::nj`. This is a distance-based NJ workflow, not a packaged reusable Newick
tree release.

Feasibility decision: the published shared-PIP analysis cannot be reproduced
from public intermediates alone. Public sources provide metadata, coordinates,
code, and raw-read accessions, but no CM027530.1 VCF/BCF, BEAGLE likelihoods,
SAF files, pairwise SFS, dxy matrix, BAM/CRAM, or reusable Newick tree was
found. Stage 1B should first seek author BAM/SAF/dxy/`shared_pips_meta`
intermediates; otherwise a narrow 22-sample regional reprocessing pilot is
required.

Outputs:

- `data/songbirds_pips/metadata/stage1a_shared_pip_manifest.tsv`
- `data/songbirds_pips/processed/stage1a_shared_pip_samples.tsv`
- `data/songbirds_pips/processed/stage1a_shared_pip_interval.tsv`
- `empirical/songbirds_pips/results/stage1a_published_workflow.md`
- `empirical/songbirds_pips/results/stage1a_reusable_intermediates.tsv`
- `empirical/songbirds_pips/results/stage1a_report.md`

Status: Stage 1A complete for feasibility/provenance; no Stage 1B inference has
been started.

## Neoaves chromosome 4 Stage 4B

Stage 4B applies four predeclared quartet weighting treatments to the frozen
chr4 topology tracks: T0 window weighting, T1 topology-run normalization, T2
exclusion of frozen structural-mask windows, and T3 fixed 0.5 weighting inside
the mask versus 1.0 outside. T3 is topology neutral. Primary results use the
primary mask; stringent and inclusive masks are sensitivity analyses.

Block normalization changes dominant topology for N61 (q2 to q1) and N62 (q1
to q2). Primary structural exclusion reduces N61 q2 support but does not change
its q2 dominance. N62 remains q1 dominant after structural exclusion, which
distinguishes its T1 result as a generic block-length effect. Columbea remains
q1 dominant under every treatment. N61 structural effects retain direction
under all three masks and every leave-one-event-out analysis remains q2
dominant, so no single frozen event drives the result.

Full species-tree inference was not run because the repository has focal
quartet summaries but lacks per-locus multi-taxon gene trees and a fixed
empirical estimator configuration. Stage 4B does not identify a true topology
or establish species-tree bias; independent validation is reserved for Stage
4C.

## Neoaves chromosome 4 Stage 4C

The primary independent reference is the genome-wide focal QQS dataset after
prospective exclusion of every chromosome-4-labelled sequence. It contains
54,373–56,247 usable non-chr4 loci per focal clade across 32 chromosome units.
No per-locus Newick trees, standalone frozen species tree, or fixed empirical
species-tree estimator configuration was found.

The non-chr4 reference topologies are q1 for Columbea, q1 for N61, and q2 for
N62. All remain unchanged in every chromosome leave-one-out replicate. N61 T0
and N62 T0 disagree with their respective references; block normalization T1
changes both to the reference topology. Structural exclusion/downweighting
moves both margins modestly toward their references but does not change either
dominant topology. Columbea T0 already agrees with q1, and all treatments remain
q1 dominant while reducing its reference margin.

Interpretation categories are C for Columbea and A for N61 and N62. Focal
empirical species-tree bias from naive dense chr4 weighting, plus block-aware
correction, is supported for N61 and N62 relative to stable non-chr4 quartet
evidence. Partial structural-aware correction is supported for N61; N62 remains
a generic block-length case because structural exclusion does not reproduce the
topology change. This is focal topology-reference validation rather than a full
multi-taxon species-tree re-estimation, and it does not establish a causal
rearrangement mechanism.

## Atlantic cod supergenes Stage 0-5 final freeze

Atlantic cod Stage 0-5 is complete and frozen for manuscript use. The analysis uses the published Matschiner et al. 250-kb local population trees, independently frozen structural arrangement predictions, and an independently frozen collinear baseline population-history tree.

Final primary effects from the exact circular-shift spatial null are: LG01 Delta_D=1.232766, p=0.009174, BH=0.021978; LG02 Delta_D=0.892332, p=0.010989, BH=0.021978; LG07 Delta_D=-0.005981, p=0.555556, BH=0.555556; LG12 Delta_D=1.043170, p=0.038462, BH=0.051282. LG01 and LG02 are classified as robust_positive, LG12 as supportive_positive with the predeclared small-quartet/orientation caveats, and LG07 as weak_or_null.

Stage 5 added a physical-coordinate-aware sensitivity null using the nominal 250-kb grid with missing windows retained as missing. The corresponding p-values are LG01=0.022222, LG02=0.013514, LG07=0.613636, and LG12=0.072727. Manuscript-ready figures, tables, result text, methods text, and the final freeze manifest are under `empirical/atlantic_cod/figures/` and `empirical/atlantic_cod/results/`.

Next manuscript task: integrate the frozen Atlantic cod figure/table and concise Results/Methods text into the manuscript without changing frozen Stage 2-5 analysis definitions.

## Stage 1R — MalariaGEN regional SNP-access and homozygote pilot

Stage 1R uses authenticated MalariaGEN Ag3 release 3.10 data from `fontaine-2015-rebuild` via `gs://vo_agam_release_master_us_central1`. The primary pilot cohort is the frozen set of 53 2La homozygotes; 19 heterokaryotypes remain in metadata but are excluded from the unphased diploid pilot.

The fixed pilot queried eight predefined 50-kb windows on 2L, used direct regional `snp_calls` and `biallelic_snp_calls`, built IBS-style diploid distances, generated NJ trees for visualization, and scored crossed gambiae/coluzzii quartet support. No raw reads were downloaded and no full 2L scan was performed.

Stage 1R decision: **PASS**. Mean pilot C inside/outside: 0.100534 / 0.0039476. Mean pilot D inside/outside: 1 / 0.335214.

## Stage 2R — full MalariaGEN 2L spatial scan

Stage 2R completed the full coordinate-1-anchored 50-kb scan across 2L using authenticated MalariaGEN Ag3 release 3.10 processed SNP calls for the frozen 53 homozygous Fontaine-associated samples. No raw reads were processed.

Windows: 988 total; 432 inside 2La, 554 outside, 2 boundary; 978 usable by the frozen SNP/missingness rule.

Primary crossed signals: `Delta_C=0.1055` with circular p=0.001012; `Delta_M=0.217` with circular p=0.001012. Haplotype validation was concordant for 6/6 deterministic validation windows.

Decision: **SUPERSEDE**. The next step is manuscript integration using the MalariaGEN sample-level analysis as the primary Anopheles result and retaining the Fontaine tree analysis as historical/external comparison.

## Stage 3R — ASTRAL species-tree sensitivity

Stage 3R reused the Stage-2R local 53-individual NJ trees and ASTRAL4's documented multi-individual mapping (`-a`) to test whether linked arrangement-dominated 2La windows alter four-species summary-tree inference. The Stage-2R phased validation text was corrected to 6/6 concordant deterministic validation windows, and the main Anopheles figure was rebuilt with matplotlib axes.

- `T_outside`: (arabiensis,quadriannulatus)|(coluzzii,gambiae)
- `T_inside`: (arabiensis,gambiae)|(coluzzii,quadriannulatus)
- `T_all`: (arabiensis,quadriannulatus)|(coluzzii,gambiae)
- classification: **SUPPORT/BRANCH-LENGTH EFFECT**

The primary Anopheles biological result remains Stage 2R's crossed MalariaGEN result inside 2La; Stage 3R is a sensitivity analysis of propagation into species-level summary-tree inference. The improved Stage 3R figure shows that 2La does not flip the combined species topology, but progressively erodes its quartet and CU support, while 2La-only windows support a different topology with near-unit ASTRAL support.

## Stage 4R — six-species species-tree sensitivity

Stage 4R completed the additive six-species ASTRAL sensitivity analysis using all 72 Fontaine-associated MalariaGEN samples as 144 phased haplotypes. This is now the primary Anopheles species-tree sensitivity analysis. Stage 2R remains the primary clean arrangement-vs-species result, and Stage 3R remains the restricted four-species homozygote sensitivity analysis.

Usable windows: 976 non-boundary windows (546 outside 2La, 430 inside 2La). Stage 4R classification: **SUPPORT/BRANCH-LENGTH EFFECT**.

`T_outside6`: `arabiensis,coluzzii,gambiae,quadriannulatus|melas,merus; arabiensis,coluzzii,gambiae|melas,merus,quadriannulatus; arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae`

`T_inside6`: `arabiensis,coluzzii,gambiae,merus|melas,quadriannulatus; arabiensis,coluzzii,gambiae|melas,merus,quadriannulatus; arabiensis,gambiae|coluzzii,melas,merus,quadriannulatus`

`T_all6`: `arabiensis,coluzzii,gambiae,quadriannulatus|melas,merus; arabiensis,coluzzii,gambiae|melas,merus,quadriannulatus; arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae`

Most inversion-sensitive outside-tree branch: `arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae`.

Stage 4R is finalized after correcting the downweighting endpoint label to 430 inside windows; no further Anopheles empirical analysis is currently required.

## Mai30 large-tree simulation study — Phase 1 audit

Mai et al. (2017) D1 provenance is authenticated through the article, the
authors' `esayyari/simphy-uyen` repository, the author project page, and Dryad
DOI `10.6076/D1RW2G`. The deposit documents 30 ingroup species, 100 replicates
per model condition, 500 gene trees per replicate, and ten D1 conditions. The
published file manifest and SHA-256 checksums are recorded under
`simulations/large_tree/mai30/datasets/provenance/`.

The D1 archive was found locally at `simulations/large_tree/D1-true.tar.gz`
and its SHA-256 matches Dryad exactly. The species-tree, command, log, and
parameter subset has been extracted. All 1,000 condition/replicate records
have recovered SimPhy parameters: generation time 1, generation units, one
global haploid effective population size per replicate, and exact commands.
The three zero-outgroup parameter tables have a shifted-column defect; the
detailed logs recover those values and the audit records that correction.
The SimPhy haploid-to-`msrc-sim` mapping is implemented and tested as
`N_msrc=N_H/2`, with integer finite-history error reported for odd values.
Four no-outgroup original gene-tree collections were selected by observed
discordance and compared with new ordinary MSC controls.

The local `msrc-sim` audit found arbitrary-tree full Newick genealogy output,
branch-specific effective population sizes, selected origin branches, forward
Wright–Fisher/Moran histories, speciation transitions, and suppression
parameters. The MSC controls passed the 500-tree, taxon-set, bifurcation,
finite-branch-length, and complete-coalescence contract. Ordinary MSC root
extension is effectively infinite in the current implementation. Runtime
benchmarks show that original-scale forward Wright–Fisher histories are
currently prohibitive for 500 genes, so the original-scale MSRC study remains
**BLOCKED** pending a validated acceleration or an explicitly labeled
coalescent-timescale-preserving approximation.

The explicitly labeled scaled-WF pilot is implemented with `c=0.001`.
Scaled MSC controls passed for the four selected Mai replicates, with mean RF
differences from the original collections below 0.005. Six of the twelve
primary MSRC collections reached 25 independent histories × 20 loci; the
remaining collections retain validated 20-locus smoke outputs because local
runtime saturated. The resumable execution and ASTRAL analysis are recorded
under `simulations/large_tree/mai30/analysis/results/`. These outputs are a
finite-population scaled approximation and do not remove the original-scale
computational blocker.

The epsilon-mixture experiment now uses the six matched complete MSRC pools
(replicates 007 and 042) and their 500-tree MSC controls at epsilon values
0–1. ASTRAL, quartet-margin, history-bootstrap, branch-length-audit, and
preliminary figure outputs are under `simulations/large_tree/mai30/`. Apparent
quartet endpoint reversals in extreme ILS are currently source-pool finite-
sample candidates with broad history-bootstrap uncertainty; no population
inconsistency claim is supported. Original-scale `wf_diffusion` configurations
are prepared only, pending the simulator's experimental validation gate.

The Mai30 pilot now also reports two direct gene-tree comparisons: discordance
of MSRC gene trees from the fixed true species tree, and distributional
distance between MSRC gene trees and the original Mai gene trees. Complete
500-tree comparisons cover replicates 042 and 007 across all three suppression
treatments; the result tables and reproducible script are under
`simulations/large_tree/mai30/analysis/`.

The Mai30 RF analysis was corrected to use strict unrooted bipartition RF,
independent DendroPy validation, duplicate/mismatched-taxon checks, explicit
unresolved-tree normalization, and condition/replicate matching with leading
zeros preserved. The corrected analysis validates all 300 original no-outgroup
collections but finds only six complete MSRC collections, so the 300-replicate
MSRC result is **BLOCKED**. See
`simulations/large_tree/mai30/analysis/results/MAI30_RF_CORRECTION_REPORT.md`.

## House-mouse t-complex three-treatment ASTRAL update

The full-data t-haplotype-only ASTRAL analysis now uses the frozen 4,046 local
5-kb gene trees and 23 retained tips (15 reconstructed t-haplotype tips plus
eight *M. spretus* outgroup tips). It recovers `Q_SPECIES` with CU = 0.202004,
localPP = 1.000000, and biological quartet frequencies
0.455385/0.267675/0.276940 for q_species/q_t_alt/q_other. This agrees with the
existing matched t-only resampling result. Standard-only remains `Q_SPECIES`
and mixed remains `Q_T_ALT`.

The seven-population t-only treatment is not estimable because SPRE has no
t-haplotype representative. The population figure therefore shows the original
seven-population standard-only and mixed trees alongside a six-population
t-haplotype-only tree with SPRE removed. A separate common-six comparison
allows direct comparison of all three treatments on the same taxon set. The
original standard-only versus mixed seven-population comparison remains
unchanged.
Updated three-treatment figures, input/output provenance, and validation checks are documented under
`empirical/house_mouse_t_complex/`.

## Atlantic cod Stage 10 validity correction and Stage 10b

The original Atlantic cod Stage 10 parameter-to-prediction analysis is
classified as exploratory and invalid for empirical MSRC parameter inference.
It treated SNAPP population-tree internal node ages as realized lineage
coalescence times, even though SNAPP integrates over gene trees. The resulting
`T < d` incompatibility affected 11,673/27,786 audited observations (42.010%).
Its lambda and m estimates should not be used in the paper; the original
`stage10_*` outputs remain frozen for provenance.

Stage 10b replaces that observation model with the validated conditional 2:2
quartet model from `simulations/quartet`, using the original rooted Stage-8b
time tree only to obtain baseline quartet intervals. It fits effective
reduced-model parameters from training quartet outcomes and predicts held-out
contiguous windows. The primary objective is a quartet composite/pseudo-
likelihood, with window-level D(w) summaries and comparators against ordinary
MSC and an LG-specific training-frequency null.

The exact conditional kernel matches the frozen simulation grid to 4.44e-16.
The analysis retains 6,406 eligible quartet outcomes overall and 127 held-out
windows with eligible outcomes; LG12 has no eligible positive baseline
intervals under the fixed scaffold. Pooled held-out MSRC D(w) MAE is 0.07968
versus 0.87533 for ordinary MSC and 0.08032 for the smoothed
training-frequency null. The training-frequency null has better mean
log-predictive probability (-0.296 versus -0.454 for MSRC), so the mechanistic
model does not outperform the non-mechanistic LG-specific predictor on the
primary predictive score.

The effective parameters are weakly identified. The shared lambda profile and
several m/rho profiles hit effective-rate boundaries; interpretation therefore
rests on predictive structure and rho diagnostics rather than precise
biological rates. Stage 10b supports only the limited statement that this
reduced conditional MSRC model captures arrangement-associated quartet
structure and outperforms ordinary MSC in this comparison. It does not provide
an independent second-modality validation, estimate meiotic recombination
rates, or identify the full biological MSRC process. It should remain
supplementary rather than serve as a main-paper quantitative MSRC estimate.

## Exact four-taxon species-tree MSRC theory

The symmetric four-taxon species-tree theorem is complete for equal sister
branches. The scalar reduction factors the exact matrix theorem through the
closed sister-pair coefficients and the existing `H_m` shared-population
kernel. The four conditional classes are 4:0, 3:1, concordant 2:2, and
discordant 2:2; the 16 patterns reduce to these classes under equal branch
symmetry. `SCALAR_FORMULAS.md` records the formulas, the unconditional forward
weight average, the `q2-q3` identity, the conditional topology-reversal
condition, and the MSC limit.

The scalar evaluator agrees with the matrix evaluator below `1e-12` across the
requested rate and duration grid and agrees with `scipy.linalg.expm` for the
pair kernel. Six selected 5,000-locus cells at `p=1/2`, including `1010` at
low and high switching rates, were compared with direct `msrc-sim` output and
showed Monte Carlo agreement. The earlier p-not-equal-to-one-half validation
bundle remains separate evidence for the prior general theory.

The remaining limitation is that the scalar theorem assumes equal sister
branch lengths, zero terminal lengths, symmetric rates, and a fixed terminal
pattern; arbitrary forward pattern weights are handled only after the
conditional theorem.

## Atlantic cod Stage 11 Phlag benchmark

The original Phlag implementation was run from upstream commit
`76597d0e4ebd0e507fdafebf7f5163ef634a36d6` on the frozen Atlantic-cod window
trees, independently for LG01, LG02, LG07, and LG12. It completed for all four
linkage groups and nine eligible focal internal branches per group. One basal
three-taxon branch was excluded because Phlag's own basal-bifurcation collapse
removes that label before focal-edge validation.

Phlag flagged no windows under its default HMM state call. Its union and every
branch-specific result therefore had F1=0 and interval IoU=0 in all four
linkage groups. Using fixed sign-rule comparators, pooled topology-only F1 was
0.7195 with mean LG IoU 0.4270, time-only F1 was 0.7516 with mean LG IoU
0.4322, and the joint rule was F1 0.7939 with mean LG IoU 0.4402. LG07 was
detected by the time-only rule (F1 0.5324, IoU 0.3035) but not by Phlag.

These results are a software benchmark and do not establish detector
superiority. Phlag assumes gene-tree-like inputs and MSC-calibrated quartet
signals, while these inputs are SNAPP MCC population-tree summaries. The
result should remain supplementary until Phlag's calibration on this input
representation is independently validated. Reproduction files and raw Phlag
outputs are under `empirical/atlantic_cod/results/phlag_benchmark/`.

## Atlantic cod Stage 12 boundary detection

Stage 12 adds a deterministic, gap-aware least-squares interval detector with a
predeclared BIC-like penalty, minimum segment size, an explicit no-change
option, and internal/left-censored/right-censored outcomes. Frozen inversion
coordinates are loaded only after primary fitting. Existing Stage 11 sign-rule
and Phlag outputs remain unchanged and are compared as separate baselines.

The primary topology detector selected LG01 9.000001–26.000000 Mb, LG02
18.250001–24.000000 Mb right-censored, LG07 left-censored through 12.000000
Mb, and LG12 left-censored through 13.250000 Mb. The primary time-A detector
selected LG01 9.250001–26.000000 Mb, LG02 18.500001–24.000000 Mb
right-censored, LG07 13.750001–23.000000 Mb, and LG12 1.500001–13.250000
Mb. The primary joint detector selected the same LG07 interval, but a poor
LG12 interval extending beyond the observed domain.

Mean full-reference IoU/F1 were 0.6929/0.7381 for boundary topology,
0.9314/0.9682 for boundary time-A, and 0.7250/0.7559 for boundary joint.
LG07 is recovered accurately by the time-A detector (IoU 0.9566, F1 0.9730),
while topology-only remains a weak and unstable localization signal. The
forced topology fit for LG07 selected a short early outlier, so its apparent
topology signal does not support the long known inversion interval.

The Stage 11 Phlag diagnostic found all 36 saved runs had all-zero paths,
`nan` posterior anomaly probabilities, `State divergence: nan`, and final
focal-edge estimates collapsed to `1e-06`. This is a numerical/calibration
degeneracy on SNAPP MCC population-tree inputs, not evidence of a calibrated
Phlag null state. Stage 12 therefore does not claim that Phlag generally fails
on inversion datasets.

The time result supports the limited statement that SNAPP-derived temporal
contrasts provide useful boundary-localization structure in this Atlantic-cod
benchmark. It does not establish a statistically significant advantage over
all alternatives, and the time statistic is not a realized gene-lineage
coalescence time. Next steps are replication on house mouse, fire ants, and
Anopheles with data types appropriate to each detector's observation model.

## Atlantic cod Stage 12b robust boundary validation

Stage 12b preserves the Stage 12 outputs and reproduces all 20 frozen primary
fit scores and statuses exactly. It adds four prespecified objectives: Model A
is the original Stage 12 joint detector; Model B uses robust median/MAD and
IQR-bounded scaling with equal topology/time block weight; Model C requires
positive inside-minus-outside shifts for both `D(w)` and the arrangement-
informed SNAPP population-tree contrast `A(w)`; and Model D combines separate
one-dimensional topology and time score profiles with fixed equal weights.

The original joint objective fails on LG12 because its highest-scoring
candidate is a right-censored 13.500001–27.250000 Mb interval, giving IoU
0.0314. Candidate profiles show that the objective is driven by a large
topology improvement for the opposite interval and a positive time improvement
for the selected interval; the unconstrained combined scan can therefore
prefer a biologically inconsistent spatial segment. Model B selects
1–13.750000 Mb on LG12 (IoU 0.9147), while the directional Model C selects
1–13.250000 Mb (IoU 0.8967). Model D retains the original wrong LG12 segment.
These alternatives are exploratory and were not tuned to LG12 coordinates.

Across the four linkage groups, mean full-reference IoU/F1 were 0.7250/0.7559
for Model A, 0.9458/0.9780 for Model B, 0.7022/0.7323 for Model C, and
0.7250/0.7559 for Model D. Model B is a useful robustness result but is not a
validated universal detector from four linkage groups. LG07 was recovered at
13.750001–23.000000 Mb by Models A, B, and D with IoU 0.9566 and F1 0.9730;
the direction-aware Model C declined to no detection. Leave-contiguous-block-
out perturbations keep the Model A LG07 interval stable apart from modest
endpoint changes, while robust-model LG07 stability is affected by a small
number of influential blocks.

The arrangement-blind analysis summarizes each published local SNAPP tree by
pairwise-distance, branch-length, and tree-height features, with total-length
normalized sensitivities. The feature selected without inversion coordinates
was mean internal branch length. It recovers LG07 at 13.750001–23.000000 Mb
(IoU 0.9566, F1 0.9730), showing that the LG07 temporal/local-tree-scale
signal is not dependent on the known arrangement labels. These remain
summaries of SNAPP population-tree/local-tree scale, not realized lineage
coalescence times.

Maximum-score detection p-values use 99 deterministic permutations of
contiguous five-window blocks and Holm adjustment within method. This is a
heuristic spatial-dependence diagnostic that accounts for the interval scan;
it is not treated as calibrated boundary significance. The four linkage groups
are too few for strong generalization. Model A remains the reference for the
main analysis; Model B and arrangement-blind branch-length summaries should
remain supplementary pending replication on house mouse, fire ants, and
Anopheles.

Stage 12b outputs, figures, methods, diagnostics, and the reproducibility
manifest are under
`empirical/atlantic_cod/results/boundary_detection_stage12b/` and
`empirical/atlantic_cod/figures/`. The original Stage 11 Phlag diagnostic and
all Stage 12 files remain unchanged.

## Fire-ant chromosome-16 patristic distance track

Stage 9 adds a new descriptive distance track from all 213 published RAxML
local trees. It uses the verified Stage-6C individual-to-group mapping and
calculates patristic distances in substitution units among the four focal
supergene groups. The primary statistic averages within each of the six focal
group pairs and then gives equal weight to the four opposite-haplotype pairs
and the two same-haplotype pairs. It is explicitly not a coalescence-time
estimate and does not assume a molecular clock.

All 213 windows had complete focal-group observations. In the 52 frozen
supergene windows, mean `A_distance` was 0.001718 substitution units, versus
-0.003447 in the 44 chromosome-16 outside windows, for a difference of
0.005165. The corresponding `A_relative` means were 0.2417 and -0.3029, a
difference of 0.5447. The individual-pair-weighted sensitivity gave a
0.004720 `A_distance` difference, and the outgroup-normalized sensitivity
gave a 0.3570 difference.

Using the existing chromosome-16 spatial-null framework, the exact circular
shift p-values were 0.0208 for `A_distance` and 0.0104 for `A_relative`.
Coordinate-membership sensitivity gave p=0.0309 and 0.0206, respectively;
continuous coordinate-length weighting gave p=0.0241 and 0.0603. The
individual-weighted and outgroup-normalized tracks retained the same positive
direction. These are spatial alignment summaries for substitution-unit
distance tracks, not evidence of a molecular-clock divergence-time shift.

Outputs are under
`empirical/fire_ants_chr16/results/stage9_distance/`, with the track figure
under `empirical/fire_ants_chr16/figures/`. The plotted points retain their
observed genomic coordinates and lines are broken across large unobserved
gaps. Existing fire-ant results and figures were not overwritten.

## Anopheles 2La and house-mouse t-complex distance tracks

New dataset-specific spatial analyses were added without changing the frozen
upstream results. Anopheles Stage 12 processes the 976 usable Stage-4R local
haplotype NJ trees and restricts the primary arrangement contrast to within-
species homozygous 2La versus 2L+a haplotypes; heterokaryotypes remain in the
arrangement-blind summaries. Mouse Stage 4 processes all 4,046 published
5-kb IQ-TREE ML trees and restricts the t-versus-standard contrast to mapped
within-subspecies groups. M. spretus is excluded from that primary contrast
and retained in the arrangement-blind all-tip audit.

Both analyses report patristic branch-length distances in the source tree
units. These are substitution-distance summaries, not calibrated coalescence
times and do not assume a molecular clock. The tracks preserve coordinates
and gaps, include normalized contrasts and sample-balanced sensitivities, and
apply frozen inversion intervals only after track construction for spatial
alignment summaries and interval evaluation. The mouse analysis also includes
an explicitly labelled multi-interval scan because its reference contains four
published intervals.

Outputs are under
`empirical/anopheles_2la/results/stage12_distance/` and
`empirical/house_mouse_t_complex/results/stage4_distance/`; figures and
detector outputs are in the corresponding `figures/` directories. These
analyses are descriptive cross-dataset extensions and should remain
supplementary until independent replication establishes how much branch-length
signal is driven by local informative-site density or sampling composition.

The Anopheles arrangement-informed interval scan selected 20.500001--42.200000
Mb for both `A_distance` and `A_relative`, giving IoU 0.9973 against the frozen
2La interval. The arrangement-blind mean-pair-distance scan was broader and
gave IoU 0.4705. The prespecified inside/outside spatial alignment permutation
gave p=0.001 for the arrangement contrast; this is an interval-alignment test,
not a change-point-existence test. The primary contrast was supported by
within-species homozygous groups in gambiae and coluzzii; heterokaryotypes were
not assigned to either arrangement.

The mouse t/standard track showed mean inside-minus-outside `A_distance` values
of 0.00201, 0.00307, 0.00186, and -0.00279 for In(17)1 through In(17)4,
respectively. A single-interval scan selected 5.780--15.865 Mb and therefore
does not represent the four published intervals; its IoU values were 0.116,
0.553, 0.146, and 0.000 for those four intervals. The separately labelled
greedy multi-interval extension is retained for this multi-interval dataset.
The spatial permutation summaries are alignment diagnostics and are not
evidence that a change point exists. These mouse results are limited by the
published 5-kb ML tree design, pseudo-t haplotype construction, and the fact
that substitution branch lengths are sensitive to sequence-scale variation.

Across the two extensions, continuous branch-length summaries add spatial
information beyond discrete topology in Anopheles, while the mouse analysis
shows that a single-interval detector is inadequate for a multi-interval
reference. Neither result is a direct coalescence-time measurement or a
general claim that distance tracks improve boundary recovery in all systems.

## Stage 13 unified spatial boundary detection

Stage 13 adds a shared robustly standardized spatial segmentation module at
`empirical/shared/spatial_boundary_detector.py` and applies it to Atlantic cod,
fire ants, Anopheles, and house mouse. The fixed detector uses an explicit
`2 log(n)` penalty per added free-mean segment and exact dynamic programming
for each fixed segment count. It preserves observed coordinates and gaps,
supports censored endpoints, and reports a separate exact multi-interval model.

Mean union IoU/F1 for topology versus continuous features were 0.627/0.678 vs
0.784/0.846 in cod, 0.530/0.693 vs 0.000/0.000 in fire ants, 0.995/0.997 vs
0.997/0.999 in Anopheles, and 0.143/0.269 vs 0.117/0.200 in house mouse.
The joint detector did not improve over the better single modality. Anopheles
is the clearest concordant case; fire ants favor topology under the common
penalty; mouse requires multi-interval treatment and remains heterogeneous.

Circular-shift results test alignment with frozen reference regions. The
separate change-point permutation diagnostic tests existence of a scanned
change point and is heuristic under strong spatial dependence. Stage 13 does
not treat patristic substitution distances as calibrated coalescence times.
The unified results are suitable for supplementary cross-system analysis;
claims of general superiority require more independent systems and null-model
validation for each spatial design.

## Stage 13 Phlag cross-system benchmark

The original Phlag implementation was audited at version 0.2.0, commit
`76597d0e4ebd0e507fdafebf7f5163ef634a36d6`, with its bundled example retained
as a finite-output check. The new audit is under
`empirical/cross_system_boundary_detection/phlag/` and does not modify the
Atlantic cod Stage 11 raw outputs.

The unchanged cod Phlag runs are classified as numerically degenerate on
SNAPP MCC population-tree inputs: posterior anomaly probabilities are
nonfinite and the state paths are not interpretable. They are not scored as
valid Phlag negatives. Fire ants, Anopheles, and house mouse are classified as
input-incompatible for the primary benchmark because their local trees contain
individual, phased haplotype, or pseudo-t-haplotype tips while the available
reference trees use collapsed group/species/population taxa. No arbitrary tip
collapse or representative selection was introduced.

Therefore Phlag has no valid empirical IoU/F1 score in this cross-system
extension. The existing unified MSRC detector scores remain available in the
comparison table for context, while Phlag invalidity is reported separately.
The result is a compatibility limitation, not evidence that Phlag generally
fails on inversion datasets. A future benchmark needs taxon-compatible
individual-level species trees or an independently validated multi-sample
mapping layer. The Phlag extension is supplementary and does not support a
cross-system superiority claim.

## Stage 13b corrected boundary detection and matched-taxon Phlag pilots

Stage 13b preserves the earlier Stage 12, 12b, and 13 outputs and restores the
validated Atlantic-cod Stage 12b Model B as the primary reference detector.
The four frozen Stage 12b IoUs were reproduced exactly: LG01 0.9676, LG02
0.9444, LG07 0.9566, and LG12 0.9147 (mean 0.9458). Stage 13 changed these
results because it used independent flank means, different scale handling,
and a different penalty; those objectives are not equivalent. Stage 13b also
uses merged interval-set IoU, so disconnected predictions are not evaluated
by their outer hull.

The corrected cross-system outputs are in
`empirical/cross_system_boundary_detection/stage13b/`. The spatial-null table
records the prespecified ordered-block design but does not report p-values
where the full detector refit per null replicate was not completed. This is a
documented computational limitation rather than evidence for significance.

The matched-taxon Phlag pilot selected one fixed representative per taxon for
each replicate and used the unmodified Phlag 0.2.0 executable. The Anopheles
and fire-ant pilots executed with exact matched taxon sets but produced
nonfinite posterior diagnostics, so they are classified as numerically
degenerate rather than negative detections. A house-mouse pilot completed with
finite posteriors and a nonconstant state path; its default state-call pilot
had IoU 0.1326 and F1 0.2342 for one focal edge and one representative
replicate. These pilot results are not a cross-system benchmark: the requested
20-replicate expansion and full edge aggregation remain future work, and the
valid mouse result is explicitly exploratory. No Phlag accuracy is assigned
to the degenerate Anopheles/fire-ant runs.

## Stage 13c final validation

Stage 13c makes the frozen Atlantic-cod Stage 12b Model B executable in the
cross-system validation path. All four predicted intervals reproduce exactly:
LG01 9.00--26.00 Mb, LG02 18.50--24.00 Mb (right-censored), LG07
13.75--23.00 Mb, and LG12 0--13.75 Mb. The previous LG02 discrepancy arose
because the Stage 13 cross-system feature matrix used `D(w)` plus `A(w)`,
whereas Stage 12b Model B uses two topology simplex coordinates plus `A(w)`
with its original robust modality scaling. These are now separate method
identifiers, and the top-20 candidate audit is in
`stage13c/results/lg02_detector_discrepancy.tsv`.

Stage 13c adds BIC-style K=0--4 model-selection records, Huber-loss fire-ant
outlier sensitivity, descriptive leave-block boundary stability, and full
ordered-block maximum-score null pilots. The null pilots use 19 replicates for
the full cross-system grid because larger replicate counts were computationally
infeasible in this execution; their p-values have coarse Monte Carlo
resolution and should not be used as manuscript significance claims.

Three house-mouse Phlag representative sets were run on focal edge N002. All
three had finite posterior probabilities and informative state paths. Their
default state-call scores were IoU 0.1326/F1 0.2342, IoU 0.1256/F1 0.2231,
and IoU 0.1394/F1 0.2446. This remains a one-edge pilot rather than a full
multi-edge, 20-replicate benchmark. Anopheles and fire-ant matched-taxon runs
remain numerically degenerate because posterior probabilities are nonfinite;
their failures are retained as diagnostic results and are not scored as zero
accuracy.

## House mouse Stage 5 unified topology and distance tracks

The new Stage 5 visualization restores the frozen ALL_TIPS quartet scan as the
primary topology signal for all 4,046 chromosome-17 windows. It defines
`D(w) = q_T(w) - q_S(w)`, the negative of the historical `delta_species_alt`
sign; historical source files remain unchanged. The TTT sampling scheme is
retained as a supplementary comparison because it uses a different
status-pattern-specific sampling design. The ALL_TIPS and TTT tracks have
Pearson correlation 0.3286, so they should not be substituted for one another.

The existing Stage 4 `A_distance` and `A_relative` values were copied without
recomputing the underlying ML-tree distances. They are patristic distances in
substitution units, not calibrated divergence times or realized coalescence
times. All topology coordinates and support checks passed; resolved support
averages 0.999707 because unresolved quartets are retained in the source
denominator. The 250-kb matched-resolution figure is recommended for the main
paper because it is readable while preserving observed genomic gaps; the
native 5-kb figure is retained for supplementary local detail. Outputs are
under `empirical/house_mouse_t_complex/results/stage5_unified_tracks/`.

The follow-up multi-resolution analysis adds pooled 25-, 50-, 100-, and
250-kb tracks with common bin boundaries, observed-window counts, base-pair
coverage, and explicit missing-bin handling. Quartet support is pooled from
the frozen ALL_TIPS counts; distance features are means of valid observations.
The 25-kb version retains more local variation but remains noisy, while 100 kb
is smoother at the cost of lower spatial detail and coverage. The 50-kb
version is the initial main-paper candidate; 25 and 100 kb are retained as
readability sensitivities. The multi-resolution comparison is a visualization
analysis only and does not rerun boundary detection. The original 5-kb and
250-kb Stage 5 outputs remain preserved as supplementary outputs.

The arrangement-blind extension adds `B(w) = mean_pair_distance_no_spretus`,
the mean patristic distance among non-*M. spretus* samples in each local ML
tree. It does not use t-versus-standard labels and remains a substitution-unit
branch-length summary rather than a calibrated time or coalescence estimate.
New native 5-kb, 100-kb, and separate 250-kb four-panel figures place B(w) alongside
ALL_TIPS D(w), A_distance, and A_relative. The 250-kb version is clearer for
the main paper; the 5-kb version preserves local spikes and genuine gaps for
supplementary inspection. The frozen Stage 4 values and prior figures were
preserved.

## House mouse Stage 6 distance-discontinuity diagnostic

The Stage 6 source-tree audit found the largest coordinate-blind change in
arrangement-blind B(w) at 14.7375 Mb. The downstream mean is 0.322949 versus
0.009128 upstream, a 35.38-fold increase. Direct recomputation from all 4,046
published Newick trees matches the frozen Stage 4 B(w), A_distance, and
A_relative values exactly within the recorded numerical checks.

The proximate cause is a single terminal branch: downstream of the change
point, every local tree has `Mmd_FRA4_18B.variant` as its longest terminal tip.
Removing that tip as a sensitivity analysis reduces the downstream B(w) mean to
0.007854, close to the upstream value. This verifies that the discontinuity is
driven by a subset of the retained taxa and terminal branch lengths, rather
than by a general increase in all pairwise distances. The exact cause in the
tree-estimation or input-processing pipeline remains unresolved because the
committed data do not include alignment length, SNP count, missing-data
fractions, or per-window IQ-TREE settings. The original distance tracks and
manuscript figures remain unchanged; Stage 6 outputs are under
`empirical/house_mouse_t_complex/results/stage6_distance_discontinuity/`.

## Anopheles Stage 6 three-topology support figure

The new Anopheles figure joins the frozen Stage 2R three-resolution quartet
supports with the frozen Stage 12 NJ distance tracks and the existing Stump et
al. recombination-rate schematic. The supports are `q_species`,
`q_arrangement`, and `q_third` from the same 2,304-informative-quartet design;
the preferred contrast is `D(w) = q_arrangement - q_species`. Valid support
fractions normalize to one and the saved D values agree with the reconstructed
contrast to approximately 1e-12. The Stage 12 distances remain NJ
substitution-distance summaries, not calibrated coalescence times.

Inside the frozen 2La interval, mean arrangement support is 0.9992 versus
0.0006 species support; outside, the corresponding means are 0.1306 and
0.7711. The revised five-panel figure therefore shows a coherent topology
transition aligned with the frozen interval, alongside the independent
recombination reference and the distance contrasts. Existing Anopheles
figures and frozen empirical results were preserved.

## Neoaves chromosome 4 genomic rearrangement and distance tracks

The new Neoaves stage uses the GalGal6 chromosome-4 coordinate axis and the
frozen QQS locus table. It retains three source-defined focal quartet designs
(Columbea, N61, and N62) separately: q1 is the species-history resolution and
q2/q3 are the two alternative resolutions. The displayed contrasts are
therefore D_q2-q1 and D_q3-q1, rather than an invented single cross-clade qH.

The coordinate-aligned local ML trees contain branch lengths, so a fixed mapped
taxon-set arrangement-blind mean patristic-distance track was recalculated for
complete observations. It uses 47 exact Jarvis-to-Stiller mappings and is a
substitution-distance summary, not a calibrated divergence time or coalescence
time. An arrangement-informed A_distance/A_relative statistic was not created
because the repository does not provide an independent binary arrangement
partition for these avian taxa. No validated local recombination-rate track is
available either.

The main output is the reduced, supported topology-plus-distance figure under
`empirical/neoaves_chr4/figures/neoaves_chr4_genomic_rearrangement_tracks.*`;
source structural regions and coordinates are preserved in the stage-6 result
manifest. The analysis passed three deterministic tests. The figure is suitable
for presenting the separate quartet designs and the arrangement-blind distance
sensitivity; a house-mouse-style arrangement-informed A(w) panel requires
additional independent avian arrangement metadata.

## Neoaves Stage 7 arrangement-informed divergence audit

Stage 7 audited the independent chromosome-4 synteny and structural-state
tables against the three frozen QQS focal designs. The structural labels are
region-specific, multistate synteny/anchor classifications and are independent
of the local ML trees. They do not provide a universal binary arrangement
partition with within-state population replication.

The analysis therefore does not claim a primary population-level A(w). It
retains exploratory same-state versus different-state substitution-distance
contrasts only for CR_01, CR_10, and CR_11, where the state table contains a
numerically usable pair design. CR_02–CR_09 are unsupported because their six
mapped structural species are all singleton state classes. The exploratory
tracks use a region-matched B_r denominator and a balanced state-pair
sensitivity analysis; they are not calibrated coalescence times.

The Stage 7 topology figure now shows normalized q1/q2/q3 support separately
for Columbea, N61, and N62, followed by D2=q2−q1 and D3=q3−q1. The global
fixed-taxon B(w) remains a supplementary reference. The primary manuscript
recommendation is to use the separate QQS topology designs and structural
annotation; a defensible arrangement-informed avian A(w) requires additional
independent, replicated arrangement-state sampling.

## House-mouse Stage 7 influential-sample sensitivity

Stage 7 recalculated the chromosome-17 distance and ALL_TIPS quartet tracks
directly from all 4,046 frozen local ML trees. The original Stage 4 values were
reproduced to floating-point precision. The primary analysis retains
`Mmd_FRA4_18B.variant`; a separate leave-one-out analysis removes it
consistently from both the arrangement-informed numerator and the
arrangement-blind denominator.

The Stage 6 discontinuity at 14.7375 Mb remains reproducible. B(w) changes
from about 0.00913 upstream to 0.32295 downstream, while the leave-one-out
means are about 0.00916 and 0.00785. The leave-one-out A_distance and
A_relative tracks were recalculated rather than changing only the denominator.
The ALL_TIPS quartet track was also recomputed separately after exclusion.

No concrete tree-label, branch-format, or processing error was found in the
available files. Alignments, per-window SNP counts, filtering summaries, and
complete IQ-TREE settings are unavailable, so biological divergence and a
region-specific tree-estimation or ascertainment effect remain unresolved.
The original-data and sensitivity figures should both be retained; the
leave-one-out figure should not be presented as a correction.
## Neoaves chromosome 4 Stage 8: phylogenetically controlled arrangement divergence

Stage 8 adds a separate residual-distance analysis under
`empirical/neoaves_chr4/results/stage8_phylogenetic_distance_adjustment/`.
It uses the six exact structural-state species mappings validated in Stage 7,
the independently supplied region-specific multistate synteny labels, and
pair-specific background medians from chromosome-4 local ML-tree patristic
distances outside each focal region. Broad and strict outside-region baselines,
relative residual contrasts, leave-one-species-out recomputation, and
descriptive concordance with the separate Columbea/N61/N62 QQS designs are
reported. The Stage 7 raw A-distance values reproduce exactly; no Stage 7
outputs were modified.

The structural-state sample could not be expanded beyond six exact local-tree
taxa without unsupported state imputation. CR_01, CR_10, and CR_11 remain the
only exploratory regions with same-state and different-state comparisons;
CR_01 has one same-state pair, while CR_10 and CR_11 reuse the same six species
across shared pairs. Accordingly, the adjusted contrasts are not independent
population-level evidence and are not direct estimates of recombination rate,
divergence time, or coalescence time. All distances remain substitution-unit
patristic distances from local ML trees. The Stage 8 result is suitable as a
conservative supplementary sensitivity analysis; the Stage 7 topology overview
remains the primary Neoaves figure until broader independent structural-state
replication is available.
## House mouse Stage 8: extreme ML branch-length validation

Completed `empirical/house_mouse_t_complex/scripts/08_branch_length_validation.py`
and the results under
`empirical/house_mouse_t_complex/results/stage8_branch_length_validation/`.
The approximately 40–53 substitution-unit terminal branches for
`Mmd_FRA4_18B.variant` are present in the frozen 4,046 local ML Newick trees;
2,748 windows exceed 1 substitution unit, 1,678 exceed 5, 309 exceed 10, 212
exceed 20, and 83 exceed 30. The all-sample B(w) mean changes from 0.00912753
upstream of 14.7375 Mb to 0.322949 downstream, while leave-one-out means are
0.00916046 and 0.00785359. Stage 4 B, A_distance, and A_relative values were
reproduced to floating-point precision.

The source archive comparison is informative: in representative windows, the
extreme branch is present in the coverage-plus-allele-ratio-filtered ML tree
product but is ordinary or absent in the PASS-SNP and coverage-only products.
This verifies a filtering-product-specific anomaly, while the precise cause
remains unresolved because per-window alignments and complete IQ-TREE command
provenance are unavailable. The original all-sample B(w) remains a
reproducible audit statistic with a filtering/provenance caveat; leave-one-out,
median, and trimmed tracks are sensitivity analyses rather than corrections.
Five Stage 8 regression tests pass. Original Stage 4–7 outputs and figures
were preserved.

## Neoaves chromosome 4 Stage 9: full-taxon background divergence

Stage 9 adds `empirical/neoaves_chr4/scripts/09_full_taxon_background.py` and
results under `empirical/neoaves_chr4/results/stage9_full_taxon_background/`.
The chromosome-4 trees contain 363 distinct observed taxon labels across
6,140 unique loci; only 232 loci contain all 363 labels. The historical
Stage-6 fixed mapped set contains 47 observed taxa and 2,707 complete loci.
Its distance values reproduce the frozen Stage-6 track to machine precision.

The new analysis reports all-observed mean patristic distance, the complete
47-taxon reference, a fixed 108-taxon set defined only by presence in at least
95% of loci, and pair-standardized residual summaries using pair-specific
strict-outside-region medians. The preferred full-taxon figure track is
`B_observed`, selected by a coverage rule rather than by structural-region
signal. All distances remain ML substitution units and are not calibrated
divergence or coalescence times. Missing taxa are retained as missing
observations, and the six-species Stage-8 arrangement-adjusted analysis remains
separate. The full-taxon audit is suitable as a supplementary/background
analysis; causal rearrangement interpretation remains limited by
taxon-composition variation and the absence of full-taxon arrangement labels.

## Neoaves chromosome 4 Stage 10: multi-resolution genomic tracks

Stage 10 adds chromosome-anchored 100-kb, 250-kb, and 500-kb aggregation of
the frozen Columbea, N61, and N62 quartet-support tracks and the Stage-9
`B_observed` background distance. Outputs are under
`empirical/neoaves_chr4/results/stage10_multiresolution_tracks/`, with new
figures in `empirical/neoaves_chr4/figures/`.

The 250-kb figure is the primary candidate. It has 361 occupied background
bins; the 100-kb version has 868 occupied of 902 grid bins, and the 500-kb
version has 181 occupied bins. Binned background tracks correlate 0.869
(100 vs 250 kb) and 0.897 (250 vs 500 kb). Broad q1/q2/q3 patterns are stable,
although N61 D3 is less stable at the finest comparison. The observed-taxon
background track remains descriptive: its 250-kb mean correlation with mean
taxon count is -0.3770. The topology-only 250-kb figure is retained as a
conservative alternative, and B_observed should be interpreted as an ML
substitution-distance summary rather than a divergence or coalescence time.

## Neoaves chromosome 4 Stage 11: Phlag comparison

Stage 11 adds a local-first Phlag benchmark under
`empirical/neoaves_chr4/results/stage11_phlag_benchmark/`. The upstream Phlag
v0.2.0 implementation was preserved at commit
`76597d0e4ebd0e507fdafebf7f5163ef634a36d6`. The repository's Stage-6 QQS
vectors are not silently treated as Phlag inputs because they are three
source-defined focal summaries rather than one universal labeled Phlag edge
sequence.

An unmodified default Phlag run completed on a matched Columbea pilot made
from 5,619 original local trees with one fixed representative per role
(Columba_livia, Phoenicopterus_ruber, Tauraco_erythrolophus, and
Gallus_gallus). Its posterior values were finite and its state path was
informative; 1,200 windows were assigned state 1 and runtime was approximately
48.94 seconds. The pilot is a valid software/input execution, but its reduced
taxon sampling is not a reproduction of the published all-taxon avian
analysis. No per-locus published Phlag prediction file was available in the
local checkout, so published-result reproduction is explicitly unavailable.

The fixed locus-level detector selected 44.13–56.82 Mb for both the Columbea
topology-only and topology-plus-`B_observed` configurations. Exploratory
post-fit metrics against the union of 11 canonical regions were IoU 0.565,
F1 0.722 for the detector and IoU 0.358, F1 0.528 for the matched Phlag pilot.
These are not independent accuracy claims because the canonical regions are
an imperfect biological reference and the pilot and detector use different
taxon summaries. The stage therefore supports a reproducible method-agreement
and feasibility benchmark, not a superiority claim for either method.

## Neoaves chromosome 4 Stage 11B: repaired Phlag benchmark

Stage 11B repairs the earlier Phlag coordinate failure without deleting its
provenance. The committed Stage 11 manifest and gene-tree input were empty,
while the retained raw output had 6,096 state/posterior positions; its zero
interval was therefore invalid and must not be interpreted as a Phlag negative.
Those files are archived under
`empirical/neoaves_chr4/results/stage11_phlag_benchmark/provenance_invalid_stage11/`.

The repaired input and manifest are written together in genomic order. They
contain 5,619 identical positions and 5,619 exact four-taxon local trees; 477
of 6,096 source Columbea loci were excluded because a selected representative
was missing. Official Phlag v0.2.0 (commit
`76597d0e4ebd0e507fdafebf7f5163ef634a36d6`) returned 5,619 finite posteriors
and 1,200 anomalous state-1 windows. With a prespecified 250-kb physical-gap
merge rule these form seven intervals: 25.02–32.30, 32.67–32.68,
33.46–34.48, 44.13–46.42, 46.77–47.54, 47.80–48.58, and 49.05–56.82 Mb.

The matched detector is explicitly a categorical four-taxon adaptation, not
the original multi-taxon QQS detector. It selected 25.02–56.82 Mb; the
unmatched original multi-taxon exploratory detector remains 44.13–56.82 Mb.
The corrected outputs report both method classes and exact interval-merging
sensitivity. The local Phlag run took 82.5 seconds and used approximately
602 MB peak child RSS on this machine. Because the matched detector and Phlag
use different statistical objectives, and the structural annotations are
overlapping synteny-derived regions rather than independent breakpoint truth,
these results do not support a superiority claim.
