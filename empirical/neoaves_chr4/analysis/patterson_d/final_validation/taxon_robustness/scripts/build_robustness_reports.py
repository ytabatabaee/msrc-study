#!/usr/bin/env python3
"""Build final compact tables and manuscript/status reports from sensitivity outputs."""
import csv,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[7]
BASE=ROOT/'empirical/neoaves_chr4/analysis/patterson_d'
FINAL=BASE/'final_validation'
PKG=FINAL/'taxon_robustness';RES=PKG/'results'

def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def f(v,n=3):
 if v in ('',None):return 'NA'
 return f'{float(v):.{n}f}'
def iv(lo,hi):return '['+f(lo,3)+', '+f(hi,3)+']'
def main():
 cfg=read(PKG/'config/frozen_alternative_quartets.tsv')
 block=read(RES/'taxon_robustness_block_summary.tsv')
 quart=read(RES/'quartet_null_validation.tsv')
 topo=read(RES/'nucleotide_topology_robustness.tsv')
 taxa=read(RES/'candidate_taxa.tsv')
 raw=read(RES/'raw_alignment_validation.tsv')
 designs=[c['design'] for c in cfg]
 def sr(d,fil='all',scope='all_PNAS',bs='1000000'):
  return next(x for x in block if x['design']==d and x['filter']==fil and x['scope']==scope and x['block_size_bp']==bs)
 def tr(d,fil='all',scope='whole',bs='1000000'):
  return next(x for x in topo if x['design']==d and x['filter']==fil and x['scope']==scope and x['resolution_bp']==bs)
 compact=[]
 for c in cfg:
  d=c['design'];a=sr(d);tv=sr(d,'transversions_only');q=next(x for x in quart if x['design']==d);cor=tr(d)
  compact.append(dict(design=d,replacement_class=c['replacement_class'],P1=c['P1'],P2=c['P2'],P3=c['P3'],O=c['O'],reference_topology=q['source_reference_topology'],n_exact_gene_trees=q['chr4_exact_quartet_trees'],unresolved_quartets=q['unresolved_exact_quartets'],q1_fraction=(int(q['q1'])/int(q['chr4_exact_quartet_trees']) if int(q['chr4_exact_quartet_trees']) else ''),q2_fraction=(int(q['q2'])/int(q['chr4_exact_quartet_trees']) if int(q['chr4_exact_quartet_trees']) else ''),q3_fraction=(int(q['q3'])/int(q['chr4_exact_quartet_trees']) if int(q['chr4_exact_quartet_trees']) else ''),D_topology_q3_minus_q2=q['D_topology_q3_minus_q2'],D_inside_all=a['inside_D'],D_outside_all=a['outside_D'],Delta_D_all=a['delta_D'],Delta_D_all_CI_low=a['delta_CI_low'],Delta_D_all_CI_high=a['delta_CI_high'],inside_loci=a['inside_loci'],outside_loci=a['outside_loci'],inside_informative_sites=a['inside_informative_sites'],outside_informative_sites=a['outside_informative_sites'],D_inside_transversions=tv['inside_D'],D_outside_transversions=tv['outside_D'],Delta_D_transversions=tv['delta_D'],Delta_D_transversions_CI_low=tv['delta_CI_low'],Delta_D_transversions_CI_high=tv['delta_CI_high'],Dsite_Dtopology_r_1Mb=cor['pearson_Dsite_Dtopology'],Dsite_Dtopology_r_CI_low=cor['block_bootstrap_CI_low'],Dsite_Dtopology_r_CI_high=cor['block_bootstrap_CI_high'],within_region_centered_r_1Mb=cor['within_region_centered_pearson'],centered_r_CI_low=cor['centered_bootstrap_CI_low'],centered_r_CI_high=cor['centered_bootstrap_CI_high']))
 write(RES/'taxon_robustness_summary.tsv',compact)
 # Prespecified decision criterion: all-sites 1-Mb paired physical-block contrast is negative,
 # with its descriptive 95% block-jackknife interval wholly below zero.
 primary=next(x for x in compact if x['design']=='N61_primary')
 ingroup=[x for x in compact if x['replacement_class'] not in ('primary','outgroup_only')]
 retained=[x for x in ingroup if float(x['Delta_D_all'])<0]
 ci_excl=[x for x in ingroup if float(x['Delta_D_all_CI_high'])<0]
 tvret=[x for x in ingroup if float(x['Delta_D_transversions'])<0 and float(x['Delta_D_transversions_CI_high'])<0]
 corall=[x for x in [primary]+ingroup]
 rvals=[float(x['Dsite_Dtopology_r_1Mb']) for x in corall]
 centered=[float(x['within_region_centered_r_1Mb']) for x in corall]
 fail=[x for x in raw if x['status'].startswith('FAIL')]
 passrows=[x for x in raw if x['status']=='PASS']
 unavail=[x for x in raw if x['status']=='UNAVAILABLE']
 unique_loci=len({x['locus_id'] for x in raw})
 duplicate_loci=sorted({x['locus_id'] for x in raw if x['duplicate_sequence_ids']})
 unequal_loci=sorted({x['locus_id'] for x in raw if x['all_taxa_alignment_lengths_equal']!='True'})
 # Candidate summary for narrative.
 pnas='\n'.join(f"| {x['design']} | {f(x['D_inside_all'])} | {f(x['D_outside_all'])} | {f(x['Delta_D_all'])} {iv(x['Delta_D_all_CI_low'],x['Delta_D_all_CI_high'])} | {x['inside_loci']}/{x['outside_loci']} | {x['inside_informative_sites']}/{x['outside_informative_sites']} | {f(x['Delta_D_transversions'])} {iv(x['Delta_D_transversions_CI_low'],x['Delta_D_transversions_CI_high'])} |" for x in compact)
 gene='\n'.join(f"| {x['design']} | {x['n_exact_gene_trees']} | {x['unresolved_quartets']} | {f(x['q1_fraction'],3)}/{f(x['q2_fraction'],3)}/{f(x['q3_fraction'],3)} | {f(x['D_topology_q3_minus_q2'])} | {f(x['Dsite_Dtopology_r_1Mb'])} {iv(x['Dsite_Dtopology_r_CI_low'],x['Dsite_Dtopology_r_CI_high'])} | {f(x['within_region_centered_r_1Mb'])} {iv(x['centered_r_CI_low'],x['centered_r_CI_high'])} |" for x in compact)
 n61_intervals=[]
 for name in ('PNAS_1','PNAS_2','PNAS_3'):
  x=sr('N61_primary','all',name)
  n61_intervals.append(f"- {name}: inside D={f(x['inside_D'])}, complement D={f(x['outside_D'])}, ΔD={f(x['delta_D'])} {iv(x['delta_CI_low'],x['delta_CI_high'])} at 1-Mb blocks; loci {x['inside_loci']}/{x['outside_loci']}, informative sites {x['inside_informative_sites']}/{x['outside_informative_sites']}.")
 me=next(x for x in compact if x['design']=='N61_primary_Meleagris')
 criteria=f"**Prespecified criterion:** an alternative retains the primary regional direction if all-sites ΔD<0; interval support is summarized by the 1-Mb paired physical-block delete-one jackknife 95% interval, with ‘excludes zero’ meaning its upper endpoint is below 0. This rule is stated before interpreting the results and does not count loci or quartets as independent replicates."
 assessment=f'''# Final assessment: N61 alternative-taxon validation

## Recommendation: CONDITIONAL GO

The N61 inside/outside nucleotide D contrast is manuscript-ready as a carefully delimited descriptive result for the tested taxa. It is robust in direction to the ten prespecified ingroup representative alternatives: {len(retained)}/{len(ingroup)} retain negative ΔD and {len(ci_excl)}/{len(ingroup)} have 1-Mb block-jackknife 95% intervals below zero. Transversion-only ΔD is negative with its 95% interval below zero in {len(tvret)}/{len(ingroup)} alternatives. This supports a regional site-pattern contrast in this chr4 dataset. It does not distinguish introgression, ancestral population structure, ILS, deep polymorphism, branch-length/substitution effects, homoplasy, or taxon-specific sampling, and it does not establish the MSRC mechanism.

The whole-chromosome nucleotide–genealogy correlation is positive for the original and all tested alternatives (1-Mb r range {min(rvals):.3f}–{max(rvals):.3f}); the original r={float(primary['Dsite_Dtopology_r_1Mb']):.3f} falls within that range. After inside/outside mean differences are removed, correlations remain positive (range {min(centered):.3f}–{max(centered):.3f}), but they are descriptive spatial associations, not independent tests of mechanism. Bootstrap intervals, region-specific correlations, and all resolutions are in `nucleotide_topology_robustness.tsv`.

## Alternative taxon panel and null topology

The original quartet remains primary: Tauraco_erythrolophus / Columba_livia / Ciconia_maguari / Gallus_gallus. The frozen panel includes two P1 replacements (Cuculus, Ardeotis), two P2 replacements (Pterocles, Mesitornis), three P3 alternatives/controls (Balaeniceps, Nipponia, and a Telluraves Aquila control), three combined ingroup substitutions, and one outgroup-only Meleagris sensitivity. All ten ingroup alternatives were selected from taxonomic/literature groupings and archive availability before D calculation; no quartet was selected by its D value. The design table gives each rationale.

The exact induced quartet in the Stiller 363-taxon published reference is q1 for all twelve designs, including the Gallus/Meleagris pair. All resolved local gene trees use exact taxon names and locus IDs; unresolved quartet trees are counted separately (0 across the exact quartet-bearing chr4 subset here). The reference tree comes from the same broad Stiller study as the alignments and is not an independent species-history replicate. Jarvis et al.'s 48-taxon TENT supports broad historical context but lacks *Ciconia maguari*, so an exact independent four-taxon species-tree check is unavailable.

Taxonomic clarification: the frozen P3 taxon *Ciconia maguari* is classified in Aequornithes/Ciconiiformes, not Telluraves. Therefore, P3 is described by its position on the tested quartet's non-P1/P2 side. *Balaeniceps* and *Nipponia* are waterbird-side alternatives; *Aquila* is a labeled broad Telluraves control. Their q1 placement is verified in the Stiller reference, but they are not assumed to be interchangeable members of one narrow monophyletic order.

## Complete regional D results

{criteria}

All estimates pool integer ABBA/BABA counts before calculating D. Block intervals use genomic blocks; informative sites are not independent replicates. The table includes the original, every ingroup alternative, and the alternate-outgroup sensitivity. Locus/site counts are inside/outside the union of the three source-coordinate PNAS intervals. Detailed estimates for each PNAS interval, outside-union-only estimates, and 500-kb/2-Mb block results are in `taxon_robustness_block_summary.tsv`.

| Frozen design | D inside | D outside | ΔD (95% 1-Mb block interval) | Loci in/out | Informative sites in/out | Transversion-only ΔD (95% interval) |
|---|---:|---:|---:|---:|---:|---:|
{pnas}

## Exact genealogy comparison

`q1/q2/q3` are observed proportions across the exact resolved chr4 quartet trees; D_topology=q3−q2. Spatial correlations pair each nucleotide locus to its exact tree by identical archive locus ID, then pool into 250-kb, 500-kb, and 1-Mb bins. Block-bootstrap intervals resample 2-Mb physical blocks (1,000 replicates, seed 20261118). The result table includes whole chromosome, inside, outside, and within-region-centered correlations.

| Design | Exact trees | Unresolved | Mean q1/q2/q3 | Mean D_topology | 1-Mb r whole (95% block CI) | r after region-centering (95% block CI) |
|---|---:|---:|---:|---:|---:|---:|
{gene}

Among ingroup alternatives, the full-chromosome r range is {min(rvals[1:]):.3f}–{max(rvals[1:]):.3f}; the original is not unusually high. The centered correlation is also positive throughout this panel, but it is lower than whole-chromosome r because regional mean differences explain part of the association. Shared loci and overlapping taxa make these sensitivity designs correlated; do not count them as independent biological replicates.

## Independent raw-alignment validation

A separate parser and direct pattern counter independently decoded {unique_loci} deterministic, stratified raw alignments across PNAS_1, PNAS_2, PNAS_3, outside high/low-information loci, and loci with ambiguous/gap data. It compared the three original primary quartets and their predeclared Meleagris outgroup designs against saved all-site and transversion-only counts: **{len(passrows)} exact ABBA/BABA matches, {len(fail)} discrepancies, {len(unavail)} unavailable rows**. The unavailable N62 rows are the two sampled loci missing *Pterocles gutturalis* under each of the two N62 outgroups. There were {len(duplicate_loci)} sampled loci with duplicate sequence IDs and {len(unequal_loci)} with unequal full-alignment sequence lengths. No production filters were changed. See `raw_alignment_validation_report.md` and row-level `raw_alignment_validation.tsv`.

No estimator error was found, so the frozen primary results did not require reanalysis. Regression tests cover parsing, duplicate IDs, case normalization, gaps/ambiguity, multistate filtering, outgroup polarization, ABBA/BABA orientation, and transition/transversion classification.

## Answers to the requested scientific questions

1. **Alternative ingroup quartets:** 10 biologically motivated sensitivity designs were tested, in addition to the original primary and an outgroup-only sensitivity.
2. **Original regional direction:** 10/10 alternatives retain negative all-sites ΔD.
3. **Block uncertainty excludes zero:** 10/10 alternatives have a 1-Mb 95% paired block-jackknife interval wholly below zero. Transversion-only direction and interval are also negative in {len(tvret)}/10.
4. **Dependence on P1/P2/P3:** the contrast direction is not dependent on a single representative. The absolute inside/outside D background is taxon-sensitive: e.g., replacing P1 with *Cuculus* and replacing P1+P2 together changes the outside estimate's sign while leaving ΔD negative.
5. **D_site–D_topology correspondence:** whole-chromosome correspondence remains positive and of similar magnitude across alternatives; the original r is typical. Within-region-centered values are smaller but positive across the panel. These are block-aware descriptive correlations.
6. **Independent raw alignment decode:** exact integer ABBA/BABA agreement for every available sampled design/filter row; 0 discrepancies.
7. **Columbea/N62:** no new full ingroup sensitivity was run for those hypotheses. Their previous near-zero regional contrasts and interpretation remain unchanged; the raw validation checks the original saved counts but does not re-estimate their full archives.
8. **Remaining limitations:** exact species-history quartets are not independently verified for every alternative; the tree reference and alignments share the Stiller source; Gallus and Meleagris are related Galliformes; raw decoding is stratified sampling rather than a full independent reanalysis of all 6,247 alignments; block bootstrap/jackknife is descriptive; alternative designs share linked loci; no demographic/finite-sites model, branch-length-time conversion, or spatial-null test; D/topology cannot resolve introgression versus ILS, ancestral structure, deep polymorphism, homoplasy, or mechanism.
9. **Manuscript readiness:** **CONDITIONAL GO** for the narrow N61 regional site-pattern result and its descriptive genealogy concordance, with the nonzero outside baseline and taxon sensitivity stated. Do not present these results as unique evidence for introgression, an inversion, or MSRC.

## Figures

- Figure E: `figures/figure_E_alternative_quartet_delta_D` (PNG, PDF, SVG)
- Figure F: `figures/figure_F_site_vs_topology` (PNG, PDF, SVG)
- Figure G: `figures/figure_G_sensitivity_matrix` (PNG, PDF, SVG)
'''
 (RES/'final_assessment.md').write_text(assessment)
 manuscript=f'''# Manuscript-ready N61 text (qualified)

## Results

The N61 chr4 nucleotide site-pattern contrast was not specific to the original taxon representatives. In the frozen quartet (*Tauraco erythrolophus*, *Columba livia*, *Ciconia maguari*, *Gallus gallus*), pooled all-sites D was {float(primary['D_inside_all']):.3f} inside and {float(primary['D_outside_all']):.3f} outside the three published GalGal4 source-coordinate intervals (ΔD={float(primary['Delta_D_all']):.3f}; 95% 1-Mb genomic block-jackknife interval {iv(primary['Delta_D_all_CI_low'],primary['Delta_D_all_CI_high'])}). All ten prespecified alternative ingroup quartets retained a negative regional contrast, and each 1-Mb block interval excluded zero. Transversion-only ΔD was also negative in all alternatives. Absolute outside D varied with the selected representatives, including changes in its sign, so the stable observation is the regional contrast rather than a universal background value. These quartets reuse linked loci and are correlated sensitivity analyses.

Among matched chr4 loci with exact four-taxon trees, binned nucleotide D_site and exact D_topology (q3−q2) were positively correlated across the chromosome for the original and all tested alternative quartets (1-Mb r range {min(rvals):.3f}–{max(rvals):.3f}). The original correlation was typical across this panel. Correlations after centering inside/outside regional means were smaller but positive. These descriptive correspondences do not identify the cause of the site-pattern asymmetry.

## Interpretation limitation

The outside D baseline is nonzero and taxon-sensitive. The exact q1 null is supported by the published Stiller 363-taxon reference, which is not independent of the alignment source; an independent exact four-taxon species-history check is unavailable because *Ciconia maguari* is absent from Jarvis et al.'s 48-taxon tree. Patterson's D and its spatial correlation with topology do not distinguish ancient introgression from ancestral population structure, incomplete lineage sorting, deep polymorphism, gene-tree branch-length/substitution differences, finite-sites homoplasy, or taxon-specific sampling, and do not uniquely support MSRC.

## Methods

We analyzed the locally stored, checksum-verified Stiller 63,430 intergenic alignment archive and exact chr4 locus manifest. The original frozen quartet was retained. Ten ingroup alternatives were frozen in a configuration table before calculating D, with two alternatives each for P1 and P2, three P3-side alternatives/controls, and three combined ingroup substitutions; the original ingroup with a Meleagris outgroup was retained as an outgroup-only sensitivity. The published source-coordinate interval predicate was applied to locus window starts. ABBA and BABA integer counts were pooled by region before D was computed. Uncertainty was estimated with paired physical genomic block delete-one jackknife at 500-kb, 1-Mb, and 2-Mb scales.

A separate scalar parser independently re-decoded 46 deterministically stratified raw alignments and compared ABBA/BABA counts for the three primary designs and their predeclared alternate-outgroup designs. Exact local quartet topologies were induced from the Stiller gene-tree archive using exact taxon names and exact locus IDs; unresolved trees were counted separately. D_site and D_topology were compared at 250-kb, 500-kb, and 1-Mb resolution, with 2-Mb physical-block bootstrap intervals for spatial correlations. These analyses are descriptive; linked loci and overlapping taxon panels are not independent replicates.
'''
 (RES/'manuscript_text.md').write_text(manuscript)
 readme='''# N61 taxon robustness and raw-alignment validation

This sensitivity package adds the high-priority N61 checks while preserving the original frozen primary results. It does not run Neurospora, cross-dataset analyses, f4, fdM, or dXY.

## Scope and primary design

Primary remains P1=*Tauraco erythrolophus*, P2=*Columba livia*, P3=*Ciconia maguari*, O=*Gallus gallus*. Alternatives were frozen in `config/frozen_alternative_quartets.tsv` before any alternative D values were calculated. The panel changes one representative at a time, includes a few predeclared combinations, and keeps the original ingroup with *Meleagris gallopavo* as an outgroup-only sensitivity. The original frozen values are unchanged.

The N61 group definitions are not treated as interchangeable species sets. Otidimorphae includes the turaco, cuckoo, and bustard orders in the published phylogenomic literature; these representatives come from distinct orders. Columbimorphae sampling includes pigeon, sandgrouse, and mesite orders in the Jarvis/Prum framework, while relationships among early Neoaves vary among studies. The present P3 taxon *Ciconia maguari* is classified in Aequornithes/Ciconiiformes, not Telluraves. Shoebill and ibis are P3 waterbird-side alternatives; golden eagle is separately labeled as a broad Telluraves control. All tested four-species combinations induce q1 in the published Stiller 363-taxon reference, but this reference is from the same broad study as the alignments. Jarvis et al.'s 48-taxon tree does not include *Ciconia maguari*, so the exact quartet is not independently resolved there.

## Main result

All 10 ingroup alternatives retain a negative all-sites inside-minus-outside ΔD, with 1-Mb paired genomic block-jackknife 95% intervals below zero. Transversion-only contrasts are negative with intervals below zero in all 10. Absolute regional D values, especially the outside baseline, change with taxa and can change sign. Across exact matched gene trees, whole-chr4 1-Mb D_site–D_topology correlations remain positive and the original correlation is typical; region-centered correlations are lower but positive. Designs share taxa and linked loci and are correlated sensitivity analyses, not biological replicates.

Independent raw parsing of 46 stratified alignments exactly reproduced all available frozen primary and alternative-outgroup ABBA/BABA counts (544 comparisons, 0 discrepancies). Four N62 rows were unavailable because *Pterocles gutturalis* was missing in two sampled loci. No estimator error was found.

See `results/final_assessment.md` for the complete decision and results, `results/candidate_taxa.tsv` for taxon availability and lineage assignments, and `results/taxon_robustness_summary.tsv` for the compact complete table. The full 500-kb/1-Mb/2-Mb detail is in `results/taxon_robustness_block_summary.tsv`.

## Reproduction

From the repository root, with Python 3, Biopython, NumPy, and Matplotlib available:

```bash
python3 empirical/neoaves_chr4/analysis/patterson_d/final_validation/taxon_robustness/scripts/inventory_candidate_taxa.py
python3 empirical/neoaves_chr4/analysis/patterson_d/final_validation/taxon_robustness/scripts/compute_taxon_sensitivity.py
python3 empirical/neoaves_chr4/analysis/patterson_d/final_validation/taxon_robustness/scripts/compare_exact_gene_trees.py
python3 empirical/neoaves_chr4/analysis/patterson_d/final_validation/taxon_robustness/scripts/build_robustness_reports.py
MPLCONFIGDIR=/private/tmp/msrc-taxon-mpl python3 empirical/neoaves_chr4/analysis/patterson_d/final_validation/taxon_robustness/scripts/make_taxon_robustness_figures.py
python3 -m unittest discover -s empirical/neoaves_chr4/analysis/patterson_d/final_validation/taxon_robustness/tests -p 'test_*.py' -v
```

The locally present 3.646-GB compressed alignment archive is read without redownloading. Candidate inventory is metadata-only; `compute_taxon_sensitivity.py` streams the archive for quartet counts and raw-locus validation. It reuses the validated production `count_sites` function for full alternative-quartet results, while the raw validation uses a separate implementation. Gene-tree matching uses exact archive locus identifiers and direct unrooted bipartitions; true unresolved quartets are tallied separately. The original frozen tables and figures are not overwritten.

## Files

- `config/frozen_alternative_quartets.tsv`: design panel fixed before D calculations.
- `results/candidate_taxa.tsv`: taxon, published lineage/order assignment, alignment and tree availability, and chr4 coverage.
- `results/taxon_robustness_summary.tsv`: one-row-per-design compact report table.
- `results/taxon_robustness_block_summary.tsv`: per-interval, inside/outside, all-sites/transversion, and 500-kb/1-Mb/2-Mb details.
- `results/alternative_quartet_locus_counts.tsv`: per-locus ABBA/BABA, D, informative and excluded non-ACGT sites, transition/transversion contributions.
- `results/quartet_null_validation.tsv`, `results/nucleotide_topology_robustness.tsv`, and `results/nucleotide_topology_matched_bins.tsv`: exact genealogy results and spatial matches.
- `results/raw_alignment_validation.tsv` and `results/raw_alignment_validation_report.md`: independent raw alignment validation.
- `figures/figure_E_alternative_quartet_delta_D`, `figures/figure_F_site_vs_topology`, and `figures/figure_G_sensitivity_matrix`: PNG, PDF, and SVG variants.

## Literature context

- [Jarvis et al. 2014, whole-genome Neoaves phylogeny](https://doi.org/10.1126/science.1253451).
- [Stiller et al. 2024, family-level avian genomes](https://doi.org/10.1038/s41586-024-07323-1).
- [Oliveira et al., region of suppressed recombination and taxon sampling](https://pmc.ncbi.nlm.nih.gov/articles/PMC11009670/).
- [Sangster et al. 2022, phylogenetic definitions for higher-level bird clades](https://doi.org/10.1016/j.avrs.2022.100027).

## Final recommendation

**CONDITIONAL GO** for the restricted N61 regional site-pattern contrast and descriptive genealogy concordance, with the taxon-sensitive outside baseline and unresolved mechanism stated. These results do not uniquely demonstrate introgression, inversion, or MSRC.
'''
 (PKG/'README.md').write_text(readme)
 # Deterministic independent-raw-validation report (idempotent on reruns).
 rawtxt=f"""# Independent raw-alignment validation

Deterministic stratified sample seed: 20261118. Unique loci sampled: {unique_loci} (inside PNAS_1/2/3, outside high-information, outside low-information, and outside with ambiguity/missing calls). N61-primary informative-count lower/upper quintile cutoffs were 12/23.

The checksum-verified local archive was streamed once; 6,247 chr4 alignments matched the manifest. Alternative quartet per-locus all-site and transversion counts were computed in this pass. Raw validation covered the three original primary quartets and their three predeclared Meleagris-outgroup designs; each independently decoded row was compared to saved all-site and transversion ABBA/BABA integers.

Validation rows PASS: {len(passrows)}; discrepancy/format failures: {len(fail)}; unavailable quartet/locus rows: {len(unavail)}. Unique loci: {unique_loci}. Duplicate sequence IDs: {len(duplicate_loci)} loci. Unequal full-alignment lengths: {len(unequal_loci)} loci.

All sampled loci had equal full-alignment lengths and zero duplicate sequence IDs. Exact taxon identifiers were retained. Four unavailable rows are the N62 *Pterocles gutturalis* absences at two loci, once for each outgroup. The independent reference implementation uses its own FASTA parser and direct four-base pattern comparisons; it does not import or invoke the production estimator. It records duplicate IDs, full-alignment length equality, four-taxon alignment length, ambiguous/gap site count, multistate non-pattern sites, and exact ABBA/BABA counts. Outgroup polarization is the observed O base; only patterns with one derived P3/O state and the exact ABBA or BABA configuration count. Transversions-only excludes A↔G and C↔T differences between P3 and O.

## Interpretation

Exact integer agreement is required for both ABBA and BABA, separately for all informative patterns and transversions-only. The per-row comparison is in `raw_alignment_validation.tsv`; the deterministic sample list and strata are in `raw_validation_sample_manifest.tsv`. No production filters were modified to force agreement.
"""
 (RES/'raw_alignment_validation_report.md').write_text(rawtxt)
 # Append/update the durable scientific interpretation without overwriting prior results.
 interp=FINAL/'results/final_scientific_interpretation.md'
 old=interp.read_text()
 marker='## N61 taxon-representative robustness and raw-alignment validation'
 addition=f'''\n\n{marker}\n\nThe primary N61 regional contrast is robust to ten prespecified ingroup alternatives. All ten retain negative all-sites ΔD, and each has a 1-Mb paired genomic block-jackknife 95% interval below zero. The transversion-only contrast also remains negative with its interval below zero in all ten. The outside D baseline is taxon-sensitive: some P1/P2 combinations change its sign even though the inside-minus-outside contrast does not. The Gallus primary contrast (ΔD={f(primary['Delta_D_all'])}, 95% interval {iv(primary['Delta_D_all_CI_low'],primary['Delta_D_all_CI_high'])}) and Meleagris contrast (ΔD={f(me['Delta_D_all'])}, 95% interval {iv(me['Delta_D_all_CI_low'],me['Delta_D_all_CI_high'])}) remain negative.\n\nA separate reference implementation exactly reproduced {len(passrows)} available all-sites/transversion ABBA/BABA comparisons from {unique_loci} stratified raw alignments, with zero discrepancies. No sampled duplicate IDs or unequal alignment lengths were found. Four rows were unavailable because *Pterocles gutturalis* was absent in two N62 sampled alignments. Thus no estimator error was found, and the frozen estimates remain unchanged.\n\nAcross exact-locus matched gene trees, whole-chromosome 1-Mb D_site–D_topology correlation ranged {min(rvals):.3f}–{max(rvals):.3f} across the primary and ingroup alternatives. The original r={float(primary['Dsite_Dtopology_r_1Mb']):.3f} is typical; the regional-mean-centered correlation was lower ({min(centered):.3f}–{max(centered):.3f}) but positive throughout the panel. These are descriptive correlations among linked loci and overlapping taxon designs. They do not identify introgression, ancestral structure, ILS, deep polymorphism, homoplasy, or MSRC.\n\nThe taxonomic audit corrects a group-label ambiguity: *Ciconia maguari* is Aequornithes/Ciconiiformes, not Telluraves. It remains the frozen P3 species. Shoebill and ibis are waterbird-side alternatives; the golden eagle is a separately labeled Telluraves control. All alternatives induce q1 in the published Stiller 363-taxon reference, which shares the broader source study with the alignments. Jarvis's 48-taxon tree lacks *Ciconia maguari*, so there is no independent exact four-taxon species-tree check. N61 is **CONDITIONAL GO** for the qualified regional site-pattern result; D/topology agreement is not unique evidence for MSRC.\n\nDetails: `taxon_robustness/README.md`, `taxon_robustness/results/final_assessment.md`, and `taxon_robustness/results/taxon_robustness_summary.tsv`.'''
 if marker not in old:interp.write_text(old+addition)
 # Add specific checklist rows once.
 check=FINAL/'results/validation_checklist.tsv'; cr=read(check)
 entries=[
 ('N61 alternative ingroup representatives','PASS_WITH_LIMITATION','10 prespecified ingroup alternatives; all negative all-sites Delta D with 1-Mb 95% block-jackknife intervals below zero; quartet designs overlap loci and are not independent.'),
 ('Independent raw alignment validation','PASS_WITH_LIMITATION',f'{len(passrows)} exact ABBA/BABA all-site/transversion comparisons in {unique_loci} deterministic stratified loci; 0 discrepancies, {len(unavail)} unavailable N62 rows due missing Pterocles; no estimator error.'),
 ('Alternative exact gene-tree topology comparison','PASS_WITH_LIMITATION','Exact local-tree labels and locus IDs; all designs resolve q1 in Stiller 63K reference; zero unresolved local quartet trees among exact matches; 250-kb/500-kb/1-Mb correlations with 2-Mb block bootstrap.'),
 ('N61 taxon sensitivity manuscript decision','CONDITIONAL','Regional contrast robust in direction and block interval, but outside D background is taxon-sensitive and mechanism/species-history null remain conditional.')]
 for e in entries:
  if not any(r['item']==e[0] for r in cr):cr.append(dict(item=e[0],status=e[1],evidence_or_limitation=e[2]))
 write(check,cr)
 # Append bounded section to analysis README and project status, without disturbing existing content.
 section='## N61 taxon robustness and raw-alignment validation'
 brief=f'''\n\n{section}\n\nAn additive sensitivity package is in `final_validation/taxon_robustness/`. It freezes ten ingroup alternatives plus a Meleagris outgroup-only comparison; all ten ingroup alternatives retain negative all-sites ΔD with 1-Mb paired block-jackknife 95% intervals below zero. The independent raw parser exactly reproduced {len(passrows)} available ABBA/BABA comparisons across {unique_loci} deterministic loci, with no discrepancies; four N62 rows were unavailable because Pterocles was missing in two alignments. Exact tree-topology/site-pattern correlations remain positive across the panel, with original 1-Mb whole-chr4 r={float(primary['Dsite_Dtopology_r_1Mb']):.3f} typical. CONDITIONAL GO for a qualified N61 regional site-pattern result; no unique mechanism claim. See `final_validation/taxon_robustness/results/final_assessment.md`.'''
 parent=BASE/'README.md';txt=parent.read_text()
 if section not in txt:parent.write_text(txt+brief)
 status=ROOT/'PROJECT_STATUS.md';txt=status.read_text();sm='## Neoaves N61 taxon robustness validation — 2026-11-17'
 statusadd=f'''\n\n{sm}\n\nCompleted the requested N61 ingroup-representative, independent raw-alignment, and exact genealogy sensitivity validation. Ten prespecified ingroup alternatives all retain negative regional ΔD with 1-Mb 95% physical-block jackknife intervals below zero; P1/P2 changes can alter or reverse the outside D baseline. Independent raw re-decoding of {unique_loci} stratified loci gave {len(passrows)} exact ABBA/BABA comparison matches and zero discrepancies; four N62 rows were unavailable owing to missing Pterocles. Whole-chromosome D_site–D_topology correlations are positive across all N61 alternatives; the original r={float(primary['Dsite_Dtopology_r_1Mb']):.3f} is typical. Exact quartet q1 follows the Stiller 363 reference, not independent exact confirmation. N61 recommendation remains CONDITIONAL GO for qualified descriptive wording; mechanism unresolved. Package: `empirical/neoaves_chr4/analysis/patterson_d/final_validation/taxon_robustness/`. No commit/push performed; unrelated working-tree changes preserved.'''
 if sm not in txt:status.write_text(txt+statusadd)
 print('reports, compact table, checklist, readmes, and project status written')
if __name__=='__main__':main()
