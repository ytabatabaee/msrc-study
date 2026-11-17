# N61 taxon robustness and raw-alignment validation

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

The locally present 3.646-GB compressed alignment archive is read without redownloading. Candidate inventory is metadata-only; `compute_taxon_sensitivity.py` streams the archive for quartet counts and raw-locus validation. It reuses the validated production `count_sites` function for full alternative-quartet results, while the raw validation uses a separate implementation. Quartets are q1=(P1,P2)|(P3,O), q2=(P1,P3)|(P2,O), and q3=(P2,P3)|(P1,O); D_topology=q3−q2, matching the ABBA−BABA sign orientation. Gene-tree matching uses exact archive locus identifiers and direct unrooted bipartitions; true unresolved quartets are tallied separately. The original frozen tables and figures are not overwritten.

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
