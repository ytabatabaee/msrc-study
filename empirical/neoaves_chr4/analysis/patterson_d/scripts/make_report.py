#!/usr/bin/env python3
"""Turn the executed Patterson-D tables into manuscript-oriented summaries."""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path

BASE=Path(__file__).resolve().parents[1]; RESULTS=BASE/"results"

def rows(results):
    with (results/"patterson_d_summary.tsv").open() as fh: return list(csv.DictReader(fh,delimiter="\t"))
def f(x): return float(x) if x not in ("",None) else None
def main(argv=None):
    global RESULTS
    ap=argparse.ArgumentParser(description="Build Patterson-D manuscript tables and report")
    ap.add_argument("--results-dir", type=Path, default=RESULTS)
    args=ap.parse_args(argv); RESULTS=args.results_dir
    rs=rows(RESULTS); primary=[r for r in rs if r["filter"]=="all" and r["block_size"]=="1000000" and r["design"].endswith("_primary")]
    all_designs=[r for r in rs if r["filter"]=="all" and r["block_size"]=="1000000"]
    table=[]
    fields=["dataset","quartet_definition","genomic_region","region_class","n_loci","informative_sites","ABBA_total","BABA_total","Patterson_D","standard_error","ci_low","ci_high","inside_minus_outside_D","inside_minus_outside_se","inside_minus_outside_ci_low","inside_minus_outside_ci_high"]
    for r in primary:
        d=f(r["D"]); se=f(r["jackknife_se"]); ci1=None if d is None or se is None else d-1.96*se; ci2=None if d is None or se is None else d+1.96*se
        diff=f(r["inside_minus_outside_D"]); dse=f(r["inside_minus_outside_se"])
        table.append({"dataset":"Stiller2024_63430_intergenic","quartet_definition":r["design"],"genomic_region":r["region"],"region_class":r["region_class"],"n_loci":r["n_loci"],"informative_sites":r["informative_sites"],"ABBA_total":r["ABBA"],"BABA_total":r["BABA"],"Patterson_D":r["D"],"standard_error":r["jackknife_se"],"ci_low":"" if ci1 is None else ci1,"ci_high":"" if ci2 is None else ci2,"inside_minus_outside_D":r["inside_minus_outside_D"],"inside_minus_outside_se":r["inside_minus_outside_se"],"inside_minus_outside_ci_low":"" if diff is None or dse is None else diff-1.96*dse,"inside_minus_outside_ci_high":"" if diff is None or dse is None else diff+1.96*dse})
    with (RESULTS/"manuscript_results_table.tsv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(table)
    # Keep the manuscript-facing primary table stable, while also exposing all
    # predeclared alternative-outgroup designs in an equivalent table.
    all_table=[]
    for r in all_designs:
        d=f(r["D"]); se=f(r["jackknife_se"]); diff=f(r["inside_minus_outside_D"]); dse=f(r["inside_minus_outside_se"])
        all_table.append({"dataset":"Stiller2024_63430_intergenic","quartet_definition":r["design"],"genomic_region":r["region"],"region_class":r["region_class"],"n_loci":r["n_loci"],"informative_sites":r["informative_sites"],"ABBA_total":r["ABBA"],"BABA_total":r["BABA"],"Patterson_D":r["D"],"standard_error":r["jackknife_se"],"ci_low":"" if d is None or se is None else d-1.96*se,"ci_high":"" if d is None or se is None else d+1.96*se,"inside_minus_outside_D":r["inside_minus_outside_D"],"inside_minus_outside_se":r["inside_minus_outside_se"],"inside_minus_outside_ci_low":"" if diff is None or dse is None else diff-1.96*dse,"inside_minus_outside_ci_high":"" if diff is None or dse is None else diff+1.96*dse})
    with (RESULTS/"manuscript_results_table_all_designs.tsv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(all_table)
    lines=["# Real nucleotide Patterson's D analysis: Neoaves chromosome 4","", "## Inputs and validation", "", "The original Stiller 63,430-intergenic-alignment archive was downloaded from the current ERDA path:", "`https://sid.erda.dk/share_redirect/ENhZODU9YE/01_alignments_and_gene_trees/intergenic_regions/63430.alns.tar.gz`", "", "SHA-256: `87a2a643083bb567022683d99b4f5446cf9d52b1da0b0d694bd52aca55f37e43`; compressed size 3,645,947,507 bytes. The Stage-4D manifest contains 63,430 loci and 6,247 chr4 loci. All 6,247 chr4 IDs matched alignment members exactly after normalizing only the archive's `.1k.start` and manifest `_1k_start` naming and removing the tree-only suffix. No positional matching was used.", "", "The three frozen PNAS intervals are the published GalGal4 window-start inclusive intervals 25,030,000--32,670,000, 33,510,000--34,470,000, and 44,130,000--56,810,000. The local locus coordinates are retained from the Stiller manifest; interval membership follows the existing repository's verified published predicate and should not be read as a nucleotide-resolution breakpoint lift.", "", "## Frozen designs", "", "Primary designs use one representative per published clade and *Gallus gallus* as outgroup: Columbea = `Columba_livia, Phoenicopterus_ruber, Ciconia_maguari`; N61 = `Tauraco_erythrolophus, Columba_livia, Ciconia_maguari`; N62 = `Columba_livia, Pterocles_gutturalis, Tauraco_erythrolophus`. `Meleagris_gallopavo` was frozen as a galliform sensitivity outgroup. These choices were written before reading nucleotide D values.", "", "## Results", "", "D was calculated directly from nucleotide site patterns with outgroup state A and derived state B: ABBA = A/B/B/A and BABA = B/A/B/A. Gaps, IUPAC ambiguity, and missing calls were excluded. Sites with more than two states or no ABBA/BABA pattern were uninformative. Regional D pools ABBA and BABA counts before taking the ratio. Uncertainty is a genomic block jackknife; the main table uses 1-Mb blocks and sensitivity tables retain 500-kb and 2-Mb results.", ""]
    for design in ("Columbea_primary","N61_primary","N62_primary"):
        inside=next(r for r in primary if r["design"]==design and r["region"]=="all_PNAS" and r["region_class"]=="inside")
        outside=next(r for r in primary if r["design"]==design and r["region"]=="all_PNAS" and r["region_class"]=="outside")
        diff=f(inside["inside_minus_outside_D"]); se=f(inside["inside_minus_outside_se"]); z=None if diff is None or not se else diff/se
        lines.append(f"* **{design.replace('_primary','')}**: inside D = {f(inside['D']):.4f} (SE {f(inside['jackknife_se']):.4f}; {inside['n_loci']} loci, {inside['informative_sites']} sites), outside D = {f(outside['D']):.4f} (SE {f(outside['jackknife_se']):.4f}; {outside['n_loci']} loci, {outside['informative_sites']} sites), inside-minus-outside = {diff:.4f} (SE {se:.4f}; descriptive z = {z:.2f}).")
    lines += ["", "The interval-specific results are in `manuscript_results_table.tsv`; `manuscript_results_table_all_designs.tsv` additionally includes the predeclared alternative-outgroup designs. Columbea is nearly unchanged after pooling all PNAS intervals; N61 shows a negative inside-minus-outside contrast; N62 is essentially unchanged. The three intervals are not equivalent in size or site count, especially PNAS_2, so the interval-specific estimates have different precision.", "", "## Sensitivity and interpretation", "", "The predeclared alternate galliform outgroup and transversion-only estimates are in `patterson_d_summary.tsv`. They are sensitivity analyses, not alternate post hoc sign choices. The signal is a nucleotide discordance statistic and is not by itself evidence for introgression, recombination suppression, or MSRC. It can also reflect incomplete lineage sorting, taxon sampling, ancestral polymorphism, alignment/quality effects, mutation-rate variation, or outgroup polarization error.", "", "The results provide a real, reproducible Patterson-D baseline for comparison with the existing local-tree tracks. They do not independently establish that the PNAS structural intervals caused the observed D differences, and no claim of causal MSRC support is made from this analysis alone.", "", "## Outputs", "", "* `results/input_audit.json`, `results/locus_alignment_mapping.tsv`, `results/patterson_d_locus_results.tsv`", "* `results/patterson_d_summary.tsv`, `results/patterson_d_1Mb_windows.tsv`, `results/manuscript_results_table.tsv`, `results/manuscript_results_table_all_designs.tsv`", "* `figures/patterson_d_spatial_chr4.{svg,png,pdf}`", "* `figures/patterson_d_inside_outside.{svg,png,pdf}`", ""]
    (RESULTS/"interpretation_report.md").write_text("\n".join(lines))
    print("wrote", RESULTS/"manuscript_results_table.tsv", RESULTS/"interpretation_report.md")
if __name__=="__main__": main()
