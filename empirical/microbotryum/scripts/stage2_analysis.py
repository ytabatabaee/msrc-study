#!/usr/bin/env python3
"""Stage 2 Microbotryum coordinate-recovery and cross-species audit.

Stage 2 reuses the Stage 1 quartet assignments.  It does not infer a
chromosome from topology or concatenate unordered MPbarcode50 contigs.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import random
import subprocess
import sys
from collections import defaultdict
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from empirical.microbotryum.scripts.stage0_audit import gtf_coordinates, parse_gene_rows

RAW = ROOT / "data/microbotryum/raw"
OUT = ROOT / "empirical/microbotryum/results"
FIG = ROOT / "empirical/microbotryum/figures"
PAIR_NAMES = {
    "cat|lag": "M. shykoffianum (C212)–M. lagerheimii",
    "cat|sup": "M. shykoffianum (C212)–M. superbum (1065)",
    "lag|sup": "M. lagerheimii–M. superbum",
}


def write_table(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def read_tsv(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def species_pair(quartet: str) -> str | None:
    taxa = quartet.split("|")
    species = sorted({
        "cat" if taxon.startswith("Mv-cat-") else "lag" if taxon.startswith("Mv-lag-") else "sup" if taxon.startswith("Mv-sup-") else "unknown"
        for taxon in taxa
    })
    if len(species) != 2 or "unknown" in species:
        return None
    return "|".join(species)


def coordinate_mapping(rows: list[dict], coord: dict[int, dict]) -> list[dict]:
    result = []
    for row in rows:
        item = coord.get(row["gene_number"])
        if item:
            contigs = sorted(item["contigs"])
            result.append({
                "gene_tree_id": row["gene_label"], "gene_number": row["gene_number"], "deposited_stratum": row["stratum"],
                "tree_taxa": ";".join(sorted(row["taxa"])), "orthogroup_id": "not retained in gene_trees",
                "annotation": "MPbarcode50_tsebra.gtf", "annotation_gene_namespace": f"Msup_*_g{row['gene_number']}",
                "contig": ";".join(contigs), "start": item["start"], "end": item["end"],
                "coordinate_type": "contig-local only", "chromosome": "NA", "gene_order_rank": "NA",
                "stratum_assignment_status": "deposited label retained; physical boundary not independently mapped",
                "mapping_confidence": "provisional_identifier_match",
                "provenance": "Dryad gene_trees + MPbarcode50_tsebra.gtf; numeric suffix match only",
            })
        else:
            result.append({
                "gene_tree_id": row["gene_label"], "gene_number": row["gene_number"], "deposited_stratum": row["stratum"],
                "tree_taxa": ";".join(sorted(row["taxa"])), "orthogroup_id": "not retained in gene_trees",
                "annotation": "MPbarcode50_tsebra.gtf", "annotation_gene_namespace": "NA", "contig": "NA", "start": "NA", "end": "NA",
                "coordinate_type": "unmapped", "chromosome": "NA", "gene_order_rank": "NA",
                "stratum_assignment_status": "deposited label retained", "mapping_confidence": "unmapped", "provenance": "Dryad gene_trees + MPbarcode50_tsebra.gtf",
            })
    return result


def cross_species_summary(assignments: list[dict]) -> list[dict]:
    groups = defaultdict(list)
    for row in assignments:
        if row["meaningful_species_mating_comparison"] == "true":
            pair = species_pair(row["quartet"])
            if pair:
                groups[(row["stratum"], pair)].append(row)
    output = []
    for (stratum, pair), group in sorted(groups.items()):
        n = len(group); q1 = sum(row["topology"] == "q1_species" for row in group) / n; q2 = sum(row["topology"] == "q2_mating_type" for row in group) / n; q3 = sum(row["topology"] == "q3_third" for row in group) / n
        output.append({"stratum": stratum, "species_pair": pair, "species_pair_label": PAIR_NAMES[pair], "n_trees": n // 1, "q1_species": q1, "q2_mating_type": q2, "q3_third": q3, "D_q2_minus_q1": q2 - q1, "mean_B": sum(float(row["B_mean_pairwise_patristic"]) for row in group) / n, "mean_internal_branch": sum(float(row["internal_branch_length"]) for row in group if row["internal_branch_length"] != "NA") / sum(row["internal_branch_length"] != "NA" for row in group)})
    return output


def bootstrap_cross_species(assignments: list[dict], replicates: int = 1000) -> list[dict]:
    rng = random.Random(20261108); groups = defaultdict(list)
    for row in assignments:
        if row["meaningful_species_mating_comparison"] == "true":
            pair = species_pair(row["quartet"])
            if pair: groups[(row["stratum"], pair)].append(row)
    result = []
    for (stratum, pair), group in sorted(groups.items()):
        stats = {key: [] for key in ("q1_species", "q2_mating_type", "q3_third", "D_q2_minus_q1", "B", "internal")}
        for _ in range(replicates):
            sample = [group[rng.randrange(len(group))] for _ in group]
            q1 = sum(row["topology"] == "q1_species" for row in sample) / len(sample); q2 = sum(row["topology"] == "q2_mating_type" for row in sample) / len(sample); q3 = sum(row["topology"] == "q3_third" for row in sample) / len(sample)
            stats["q1_species"].append(q1); stats["q2_mating_type"].append(q2); stats["q3_third"].append(q3); stats["D_q2_minus_q1"].append(q2 - q1)
            stats["B"].append(sum(float(row["B_mean_pairwise_patristic"]) for row in sample) / len(sample))
            internal = [float(row["internal_branch_length"]) for row in sample if row["internal_branch_length"] != "NA"]; stats["internal"].append(sum(internal) / len(internal))
        for statistic, values in stats.items():
            ordered = sorted(values); estimate = sum(values) / len(values)
            result.append({"stratum": stratum, "species_pair": pair, "species_pair_label": PAIR_NAMES[pair], "n_trees": len(group), "statistic": statistic, "estimate": estimate, "ci_low": ordered[int(.025 * replicates)], "ci_high": ordered[int(.975 * replicates) - 1], "bootstrap_replicates": replicates, "note": "gene-level descriptive bootstrap; physical linkage remains"})
    return result


def stratum_order_metadata(summary: list[dict]) -> list[dict]:
    descriptions = {
        "PARs": ("pseudo-autosomal/recombining control", "control", "small PARs at chromosome edges; categorical evidence"),
        "black": ("ancestral shared nonrecombining stratum", "older; predates M. superbum–M. shykoffianum divergence", "shared stratum; published trans-specific-polymorphism interpretation"),
        "turquoise": ("M. superbum-specific young extension", "younger; post-speciation", "young stratum inferred from dS/changepoint and lack of trans-specific polymorphism"),
        "dark_ruby": ("young species-specific ruby sub-stratum", "younger; exact relative order not converted here", "published dS/changepoint interpretation"),
        "light_ruby": ("young species-specific ruby sub-stratum", "younger; exact relative order not converted here", "published dS/changepoint interpretation"),
        "blue": ("ancient stratum near HD locus", "ancient/shared", "published ancestral stratum"),
        "orange": ("ancient shared stratum", "ancient/shared", "published ancestral stratum"),
        "purple": ("ancient stratum near PR locus", "ancient/shared", "published ancestral stratum"),
    }
    result = []
    counts = {row["stratum"]: row["n_trees"] for row in summary}
    for stratum in sorted(counts):
        desc, age, evidence = descriptions.get(stratum, ("author-labeled stratum", "not independently assigned", "no recovered physical boundary"))
        result.append({"stratum": stratum, "biological_interpretation": desc, "published_age_class": age, "published_physical_order": "not recovered as MPbarcode50 chromosome coordinates", "gene_order_proxy": "M. lagerheimii ancestral order described in paper; not present as a retained mapping table", "recombination_regime": "recombining control" if stratum == "PARs" else "suppressed/stratum-labeled", "n_tree_records": counts[stratum], "boundary_status": "not independently mapped", "evidence": evidence, "reference": "Lucotte et al. 2025 DOI 10.1038/s41467-025-60222-5"})
    return result


def boundary_summary(summary: list[dict]) -> list[dict]:
    by = {(row["stratum"], row["species_pair"]): row for row in summary}
    result = []
    for pair in PAIR_NAMES:
        for left, right in (("PARs", "black"), ("black", "turquoise")):
            a, b = by[(left, pair)], by[(right, pair)]
            result.append({"species_pair": pair, "left_stratum": left, "right_stratum": right, "boundary_status": "not_spatially_testable; no validated coordinate boundary", "left_n": a["n_trees"], "right_n": b["n_trees"], "delta_q1": float(b["q1_species"]) - float(a["q1_species"]), "delta_q2": float(b["q2_mating_type"]) - float(a["q2_mating_type"]), "delta_D": float(b["D_q2_minus_q1"]) - float(a["D_q2_minus_q1"]), "delta_B": float(b["mean_B"]) - float(a["mean_B"]), "delta_internal": float(b["mean_internal_branch"]) - float(a["mean_internal_branch"]), "interpretation": "stratum contrast only; not a detected genomic change point"})
    return result


def write_coordinate_report(mapping: list[dict]) -> None:
    chrom = sum(row["chromosome"] != "NA" for row in mapping); contig = sum(row["coordinate_type"] == "contig-local only" for row in mapping); unmapped = sum(row["coordinate_type"] == "unmapped" for row in mapping); orthogroup = sum(row["orthogroup_id"] != "not retained in gene_trees" for row in mapping)
    lines = ["# Microbotryum Stage 2 coordinate validation", "", "## Mapping result", "", f"The deposited two-column `gene_trees` file contains 260 records with prefixes such as `black-g2118`; the Newick tips retain only six taxon labels. No orthogroup identifier is retained in these records. The local `Orthogroups_filtered.txt` contains full sample-specific gene identifiers and orthogroup rows, but there is no deposited table linking the `gNNNN` tree prefix to an `OG` row. Therefore orthogroup recovery is **not verified**.", "", f"The numeric suffixes in all {len(mapping)} tree records map to an `Msup_*_gNNN` identifier in `MPbarcode50_tsebra.gtf`; this is an identifier-suffix match, not proof of orthology or sample identity. Mapping counts are: {len(mapping)-unmapped} contig-local provisional matches, {chrom} chromosome-scale coordinates, {orthogroup} orthogroup-resolved records, and {unmapped} unmapped records.", "", "## Assembly/sample validation", "", "The paper identifies the focal phased assemblies as M. superbum strain 1065 and M. shykoffianum strain c212. The local README describes `MPbarcode50.fa` only as a haploid Flye ONT assembly; its anonymous `contig_*` headers and `Msup_tig...` annotation namespace are compatible with an M. superbum annotation, but do not identify the isolate or mating type. The article describes a submitted de novo annotation for M. superbum strain 6P, whereas the gene-tree taxa are Mv-sup-1065-A1/A2 and Mv-cat-C212-A1/A2. This prevents treating MPbarcode50 as the 1065 tree assembly without an explicit submitted-file cross-reference.", "", "## Physical order and boundaries", "", "The publication uses M. lagerheimii as a proxy for ancestral mating-type-chromosome gene order and documents the PR chromosome fusion, edge PARs, and colored strata. It does not provide a chromosome-coordinate transform in the local deposited files. The current MPbarcode50 matches span several anonymous contigs, including multiple strata sharing contigs. No chromosome-scale order, validated gene-order rank, or base-pair stratum boundary was recovered. No spatial genealogy tracks are generated and no contigs are concatenated.", "", "## Consequence", "", "The Stage 2 spatial question remains blocked. The valid empirical comparison is stratum-level, with the deposited stratum label as the independent regime annotation and the PAR category as a categorical control. The complete one-row-per-tree mapping table is retained for future resolution.", "", "References: [Lucotte et al. 2025](https://doi.org/10.1038/s41467-025-60222-5), [Dryad](https://doi.org/10.5061/dryad.0k6djhb89), [authors' repository](https://github.com/eliselucotte/MicrobotryumOnDianthus)."]
    (OUT / "stage2_coordinate_validation_report.md").write_text("\n".join(lines) + "\n")


def convert(svg: Path) -> None:
    for fmt in ("png", "pdf"):
        try: subprocess.run(["rsvg-convert", "-f", fmt, "-o", str(svg.with_suffix("." + fmt)), str(svg)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (OSError, subprocess.CalledProcessError): pass


def make_figures(cross: list[dict], order: list[dict]) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    strata = [row["stratum"] for row in order]; pairs = list(PAIR_NAMES)
    width, height, left, top, right, bottom = 1450, 820, 150, 90, 40, 190; plot_w, plot_h = width-left-right, height-top-bottom; xstep = plot_w / max(1, len(strata)-1); colors = {"cat|lag":"#4c78a8", "cat|sup":"#f58518", "lag|sup":"#54a24b"}
    def y(v): return top + plot_h * (1-v)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">','<rect width="100%" height="100%" fill="white"/>',f'<text x="{width/2}" y="34" text-anchor="middle" font-family="sans-serif" font-size="22">Cross-species-pair mating-type support by stratum</text>', '<text x="25" y="400" transform="rotate(-90 25 400)" text-anchor="middle" font-family="sans-serif" font-size="14">q2 frequency</text>']
    for j, stratum in enumerate(strata):
        x=left+j*xstep; parts.append(f'<text x="{x}" y="{top+plot_h+28}" transform="rotate(45 {x} {top+plot_h+28})" text-anchor="start" font-family="sans-serif" font-size="12">{escape(stratum)}</text>')
    for pair in pairs:
        vals=[]
        for stratum in strata:
            item=next(row for row in cross if row["stratum"]==stratum and row["species_pair"]==pair); vals.append(float(item["q2_mating_type"]))
        points=[]
        for j,v in enumerate(vals):
            x=left+j*xstep; points.append(f"{x},{y(v)}"); parts.append(f'<circle cx="{x}" cy="{y(v)}" r="5" fill="{colors[pair]}"/>')
        parts.append(f'<polyline points="{" ".join(points)}" fill="none" stroke="{colors[pair]}" stroke-width="3"/>')
    for j,pair in enumerate(pairs):
        x=left+j*360; parts.append(f'<line x1="{x}" y1="58" x2="{x+18}" y2="58" stroke="{colors[pair]}" stroke-width="3"/><text x="{x+25}" y="63" font-family="sans-serif" font-size="13">{escape(PAIR_NAMES[pair])}</text>')
    parts.append("</svg>"); svg=FIG/"microbotryum_stage2_cross_species_quartets.svg"; svg.write_text("\n".join(parts)); convert(svg)
    # Evidence-status organization diagram, explicitly not a physical map.
    parts=['<svg xmlns="http://www.w3.org/2000/svg" width="1450" height="650" viewBox="0 0 1450 650">','<rect width="100%" height="100%" fill="white"/>','<text x="725" y="35" text-anchor="middle" font-family="sans-serif" font-size="22">Documented Microbotryum mating-type organization and mapping status</text>','<text x="725" y="63" text-anchor="middle" font-family="sans-serif" font-size="13">Gene-order schematic only; no MPbarcode50 chromosome coordinates or boundaries are asserted.</text>']
    boxes=[("PAR", "recombining edge regions", "published; categorical control", "#bdbdbd"),("turquoise", "young species-specific extension", "published age class; position not recovered", "#17becf"),("light/dark ruby", "young species-specific strata", "published; exact order not recovered", "#8c564b"),("black", "older shared stratum", "predates superbum–shykoffianum split", "#222222"),("blue/orange/purple", "ancient shared strata", "published ancestral strata", "#6b6ecf")]
    y=115
    for name,desc,status,color in boxes:
        parts += [f'<rect x="130" y="{y}" width="280" height="64" rx="8" fill="{color}" opacity="0.9"/><text x="270" y="{y+26}" text-anchor="middle" font-family="sans-serif" font-size="16" fill="white">{escape(name)}</text>',f'<text x="450" y="{y+24}" font-family="sans-serif" font-size="15">{escape(desc)}</text>',f'<text x="450" y="{y+46}" font-family="sans-serif" font-size="13">{escape(status)}</text>']; y+=85
    parts += ['<line x1="90" y1="115" x2="90" y2="535" stroke="#555" stroke-width="3"/><text x="270" y="570" text-anchor="middle" font-family="sans-serif" font-size="13">published regime labels, not physical scale</text>','<text x="950" y="180" font-family="sans-serif" font-size="16" font-weight="bold">Validated in this repository</text>','<text x="950" y="215" font-family="sans-serif" font-size="14">260/260 tree records parsed</text>','<text x="950" y="245" font-family="sans-serif" font-size="14">260/260 provisional contig-ID matches</text>','<text x="950" y="275" font-family="sans-serif" font-size="14">0 chromosome-scale positions</text>','<text x="950" y="305" font-family="sans-serif" font-size="14">0 recovered orthogroup-to-tree links</text>','<text x="950" y="350" font-family="sans-serif" font-size="16" font-weight="bold">Published structural context</text>','<text x="950" y="385" font-family="sans-serif" font-size="14">PR chromosome fused with part of ancestral HD</text>','<text x="950" y="415" font-family="sans-serif" font-size="14">small PARs at chromosome edges</text>','<text x="950" y="445" font-family="sans-serif" font-size="14">M. lagerheimii used as ancestral-order proxy</text>','</svg>']; svg=FIG/"microbotryum_stage2_genomic_organization_status.svg"; svg.write_text("\n".join(parts)); convert(svg)


def write_report(mapping: list[dict], cross: list[dict], boundaries: list[dict], order: list[dict]) -> None:
    black = next(row for row in cross if row["stratum"]=="black" and row["species_pair"]=="cat|sup"); pars = next(row for row in cross if row["stratum"]=="PARs" and row["species_pair"]=="cat|sup"); turquoise = next(row for row in cross if row["stratum"]=="turquoise" and row["species_pair"]=="cat|sup")
    lines=["# Microbotryum Stage 2 analysis report","","## Recommendation","","**CONDITIONAL GO — strong but taxon-pair-specific stratum-level evidence, with spatial validation incomplete.** Stage 2 recovered no defensible chromosome-scale coordinates or orthogroup-to-tree mapping. It therefore does not claim a spatial change point or generate chromosome tracks.","","## Assembly and source validation","", "The paper identifies the phased biological assemblies as M. superbum strain 1065 and M. shykoffianum strain c212. The local `MPbarcode50.fa` is described only as a haploid Flye ONT assembly with anonymous contigs, and its GTF uses an `Msup_tig...` namespace. The paper separately describes a submitted M. superbum strain 6P annotation. The tree tips are Mv-sup-1065 and Mv-cat-C212. Because no local file links MPbarcode50 to 1065, 6P, or a mating type, the GTF match is provisional.","","The source paper assigns strain C212 to M. shykoffianum based on barcode genes, despite the deposited tree prefix `Mv-cat`. Stage 2 therefore labels the C212 taxon as M. shykoffianum (C212) and retains `Mv-cat-C212` as the literal sample identifier. The `Mv-cat` prefix is not treated as an independent species assignment.","","The authors' paper documents a PR chromosome fused with part of the ancestral HD chromosome, small PARs at chromosome edges, and M. lagerheimii as a proxy for ancestral gene order. It also describes black as shared/older and turquoise/ruby as young species-specific strata. These are independent biological annotations, but the submitted files currently available here do not provide a chromosome-coordinate transform.","","## Complete tree-to-coordinate audit","",f"The one-row-per-tree mapping contains {len(mapping)} records: {sum(row['coordinate_type']=='contig-local only' for row in mapping)} provisional contig-local matches, {sum(row['coordinate_type']=='unmapped' for row in mapping)} unmapped records, 0 chromosome-scale coordinates, 0 recovered orthogroup IDs, and 260 deposited stratum labels. The numeric suffix match to MPbarcode50 is not treated as proof of biological orthology.","","The six taxa are retained as Mv-cat-C212, Mv-lag-1253, and Mv-sup-1065, each with A1/A2 labels. The published biological assignment used in this report is C212 = M. shykoffianum, 1065 = M. superbum, and the `Mv-lag` sample as M. lagerheimii; sample prefixes remain preserved in all machine-readable tables.","","## Cross-species-pair results","","The three complete species-pair quartets were analyzed from the existing Stage 1 assignments. q1 is species clustering and q2 is mating-type clustering, with fixed labels across strata.","","| stratum | species pair | n | q1 | q2 | q3 | D | B | internal |","|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for row in cross:
        lines.append(f"| {row['stratum']} | {row['species_pair']} | {row['n_trees']} | {float(row['q1_species']):.3f} | {float(row['q2_mating_type']):.3f} | {float(row['q3_third']):.3f} | {float(row['D_q2_minus_q1']):.3f} | {float(row['mean_B']):.5f} | {float(row['mean_internal_branch']):.5f} |")
    lines += ["",f"For the primary C212–1065 pair, q2 changes from {float(pars['q2_mating_type']):.3f} in PAR to {float(black['q2_mating_type']):.3f} in black and {float(turquoise['q2_mating_type']):.3f} in turquoise. This reversal is **not shared** by the other two pairings: C212–lagerheimii and lagerheimii–1065 remain q1 = 1.000 in black. The black signal is therefore taxon-pair-specific in these trees, not evidence of a uniform three-species reversal.","","This restriction weakens a broad ancestral-process interpretation and makes introgression, species-specific history, or taxon-specific gene-tree effects important alternatives. It remains a stratum contrast rather than a boundary-localized test.","","## Boundary and age analysis","", "PAR–black and black–turquoise effects are reported in `stage2_boundary_transition_summary.tsv`, but every row is marked `not_spatially_testable`: no validated physical or gene-order boundary was recovered. These deltas must not be interpreted as topology-derived boundary discovery.","","The age pattern is compatible with the published chronology: black is older/shared and turquoise is younger/species-specific; the paper reports turquoise and ruby as lower-dS young strata relative to black and identifies light/dark ruby sub-strata. The gene-tree data show the strongest mating-type support in black only for the C212–1065 pair, not a monotonic or universal age trend. Small strata such as orange, purple, and red are too sparsely sampled for strong age-trend claims.","","## Spatial analysis status","","No q1/q2/q3, D, A, or B chromosome tracks were generated. The available GTF matches span multiple anonymous contigs and the tree prefix cannot be linked to an orthogroup or published ancestral gene-order rank. Contig-local plotting would not establish adjacency of strata, so the spatial analysis remains blocked.","","## Biological interpretation and limitations","","The C212–1065 black-stratum reversal is compatible with an older mating-type-associated genealogy maintained under recombination suppression, but the lack of the same reversal in the other two species-pair quartets is a major limitation. Selection, gene conversion, introgression, branch-length heterogeneity, alignment quality, and gene-tree estimation error remain alternatives. In particular, the fungal mating-type chromosome history includes chromosome fusion, rearrangement, selfing/intra-tetrad mating, and loss of HD mating function, which differs from an idealized two-background model.","","## Files and next blocker","","- `stage2_coordinate_mapping.tsv`: complete one-row-per-tree mapping audit.","- `stage2_stratum_order_metadata.tsv`: published regime and age annotations without invented coordinates.","- `cross_species_pair_summary.tsv` and `cross_species_pair_bootstrap.tsv`: all three complete species-pair quartet comparisons.","- `stage2_boundary_transition_summary.tsv`: descriptive stratum contrasts explicitly marked non-spatial.","- `microbotryum_stage2_cross_species_quartets.*`: cross-species quartet figure.","- `microbotryum_stage2_genomic_organization_status.*`: evidence-status organization diagram, not a physical map.","","The concrete next step is to obtain the authors' source-data/table containing the M. lagerheimii ancestral-order ranks and the exact tree-prefix-to-orthogroup conversion. With those files, the existing assignments can be projected onto gene-order bins without reconstructing trees.","","References: [Lucotte et al. 2025](https://doi.org/10.1038/s41467-025-60222-5), [Dryad](https://doi.org/10.5061/dryad.0k6djhb89), [authors' repository](https://github.com/eliselucotte/MicrobotryumOnDianthus)."]
    (OUT/"stage2_analysis_report.md").write_text("\n".join(lines)+"\n")


def main() -> None:
    parser=argparse.ArgumentParser(); parser.add_argument("--raw-dir",type=Path,default=RAW); parser.add_argument("--bootstrap-replicates",type=int,default=1000); args=parser.parse_args()
    rows=[row for row in parse_gene_rows(args.raw_dir/"gene_trees") if row["root"] is not None]
    assignments=read_tsv(OUT/"all_quartet_assignments.tsv"); support=read_tsv(OUT/"preliminary_quartet_support.tsv")
    coord=gtf_coordinates(args.raw_dir/"MPbarcode50_tsebra.gtf"); mapping=coordinate_mapping(rows,coord); cross=cross_species_summary(assignments); boot=bootstrap_cross_species(assignments,args.bootstrap_replicates); order=stratum_order_metadata(support); boundaries=boundary_summary(cross)
    write_table(OUT/"stage2_coordinate_mapping.tsv",mapping); write_table(OUT/"stage2_stratum_order_metadata.tsv",order); write_table(OUT/"cross_species_pair_summary.tsv",cross); write_table(OUT/"cross_species_pair_bootstrap.tsv",boot); write_table(OUT/"stage2_boundary_transition_summary.tsv",boundaries); write_coordinate_report(mapping); make_figures(cross,order); write_report(mapping,cross,boundaries,order)
    print(f"Stage 2: {len(mapping)} mapped tree records; {len(cross)} cross-species stratum summaries; chromosome tracks blocked")


if __name__ == "__main__": main()
