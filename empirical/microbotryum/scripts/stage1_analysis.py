#!/usr/bin/env python3
"""Stage 1 stratum-level Microbotryum analysis.

This script deliberately does not overwrite Stage 0 outputs.  It consumes the
local deposited gene trees and produces all-quartet assignments, descriptive
gene-level bootstrap intervals, stratum summaries, and figures.  Chromosome-
wide tracks are not produced because the current MPbarcode50 mapping is not
independently validated as a chromosome-scale mapping for the tree source.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import random
import subprocess
import sys
from collections import Counter, defaultdict
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

try:
    from empirical.microbotryum.scripts.stage0_audit import (
        FOCAL, OUT as STAGE0_OUT, RAW, distances, gtf_coordinates,
        internal_branch, parse_gene_rows, split_for_quartet,
    )
except ModuleNotFoundError:
    from stage0_audit import (
        FOCAL, OUT as STAGE0_OUT, RAW, distances, gtf_coordinates,
        internal_branch, parse_gene_rows, split_for_quartet,
    )

OUT = ROOT / "empirical/microbotryum/results"
FIG = ROOT / "empirical/microbotryum/figures"
QUARTET = tuple(FOCAL.values())
TAXA_ORDER = (
    "Mv-cat-C212-A1", "Mv-cat-C212-A2",
    "Mv-lag-1253-A1", "Mv-lag-1253-A2",
    "Mv-sup-1065-A1", "Mv-sup-1065-A2",
)
SPECIES = {
    "Mv-cat-C212-A1": "cat", "Mv-cat-C212-A2": "cat",
    "Mv-lag-1253-A1": "lag", "Mv-lag-1253-A2": "lag",
    "Mv-sup-1065-A1": "sup", "Mv-sup-1065-A2": "sup",
}
MATING = {taxon: taxon.rsplit("-", 1)[1] for taxon in TAXA_ORDER}
MEANINGFUL = {
    tuple(sorted(pair))
    for pair in itertools.combinations(("cat", "lag", "sup"), 2)
}


def write_table(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def resolutions(quartet: tuple[str, ...]) -> list[frozenset]:
    a, b, c, d = quartet
    return [
        frozenset((frozenset((a, b)), frozenset((c, d)))),
        frozenset((frozenset((a, c)), frozenset((b, d)))),
        frozenset((frozenset((a, d)), frozenset((b, c)))),
    ]


def resolution_index(split, quartet: tuple[str, ...]) -> int | None:
    if split is None:
        return None
    for index, candidate in enumerate(resolutions(quartet), 1):
        if split == candidate:
            return index
    return None


def meaningful_labels(quartet: tuple[str, ...]) -> tuple[bool, list[str]]:
    species = sorted({SPECIES[taxon] for taxon in quartet})
    meaningful = len(species) == 2 and all(sum(SPECIES[t] == species_name for t in quartet) == 2 for species_name in species)
    if not meaningful or tuple(species) not in MEANINGFUL:
        return False, ["resolution_1", "resolution_2", "resolution_3"]
    members = {species_name: sorted(t for t in quartet if SPECIES[t] == species_name) for species_name in species}
    species_split = frozenset((frozenset(members[species[0]]), frozenset(members[species[1]])))
    same_mating = frozenset((
        frozenset(t for t in quartet if MATING[t] == "A1"),
        frozenset(t for t in quartet if MATING[t] == "A2"),
    ))
    third = next(candidate for candidate in [*resolutions(quartet)] if candidate not in {species_split, same_mating})
    labels = []
    for candidate in resolutions(quartet):
        labels.append("q1_species" if candidate == species_split else "q2_mating_type" if candidate == same_mating else "q3_third" if candidate == third else "unresolved")
    return True, labels


def primary_split_labels(split) -> str:
    if split is None:
        return "unresolved"
    species = frozenset((frozenset((FOCAL["cat_a1"], FOCAL["cat_a2"])), frozenset((FOCAL["sup_a1"], FOCAL["sup_a2"]))))
    mating = frozenset((frozenset((FOCAL["cat_a1"], FOCAL["sup_a1"])), frozenset((FOCAL["cat_a2"], FOCAL["sup_a2"]))))
    third = frozenset((frozenset((FOCAL["cat_a1"], FOCAL["sup_a2"])), frozenset((FOCAL["cat_a2"], FOCAL["sup_a1"]))))
    return "q1_species" if split == species else "q2_mating_type" if split == mating else "q3_third" if split == third else "other"


def quartet_distance_stats(root: object, quartet: tuple[str, ...]) -> tuple[float, float | None]:
    weighted = distances(root, weighted=True)
    values = [weighted[tuple(sorted(pair))] for pair in itertools.combinations(quartet, 2)]
    split = split_for_quartet(root, quartet)
    return sum(values) / len(values), internal_branch(root, split)


def primary_metrics(root: object) -> dict[str, float | None]:
    weighted = distances(root, weighted=True)
    def pair(a: str, b: str) -> float:
        return weighted[tuple(sorted((a, b)))]
    same_mating = (pair(FOCAL["cat_a1"], FOCAL["sup_a1"]) + pair(FOCAL["cat_a2"], FOCAL["sup_a2"])) / 2
    opposite_mating = (pair(FOCAL["cat_a1"], FOCAL["sup_a2"]) + pair(FOCAL["cat_a2"], FOCAL["sup_a1"])) / 2
    b, internal = quartet_distance_stats(root, QUARTET)
    return {"B": b, "A_opposite_minus_same": opposite_mating - same_mating, "internal_branch": internal}


def all_assignments(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    assignment_rows: list[dict] = []
    summary_rows: list[dict] = []
    all_quartets = list(itertools.combinations(TAXA_ORDER, 4))
    groups: dict[tuple[str, tuple[str, ...]], list[dict]] = defaultdict(list)
    for row in rows:
        if row["root"] is None or not set(TAXA_ORDER).issubset(row["taxa"]):
            continue
        for quartet in all_quartets:
            split = split_for_quartet(row["root"], quartet)
            index = resolution_index(split, quartet)
            meaningful, labels = meaningful_labels(quartet)
            stats = quartet_distance_stats(row["root"], quartet)
            record = {
                "stratum": row["stratum"], "gene_label": row["gene_label"], "gene_number": row["gene_number"],
                "quartet": "|".join(quartet), "meaningful_species_mating_comparison": str(meaningful).lower(),
                "topology": labels[index - 1] if index else "unresolved", "resolution_index": index or "NA",
                "B_mean_pairwise_patristic": stats[0], "internal_branch_length": stats[1] if stats[1] is not None else "NA",
            }
            assignment_rows.append(record)
            groups[(row["stratum"], quartet)].append(record)
    for (stratum, quartet), group in sorted(groups.items()):
        meaningful = group[0]["meaningful_species_mating_comparison"] == "true"
        labels = ["q1_species", "q2_mating_type", "q3_third"] if meaningful else ["resolution_1", "resolution_2", "resolution_3"]
        counts = Counter(item["topology"] for item in group)
        numeric_internal = [float(item["internal_branch_length"]) for item in group if item["internal_branch_length"] != "NA"]
        summary_rows.append({
            "stratum": stratum, "quartet": "|".join(quartet), "meaningful_species_mating_comparison": str(meaningful).lower(),
            "n_trees": len(group), "q1_or_resolution1": counts[labels[0]] / len(group),
            "q2_or_resolution2": counts[labels[1]] / len(group), "q3_or_resolution3": counts[labels[2]] / len(group),
            "mean_B": sum(float(item["B_mean_pairwise_patristic"]) for item in group) / len(group),
            "mean_internal_branch": sum(numeric_internal) / len(numeric_internal) if numeric_internal else "NA",
        })
    return assignment_rows, summary_rows


def bootstrap(rows: list[dict], replicates: int = 1000) -> list[dict]:
    rng = random.Random(20261108)
    results = []
    for stratum in sorted({row["stratum"] for row in rows}):
        group = [row for row in rows if row["stratum"] == stratum and row["root"] is not None]
        values = {"q1_species": [], "q2_mating_type": [], "q3_third": [], "B": [], "A_opposite_minus_same": [], "internal_branch": []}
        labels = [primary_split_labels(split_for_quartet(row["root"], QUARTET)) for row in group]
        metrics = [primary_metrics(row["root"]) for row in group]
        for _ in range(replicates):
            sample = [rng.randrange(len(group)) for _ in group]
            for label in ("q1_species", "q2_mating_type", "q3_third"):
                values[label].append(sum(labels[i] == label for i in sample) / len(sample))
            for key in ("B", "A_opposite_minus_same", "internal_branch"):
                chosen = [metrics[i][key] for i in sample if metrics[i][key] is not None]
                values[key].append(sum(chosen) / len(chosen))
        estimate = {key: sum(value) / len(value) for key, value in values.items()}
        observed_labels = Counter(labels)
        observed_metrics = {key: [metric[key] for metric in metrics if metric[key] is not None] for key in ("B", "A_opposite_minus_same", "internal_branch")}
        observed = {"q1_species": observed_labels["q1_species"] / len(group), "q2_mating_type": observed_labels["q2_mating_type"] / len(group), "q3_third": observed_labels["q3_third"] / len(group)}
        observed.update({key: sum(vals) / len(vals) for key, vals in observed_metrics.items()})
        for key, samples in values.items():
            sorted_values = sorted(samples)
            results.append({"stratum": stratum, "n_trees": len(group), "statistic": key, "estimate": observed[key], "bootstrap_mean": estimate[key], "ci_low": sorted_values[int(0.025 * replicates)], "ci_high": sorted_values[int(0.975 * replicates) - 1], "bootstrap_replicates": replicates, "note": "gene-level descriptive bootstrap; linked genes are not independent physical replicates"})
    return results


def write_coordinate_report(rows: list[dict], coord: dict[int, dict]) -> None:
    ids = {row["gene_number"] for row in rows if row["root"] is not None}
    lines = [
        "# Microbotryum coordinate validation",
        "",
        "## Finding",
        "",
        f"All {len(ids)} unique gene identifiers in the {len(rows)} deposited gene-tree records map numerically to `MPbarcode50_tsebra.gtf` ({len(coord)} unique GTF identifiers). This is useful identifier continuity, but it is not by itself proof that the annotation and the deposited trees use the same isolate, mating type, assembly, or orthology release.",
        "",
        "The current files support **contig-local order only**. No chromosome-scale placement of MPbarcode50 contigs, no validated mapping to the M. superbum or M. shykoffianum mating-type chromosomes, and no independently verified base-pair boundaries for each color stratum were recovered from the local Dryad files or the inspected author scripts. No spatial tracks are therefore generated.",
        "",
        "## Assembly and sample identity",
        "",
        "`MPbarcode50.fa` has anonymous `contig_*` headers. Its annotation uses `Msup_tig..._g...` identifiers, indicating an M. superbum-oriented annotation namespace, but the local README does not state the exact isolate or mating type. The tree tips include `Mv-sup-1065-A1/A2`, `Mv-cat-C212-A1/A2`, and `Mv-lag-1253-A1/A2`; the MPbarcode50 sample correspondence remains unverified. The separate Mv-sup and Mv-cat FASTA resources use different sample-specific contig namespaces, so numerical gene-ID matching must not be treated as a chromosome mapping.",
        "",
        "The published study states that mating-type chromosome contigs were identified using PR/HD homology and ortholog links, and that PARs occur at chromosome edges. It also states that gene trees were inferred for one-to-one orthologs and deposited with the species tree. Those statements independently support the stratum labels and biological interpretation, but do not supply a validated MPbarcode50 contig order in the files currently available here.",
        "",
        "## Mapping confidence",
        "",
        "| component | status | confidence |",
        "|---|---|---|",
        "| gene-tree Newick parsing | all 260 parse with six taxa and branch lengths | high |",
        "| tree gene identifier to MPbarcode50 GTF identifier | 260/260 unique tree IDs represented | moderate for identifier continuity; low for biological sample correspondence |",
        "| physical coordinates within an MPbarcode50 contig | available from GTF | moderate, annotation-dependent |",
        "| chromosome-scale ordering | not recovered | none |",
        "| stratum-to-boundary base-pair mapping | not independently verified | none |",
        "",
        "## Consequence",
        "",
        "The valid analysis unit is evolutionary stratum, not a continuous chromosome window. The PARs category is used as a categorical recombining control. Contig-local coordinates are retained in `coordinate_mapping_audit.tsv` but are not concatenated or plotted as a chromosome.",
        "",
        "References: [Lucotte et al. 2025](https://doi.org/10.1038/s41467-025-60222-5), [Dryad](https://doi.org/10.5061/dryad.0k6djhb89), [authors' repository](https://github.com/eliselucotte/MicrobotryumOnDianthus).",
    ]
    (OUT / "coordinate_validation_report.md").write_text("\n".join(lines) + "\n")


def write_stratum_metadata(support: list[dict]) -> None:
    interpretations = {
        "PARs": ("pseudo-autosomal/recombining control", "control", "published PAR category; categorical control"),
        "black": ("older ancestral nonrecombining stratum with trans-specific polymorphism", "older", "published evolutionary-stratum interpretation"),
        "turquoise": ("younger species-specific nonrecombining stratum", "younger", "published evolutionary-stratum interpretation"),
        "in_HD": ("HD-associated subset", "unspecified", "author stratum label; no independent age assignment"),
        "HD_RR": ("HD/recombining-region associated subset", "unspecified", "author stratum label; no quantitative rate"),
    }
    rows = []
    for item in support:
        name = item["stratum"]
        interpretation, age, evidence = interpretations.get(name, ("color-coded evolutionary stratum; exact age not assigned here", "unspecified", "author stratum label; physical boundary not independently recovered"))
        rows.append({"stratum": name, "biological_interpretation": interpretation, "age_class": age, "physical_location": "not chromosome-validated; contig-local coordinates only", "recombination_evidence": evidence, "n_trees": item["n_trees"], "q1_species": item["q1_species_frequency"], "q2_mating_type": item["q2_mating_type_frequency"], "q3_third": item["q3_third_frequency"], "reference": "Lucotte et al. 2025 DOI 10.1038/s41467-025-60222-5"})
    write_table(OUT / "stratum_metadata.tsv", rows)


def svg_convert(svg_path: Path) -> None:
    for fmt in ("png", "pdf"):
        try:
            subprocess.run(["rsvg-convert", "-f", fmt, "-o", str(svg_path.with_suffix("." + fmt)), str(svg_path)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (OSError, subprocess.CalledProcessError):
            pass


def svg_text_plot(path: Path, title: str, series: list[tuple[str, str, list[float], list[tuple[float, float]] | None]], labels: list[str], y_label: str, y_min: float, y_max: float) -> None:
    width, height = 1300, 760; left, right, top, bottom = 100, 45, 75, 210; plot_w, plot_h = width - left - right, height - top - bottom
    x_step = plot_w / max(1, len(labels) - 1)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>', f'<text x="{width/2}" y="32" text-anchor="middle" font-family="sans-serif" font-size="22">{escape(title)}</text>', f'<text x="24" y="{top + plot_h/2}" transform="rotate(-90 24 {top + plot_h/2})" text-anchor="middle" font-family="sans-serif" font-size="14">{escape(y_label)}</text>', f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="black"/>', f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="black"/>']
    def y(value): return top + plot_h * (1 - (value - y_min) / (y_max - y_min or 1))
    for i, label in enumerate(labels):
        x = left + i * x_step; parts.append(f'<line x1="{x}" y1="{top + plot_h}" x2="{x}" y2="{top + plot_h + 5}" stroke="black"/>'); parts.append(f'<text x="{x}" y="{top + plot_h + 25}" transform="rotate(45 {x} {top + plot_h + 25})" text-anchor="start" font-family="sans-serif" font-size="12">{escape(label)}</text>')
    for name, color, values, intervals in series:
        points = []
        for i, value in enumerate(values):
            x = left + i * x_step; yy = y(value); points.append(f"{x},{yy}")
            if intervals:
                low, high = intervals[i]; parts.append(f'<line x1="{x}" y1="{y(low)}" x2="{x}" y2="{y(high)}" stroke="{color}" stroke-width="2"/><line x1="{x-4}" y1="{y(low)}" x2="{x+4}" y2="{y(low)}" stroke="{color}"/><line x1="{x-4}" y1="{y(high)}" x2="{x+4}" y2="{y(high)}" stroke="{color}"/>')
        parts.append(f'<polyline points="{" ".join(points)}" fill="none" stroke="{color}" stroke-width="3"/><circle cx="{left}" cy="{top}" r="0" fill="none"/>')
        for i, value in enumerate(values):
            parts.append(f'<circle cx="{left + i*x_step}" cy="{y(value)}" r="4" fill="{color}"/>')
    for i, (name, color, _, _) in enumerate(series):
        x = left + i * 225; parts.append(f'<line x1="{x}" y1="54" x2="{x+18}" y2="54" stroke="{color}" stroke-width="3"/><text x="{x+24}" y="59" font-family="sans-serif" font-size="13">{escape(name)}</text>')
    parts.append("</svg>"); path.write_text("\n".join(parts)); svg_convert(path)


def make_figures(support: list[dict], bootstrap_rows: list[dict], all_summary: list[dict], rows: list[dict]) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    strata = [item["stratum"] for item in support]; boot = {(item["stratum"], item["statistic"]): item for item in bootstrap_rows}
    series = []
    for key, name, color in [("q1_species", "q1 species", "#4c78a8"), ("q2_mating_type", "q2 mating type", "#f58518"), ("q3_third", "q3 third", "#54a24b")]:
        values = [next(float(item[k]) for item in support if item["stratum"] == s) for s, k in [(s, {"q1_species":"q1_species_frequency", "q2_mating_type":"q2_mating_type_frequency", "q3_third":"q3_third_frequency"}[key]) for s in strata]]
        intervals = [(boot[(s, key)]["ci_low"], boot[(s, key)]["ci_high"]) for s in strata]
        series.append((name, color, values, intervals))
    svg_text_plot(FIG / "microbotryum_stage1_topology_by_stratum.svg", "Microbotryum focal quartet topology by evolutionary stratum", series, strata, "frequency across genes", 0, 1)
    distance_series = []
    for key, name, color in [("B", "B mean pairwise distance", "#4c78a8"), ("A_opposite_minus_same", "A opposite - same mating type", "#f58518"), ("internal_branch", "quartet internal branch", "#54a24b")]:
        values = [boot[(s, key)]["estimate"] for s in strata]; intervals = [(boot[(s, key)]["ci_low"], boot[(s, key)]["ci_high"]) for s in strata]; distance_series.append((name, color, values, intervals))
    all_dist = [value for _, _, values, intervals in distance_series for value in values + [x for pair in intervals for x in pair]]; span = max(all_dist) - min(all_dist); svg_text_plot(FIG / "microbotryum_stage1_distance_by_stratum.svg", "Microbotryum focal quartet divergence by evolutionary stratum", distance_series, strata, "substitution-tree units", min(all_dist) - max(span * .1, .001), max(all_dist) + max(span * .1, .001))
    # Heatmap: fixed resolution-2 support for every quartet; meaningful rows are explicitly tagged.
    lookup = {(item["stratum"], item["quartet"]): item for item in all_summary}; quartets = sorted({item["quartet"] for item in all_summary})
    width, height, left, top, cell_w, cell_h = 1450, 900, 420, 100, 85, 34
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>', f'<text x="{width/2}" y="35" text-anchor="middle" font-family="sans-serif" font-size="22">Alternative-resolution support across all 15 quartets</text>', '<text x="20" y="60" font-family="sans-serif" font-size="12">Cells show fixed resolution-2 frequency; primary meaningful quartets are marked *</text>']
    for j, stratum in enumerate(strata): parts.append(f'<text x="{left+j*cell_w+cell_w/2}" y="82" transform="rotate(45 {left+j*cell_w+cell_w/2} 82)" text-anchor="start" font-family="sans-serif" font-size="12">{escape(stratum)}</text>')
    for i, quartet in enumerate(quartets):
        meaningful = next(item for item in all_summary if item["quartet"] == quartet)["meaningful_species_mating_comparison"] == "true"
        short_names = {"Mv-cat-C212-A1": "C1", "Mv-cat-C212-A2": "C2", "Mv-lag-1253-A1": "L1", "Mv-lag-1253-A2": "L2", "Mv-sup-1065-A1": "S1", "Mv-sup-1065-A2": "S2"}
        label = ("* " if meaningful else "") + "/".join(short_names[taxon] for taxon in quartet.split("|"))
        parts.append(f'<text x="{left-8}" y="{top+i*cell_h+22}" text-anchor="end" font-family="sans-serif" font-size="10">{escape(label)}</text>')
        for j, stratum in enumerate(strata):
            value = float(lookup[(stratum, quartet)]["q2_or_resolution2"]); shade = int(255 - 190 * value); color = f"rgb(255,{shade},{shade})"; parts.append(f'<rect x="{left+j*cell_w}" y="{top+i*cell_h}" width="{cell_w-2}" height="{cell_h-2}" fill="{color}" stroke="white"/><text x="{left+j*cell_w+cell_w/2}" y="{top+i*cell_h+20}" text-anchor="middle" font-family="sans-serif" font-size="10">{value:.2f}</text>')
    parts.append("</svg>"); heatmap = FIG / "microbotryum_stage1_all_quartet_heatmap.svg"; heatmap.write_text("\n".join(parts)); svg_convert(heatmap)
    # Representative trees: exact real Newick strings are included as text, avoiding claims of scaled branch drawing.
    selected = {}
    for row in rows:
        if row["stratum"] in {"PARs", "black", "turquoise"} and row["stratum"] not in selected:
            selected[row["stratum"]] = row
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="1500" height="620" viewBox="0 0 1500 620">', '<rect width="100%" height="100%" fill="white"/>', '<text x="750" y="32" text-anchor="middle" font-family="sans-serif" font-size="22">Representative deposited Microbotryum gene trees</text>', '<text x="30" y="58" font-family="sans-serif" font-size="12">Newick strings are reproduced verbatim; branch lengths are substitution-tree units and are not drawn to scale.</text>']
    y = 105
    for stratum in ("PARs", "black", "turquoise"):
        row = selected.get(stratum)
        if not row: continue
        label = f"{stratum} — {row['gene_label']} — {primary_split_labels(split_for_quartet(row['root'], QUARTET))}"; parts.append(f'<text x="30" y="{y}" font-family="sans-serif" font-size="15" font-weight="bold">{escape(label)}</text>')
        text = row["tree"]
        for chunk in [text[i:i+145] for i in range(0, len(text), 145)]:
            y += 18; parts.append(f'<text x="45" y="{y}" font-family="monospace" font-size="11">{escape(chunk)}</text>')
        y += 36
    parts.append("</svg>"); trees = FIG / "microbotryum_stage1_representative_trees.svg"; trees.write_text("\n".join(parts)); svg_convert(trees)


def write_report(support: list[dict], all_summary: list[dict], bootstrap_rows: list[dict], rows: list[dict]) -> None:
    primary = next(item for item in support if item["stratum"] == "black")
    turquoise = next(item for item in support if item["stratum"] == "turquoise")
    pars = next(item for item in support if item["stratum"] == "PARs")
    n_meaningful = len({item["quartet"] for item in all_summary if item["meaningful_species_mating_comparison"] == "true"})
    lines = [
        "# Microbotryum Stage 1 empirical analysis",
        "",
        "## Recommendation",
        "",
        "**CONDITIONAL GO — strong stratum-level signal, but incomplete spatial and assembly validation.** The deposited trees provide a compelling descriptive transition between the PAR control, the older black stratum, and the younger turquoise stratum. A chromosome-coordinate MSRC analysis is not yet justified because the MPbarcode50 annotation-to-tree sample/assembly correspondence and chromosome-scale contig order remain unvalidated.",
        "",
        "## Data and coordinate status",
        "",
        f"The analysis uses the local Dryad `gene_trees` file and {len(rows)} parseable records. Each record contains six taxa and complete branch lengths. The 260 trees span 11 deposited stratum labels. The focal quartet is `Mv-cat-C212-A1/A2` plus `Mv-sup-1065-A1/A2`; the third species pair is `Mv-lag-1253-A1/A2`.",
        "",
        "The separate `coordinate_validation_report.md` is the authoritative mapping assessment. All 260 tree gene IDs map numerically to `MPbarcode50_tsebra.gtf`, but the matching ID namespace does not establish the same sample or assembly. Contig-local coordinates are retained; contigs are not concatenated into an invented chromosome. Consequently no spatial track is generated.",
        "",
        "## Primary quartet",
        "",
        "Topology labels were fixed independently of strata: q1 is species clustering, q2 is same-mating-type clustering across species, and q3 is the third resolution. Frequencies are across genes within a stratum, not probabilities from one tree.",
        "",
        "| stratum | n | q1 | q2 | q3 | B | A (opposite - same mating type) | internal |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in support:
        stat = {row["statistic"]: row for row in bootstrap_rows if row["stratum"] == item["stratum"]}
        lines.append(f"| {item['stratum']} | {item['n_trees']} | {item['q1_species_frequency']:.3f} [{stat['q1_species']['ci_low']:.3f},{stat['q1_species']['ci_high']:.3f}] | {item['q2_mating_type_frequency']:.3f} [{stat['q2_mating_type']['ci_low']:.3f},{stat['q2_mating_type']['ci_high']:.3f}] | {item['q3_third_frequency']:.3f} [{stat['q3_third']['ci_low']:.3f},{stat['q3_third']['ci_high']:.3f}] | {stat['B']['estimate']:.5f} | {stat['A_opposite_minus_same']['estimate']:.5f} | {stat['internal_branch']['estimate']:.5f} |")
    lines += [
        "",
        f"The black stratum has q2 = {primary['q2_mating_type_frequency']:.3f} ({int(primary['q2_mating_type_count'])}/{primary['n_trees']}) and q1 = {primary['q1_species_frequency']:.3f}; the PAR control has q1 = {pars['q1_species_frequency']:.3f} ({pars['n_trees']}/{pars['n_trees']}) and q2 = 0; turquoise has q1 = {turquoise['q1_species_frequency']:.3f} and q2 = {turquoise['q2_mating_type_frequency']:.3f}. The black-to-turquoise contrast is therefore not a generic effect of every nonrecombining label: the younger turquoise stratum retains mostly species clustering.",
        "",
        "The corrected A statistic is the mean distance for opposite mating-type cross-species pairs minus the mean distance for same-mating-type cross-species pairs. B and all branch lengths are in the source tree's substitution-distance units, not calibrated coalescence times. The gene-level bootstrap intervals are descriptive; linked genes do not provide independent physical replicates.",
        "",
        "## All 15 quartets",
        "",
        f"All {len(set(item['quartet'] for item in all_summary))} possible four-taxon subsets are enumerated for the complete six-taxon trees. {n_meaningful} quartet identities (the three species-pair comparisons) have a biologically meaningful species-versus-mating-type labeling; the remaining 12 quartets retain fixed generic unrooted resolution labels and are not assigned an arrangement interpretation. The full assignment and summary tables preserve every quartet and every gene.",
        "",
        "The all-quartet heatmap uses the fixed second resolution for each quartet, not a post hoc choice of the largest local signal. The focal cat–sup quartet is the strongest and most directly interpretable contrast; the other species-pair quartets should be assessed as replication only after their species and mating-type metadata are independently confirmed.",
        "",
        "## Recombination regimes and alternatives",
        "",
        "The PARs are a categorical recombining control, not a quantitative crossover-rate track. The paper supports an older ancestral nonrecombining black region with trans-specific polymorphism and a younger turquoise species-specific stratum; it also reports additional colored strata, but this audit does not force an age ordering onto those labels without recovered physical boundaries. The observed black signal is compatible with long-lived mating-type-associated genealogies, but mutation-rate differences, alignment quality, gene-tree estimation error, residual gene conversion, selection, and introgression remain alternative explanations. The data do not make MSRC a unique causal explanation.",
        "",
        "## Spatial feasibility and blockers",
        "",
        "No chromosome-wide tracks, gene-order windows, or boundary tests are reported. A valid next step is to recover the authors' ortholog/contig placement table or an assembly-to-chromosome map, then repeat the analysis with block-aware windows and matched local controls. Until then, stratum-level comparisons are the defensible unit.",
        "",
        "## Files",
        "",
        "- `all_quartet_assignments.tsv` and `all_quartet_stratum_summary.tsv` — complete six-taxon quartet results.",
        "- `primary_stratum_bootstrap.tsv` — deterministic descriptive gene-level intervals.",
        "- `coordinate_validation_report.md` and `stratum_metadata.tsv` — mapping and biological regime documentation.",
        "- `microbotryum_stage1_topology_by_stratum.*`, `microbotryum_stage1_distance_by_stratum.*`, `microbotryum_stage1_all_quartet_heatmap.*`, and `microbotryum_stage1_representative_trees.*` — Stage 1 figures.",
        "",
        "References: [Lucotte et al. 2025](https://doi.org/10.1038/s41467-025-60222-5), [Dryad](https://doi.org/10.5061/dryad.0k6djhb89), [authors' analysis repository](https://github.com/eliselucotte/MicrobotryumOnDianthus).",
    ]
    (OUT / "stage1_analysis_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--raw-dir", type=Path, default=RAW); parser.add_argument("--bootstrap-replicates", type=int, default=1000); args = parser.parse_args()
    raw = args.raw_dir; rows = [row for row in parse_gene_rows(raw / "gene_trees") if row["root"] is not None]
    if not rows: raise SystemExit("No parseable local gene trees found")
    support_rows = list(csv.DictReader((STAGE0_OUT / "preliminary_quartet_support.tsv").open(), delimiter="\t"))
    support = [{key: (float(value) if key not in {"stratum"} and key != "n_trees" else (int(value) if key == "n_trees" else value)) for key, value in item.items()} for item in support_rows]
    for item in support:
        for key in ("q1_species_count", "q2_mating_type_count", "q3_third_count"):
            item[key] = int(item[key])
    assignments, summaries = all_assignments(rows); write_table(OUT / "all_quartet_assignments.tsv", assignments); write_table(OUT / "all_quartet_stratum_summary.tsv", summaries)
    bootstrap_rows = bootstrap(rows, args.bootstrap_replicates); write_table(OUT / "primary_stratum_bootstrap.tsv", bootstrap_rows)
    coord = gtf_coordinates(raw / "MPbarcode50_tsebra.gtf"); write_coordinate_report(rows, coord); write_stratum_metadata(support); make_figures(support, bootstrap_rows, summaries, rows); write_report(support, summaries, bootstrap_rows, rows)
    print(f"Stage 1: {len(rows)} trees; {len(summaries)} stratum-quartet summaries; {len(bootstrap_rows)} bootstrap rows; spatial tracks blocked by coordinate validation")


if __name__ == "__main__":
    main()
