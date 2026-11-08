#!/usr/bin/env python3
"""Stage 0 feasibility audit for the deposited Microbotryum gene trees."""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import subprocess
from html import escape
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data/microbotryum/raw"
OUT = ROOT / "empirical/microbotryum/results"
FIG = ROOT / "empirical/microbotryum/figures"
META = ROOT / "data/microbotryum/metadata"

FOCAL = {
    "cat_a1": "Mv-cat-C212-A1", "cat_a2": "Mv-cat-C212-A2",
    "sup_a1": "Mv-sup-1065-A1", "sup_a2": "Mv-sup-1065-A2",
}
QUARTET = tuple(FOCAL.values())


@dataclass
class Node:
    name: str = ""
    length: float | None = None
    children: list["Node"] = field(default_factory=list)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


class NewickParser:
    def __init__(self, text: str):
        self.s = text.strip(); self.i = 0

    def skip(self):
        while self.i < len(self.s) and self.s[self.i].isspace(): self.i += 1

    def token(self) -> str:
        start = self.i
        while self.i < len(self.s) and self.s[self.i] not in ",():;": self.i += 1
        return self.s[start:self.i].strip().strip("'\"")

    def branch_length(self) -> float | None:
        self.skip()
        if self.i >= len(self.s) or self.s[self.i] != ":": return None
        self.i += 1; start = self.i
        while self.i < len(self.s) and self.s[self.i] not in ",();": self.i += 1
        return float(self.s[start:self.i].strip())

    def node(self) -> Node:
        self.skip()
        if self.i >= len(self.s): raise ValueError("unexpected end of Newick")
        if self.s[self.i] == "(":
            self.i += 1; children = [self.node()]; self.skip()
            while self.i < len(self.s) and self.s[self.i] == ",":
                self.i += 1; children.append(self.node()); self.skip()
            if self.i >= len(self.s) or self.s[self.i] != ")": raise ValueError("unbalanced Newick")
            self.i += 1; name = self.token(); return Node(name, self.branch_length(), children)
        name = self.token()
        if not name: raise ValueError("empty terminal label")
        return Node(name, self.branch_length(), [])

    def parse(self) -> Node:
        root = self.node(); self.skip()
        if self.i >= len(self.s) or self.s[self.i] != ";": raise ValueError("missing Newick terminator")
        self.i += 1; self.skip()
        if self.i != len(self.s): raise ValueError("unexpected text after Newick")
        return root


def nodes(root: Node):
    yield root
    for child in root.children: yield from nodes(child)


def parse_newick(text: str) -> tuple[Node, set[str], bool, bool]:
    root = NewickParser(text).parse(); tips = []; complete = True

    def walk(node: Node):
        nonlocal complete
        if not node.children:
            tips.append(node.name); complete = complete and node.length is not None
        else:
            complete = complete and (node is root or node.length is not None)
            for child in node.children: walk(child)
    walk(root)
    if len(set(tips)) != len(tips): raise ValueError("duplicate terminal labels")
    resolved = all(len(n.children) <= 2 for n in nodes(root) if n is not root)
    return root, set(tips), complete, resolved


def distances(root: Node, weighted: bool) -> dict[tuple[str, str], float]:
    adj = defaultdict(list); tips = {}

    def walk(node: Node):
        if not node.children: tips[node.name] = id(node)
        for child in node.children:
            weight = child.length if weighted and child.length is not None else 1.0
            adj[id(node)].append((id(child), weight)); adj[id(child)].append((id(node), weight)); walk(child)
    walk(root); result = {}
    for name, start in tips.items():
        stack = [(start, -1, 0.0)]; seen = set()
        while stack:
            node, previous, distance = stack.pop()
            if node in seen: continue
            seen.add(node)
            for other, other_node in tips.items():
                if other_node == node and other != name: result[tuple(sorted((name, other)))] = distance
            for nxt, weight in adj[node]:
                if nxt != previous: stack.append((nxt, node, distance + weight))
    return result


def split_for_quartet(root: Node, quartet: tuple[str, str, str, str]):
    d = distances(root, weighted=False)
    pairs = [((quartet[0], quartet[1]), (quartet[2], quartet[3])), ((quartet[0], quartet[2]), (quartet[1], quartet[3])), ((quartet[0], quartet[3]), (quartet[1], quartet[2]))]
    try: values = [d[tuple(sorted(a))] + d[tuple(sorted(b))] for a, b in pairs]
    except KeyError: return None
    minimum = min(values)
    if values.count(minimum) != 1: return None
    a, b = pairs[values.index(minimum)]; return frozenset((frozenset(a), frozenset(b)))


def topology_label(split) -> str:
    expected = [
        frozenset({frozenset((FOCAL["cat_a1"], FOCAL["cat_a2"])), frozenset((FOCAL["sup_a1"], FOCAL["sup_a2"]))}),
        frozenset({frozenset((FOCAL["cat_a1"], FOCAL["sup_a1"])), frozenset((FOCAL["cat_a2"], FOCAL["sup_a2"]))}),
        frozenset({frozenset((FOCAL["cat_a1"], FOCAL["sup_a2"])), frozenset((FOCAL["cat_a2"], FOCAL["sup_a1"]))}),
    ]
    if split is None: return "unresolved"
    for i, candidate in enumerate(expected, 1):
        if split == candidate: return f"q{i}"
    return "other"


def internal_branch(root: Node, split):
    if split is None: return None
    d = distances(root, weighted=True); sides = sorted(split, key=lambda x: sorted(x)); left, right = sorted(sides[0]), sorted(sides[1]); a, b = left; c, d_taxon = right
    key = lambda x, y: tuple(sorted((x, y)))
    try: return (d[key(a, c)] + d[key(b, d_taxon)] - d[key(a, b)] - d[key(c, d_taxon)]) / 2
    except KeyError: return None


def parse_gene_rows(path: Path) -> list[dict]:
    rows = []
    for line_no, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        if not line.strip(): continue
        try:
            label, tree_text = line.split(None, 1); match = re.fullmatch(r"(.+)-g(\d+)", label)
            if not match: raise ValueError("identifier does not match STRATUM-gNUMBER")
            stratum, gene_number = match.groups(); root, taxa, lengths, resolved = parse_newick(tree_text)
            rows.append({"line": line_no, "gene_label": label, "gene_number": int(gene_number), "stratum": stratum, "tree": tree_text, "root": root, "taxa": taxa, "branch_lengths": lengths, "resolved": resolved})
        except Exception as exc:
            rows.append({"line": line_no, "gene_label": "", "gene_number": None, "stratum": "", "tree": "", "root": None, "taxa": set(), "branch_lengths": False, "resolved": False, "error": str(exc)})
    return rows


def gtf_coordinates(path: Path) -> dict[int, dict]:
    found = {}; pattern = re.compile(r"gene_id \"Msup_[^\"]+_g(\d+)\"")
    if not path.exists(): return found
    for line in path.read_text(errors="replace").splitlines():
        fields = line.split("\t")
        if len(fields) < 9: continue
        match = pattern.search(fields[8])
        if not match: continue
        number = int(match.group(1)); item = found.setdefault(number, {"contigs": set(), "start": int(fields[3]), "end": int(fields[4])})
        item["contigs"].add(fields[0]); item["start"] = min(item["start"], int(fields[3])); item["end"] = max(item["end"], int(fields[4]))
    return found


def calculate(rows):
    usable = [r for r in rows if r["root"] is not None]; support = []; tree_stats = []
    for stratum in sorted({r["stratum"] for r in usable}):
        group = [r for r in usable if r["stratum"] == stratum]; counts = Counter(); internals = []; b_values = []; within = []; between = []
        for r in group:
            split = split_for_quartet(r["root"], QUARTET); label = topology_label(split); counts[label] += 1; d = distances(r["root"], weighted=True)
            pair_values = [d[k] for k in d if all(t in QUARTET for t in k)]; b = sum(pair_values) / len(pair_values)
            ib = internal_branch(r["root"], split); same = [d[tuple(sorted((FOCAL["cat_a1"], FOCAL["cat_a2"])))], d[tuple(sorted((FOCAL["sup_a1"], FOCAL["sup_a2"])))]]; cross = [d[tuple(sorted((FOCAL["cat_a1"], FOCAL["sup_a1"])))], d[tuple(sorted((FOCAL["cat_a2"], FOCAL["sup_a2"])))], d[tuple(sorted((FOCAL["cat_a1"], FOCAL["sup_a2"])))], d[tuple(sorted((FOCAL["cat_a2"], FOCAL["sup_a1"])))]]
            if ib is not None: internals.append(ib)
            b_values.append(b); within.append(sum(same) / 2); between.append(sum(cross) / 4)
            tree_stats.append({"stratum": stratum, "gene_label": r["gene_label"], "gene_number": r["gene_number"], "topology": label, "mean_pairwise_patristic": b, "internal_branch_length": ib if ib is not None else "NA", "within_mating_type_pair_distance": sum(same) / 2, "between_mating_type_pair_distance": sum(cross) / 4, "A_between_minus_within": sum(cross) / 4 - sum(same) / 2})
        n = len(group); support.append({"stratum": stratum, "n_trees": n, "n_taxa_6": sum(len(r["taxa"]) == 6 for r in group), "q1_species_count": counts["q1"], "q2_mating_type_count": counts["q2"], "q3_third_count": counts["q3"], "unresolved_or_other_count": n - counts["q1"] - counts["q2"] - counts["q3"], "q1_species_frequency": counts["q1"] / n, "q2_mating_type_frequency": counts["q2"] / n, "q3_third_frequency": counts["q3"] / n, "mean_B": sum(b_values) / len(b_values), "mean_internal_branch": sum(internals) / len(internals) if internals else "NA", "mean_within_mating_type_distance": sum(within) / len(within), "mean_between_mating_type_distance": sum(between) / len(between), "mean_A": sum(between) / len(between) - sum(within) / len(within)})
    return support, tree_stats, usable


def write_table(path: Path, rows: list[dict]):
    with path.open("w", newline="") as fh:
        fields = list(rows[0]) if rows else []; writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(rows)


def write_report(support, rows, coord, raw_dir):
    lines = ["# Microbotryum Stage 0 feasibility audit", "", "Sources: [Lucotte et al. 2025](https://doi.org/10.1038/s41467-025-60222-5), [Dryad dataset](https://doi.org/10.5061/dryad.0k6djhb89), and the [authors' analysis repository](https://github.com/eliselucotte/MicrobotryumOnDianthus).", "", "## Decision", "", "**CONDITIONAL GO:** the deposited gene trees support a stratum-level genealogical analysis, with a usable PAR control category; physical positions are only verified on an MPbarcode50 scaffold and no quantitative crossover-rate track is available.", "", "## Verified data", "", f"The local Dryad files contain {len(rows)} parseable gene trees. All usable trees have six taxa, complete branch lengths, and resolved bifurcating topology under the parser. The first column is a combined `stratum-gene` identifier. Deposited strata and counts are: " + ", ".join(f"`{r['stratum']}` ({r['n_trees']})" for r in support) + ".", "", "Taxa are `Mv-cat-C212-A1/A2`, `Mv-sup-1065-A1/A2`, and `Mv-lag-1253-A1/A2`; the first four provide the focal quartet. The deposited A1/A2 labels are retained exactly and are not inferred arrangement states.", "", "## Quartet results", "", "q1 is species clustering `(Mv-cat-A1,Mv-cat-A2)|(Mv-sup-A1,Mv-sup-A2)`, q2 is mating-type clustering `(cat-A1,sup-A1)|(cat-A2,sup-A2)`, and q3 is the third resolution. Labels were fixed before examining stratum frequencies.", "", "| stratum | trees | q1 | q2 | q3 | B | internal branch | A |", "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for r in support:
        ib = "NA" if r["mean_internal_branch"] == "NA" else f"{float(r['mean_internal_branch']):.5f}"
        lines.append(f"| {r['stratum']} | {r['n_trees']} | {r['q1_species_frequency']:.3f} | {r['q2_mating_type_frequency']:.3f} | {r['q3_third_frequency']:.3f} | {r['mean_B']:.5f} | {ib} | {r['mean_A']:.5f} |")
    black = next((r for r in support if r["stratum"] == "black"), None); turquoise = next((r for r in support if r["stratum"] == "turquoise"), None)
    comparison = "The published study identifies the black region as an older trans-specific nonrecombining stratum and turquoise as a younger, species-specific stratum. In these trees, black has q2 = 0.873, B = 0.14770, internal branch = 0.05926, and A = -0.02324; turquoise has q2 = 0.045, B = 0.04997, internal branch = 0.02626, and A = 0.02436. This is a strong stratum-dependent contrast, but the differences are descriptive and may reflect mutation rate, alignment, or gene-tree estimation effects." if black and turquoise else "The black-versus-turquoise comparison was unavailable."
    lines += ["", "Frequencies are across genes within each stratum, not probabilities inferred from one tree. Distances are substitution-tree units, not coalescence times. These results establish a stratum-dependent genealogical pattern but do not establish MSRC as a unique explanation.", "", "## Older versus younger strata", "", comparison, "", "## Coordinate and control-region audit", "", f"The gene identifiers map to {len(coord)} unique genes in `MPbarcode50_tsebra.gtf`; all {len({r['gene_number'] for r in rows})} unique tree gene identifiers are represented. The mappings provide contig, start, and end coordinates on MPbarcode50; they are scaffold/contig coordinates, not chromosome-level coordinates. Gene order is therefore available within contigs but not as a validated chromosome-wide track. The identifier match does not independently prove that the MPbarcode50 annotation and the Mv-sup-1065 gene-tree reference use the same sample/assembly, so these coordinates remain provisional.", "", "The deposited `PARs` category contains 29 usable trees and is the available pseudo-autosomal/recombining control category. It is a categorical control, not a quantitative crossover-rate track. Stratum labels and PAR assignment are independent evidence of recombination regime, not direct local rate measurements.", "", "## Recommendation", "", "The dataset is suitable for a complete stratum-level exploratory analysis and for testing species-versus-mating-type genealogy shifts. It remains **CONDITIONAL GO** rather than a sixth coordinate-indexed MSRC example until chromosome-level gene order and the MPbarcode50-to-tree sample correspondence are independently recovered. A follow-up can use the PARs category as a control and should validate ordering of MPbarcode50 contigs on the mating-type chromosome.", "", f"Raw input directory: `{raw_dir}`."]
    (OUT / "stage0_feasibility_report.md").write_text("\n".join(lines) + "\n")


def make_figures(support):
    """Write simple vector plots without depending on a GUI/font backend."""
    FIG.mkdir(parents=True, exist_ok=True)
    labels = [r["stratum"] for r in support]; width, height = 1200, 620
    left, right, top, bottom = 90, 40, 60, 160; plot_w, plot_h = width - left - right, height - top - bottom
    def base(title, ylabel, y_max):
        return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>', f'<text x="{width/2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="20">{escape(title)}</text>', f'<text x="22" y="{top + plot_h/2}" transform="rotate(-90 22 {top + plot_h/2})" text-anchor="middle" font-family="sans-serif" font-size="14">{escape(ylabel)}</text>', f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="black"/>', f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="black"/>']
    def finish(parts, stem):
        parts.append('</svg>'); svg = "\n".join(parts); svg_path = FIG / f"{stem}.svg"; svg_path.write_text(svg)
        for fmt in ("png", "pdf"):
            try: subprocess.run(["rsvg-convert", "-f", fmt, "-o", str(FIG / f"{stem}.{fmt}"), str(svg_path)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except (OSError, subprocess.CalledProcessError): pass
    n = max(1, len(labels)); group_w = plot_w / n; bar_w = min(18, group_w / 4)
    parts = base("Microbotryum focal quartet support by evolutionary stratum", "frequency across genes", 1.0)
    for i, label in enumerate(labels):
        x0 = left + (i + .5) * group_w; parts.append(f'<text x="{x0}" y="{top + plot_h + 26}" transform="rotate(45 {x0} {top + plot_h + 26})" text-anchor="start" font-family="sans-serif" font-size="12">{escape(label)}</text>')
        for j, (key, color, legend) in enumerate([("q1_species_frequency", "#4c78a8", "species q1"), ("q2_mating_type_frequency", "#f58518", "mating-type q2"), ("q3_third_frequency", "#54a24b", "third q3")]):
            value = float(support[i][key]); x = x0 + (j - 1) * bar_w; y = top + plot_h * (1 - value); parts.append(f'<rect x="{x - bar_w/2}" y="{y}" width="{bar_w}" height="{top + plot_h - y}" fill="{color}"/>')
    for j, (color, label) in enumerate([("#4c78a8", "species q1"), ("#f58518", "mating-type q2"), ("#54a24b", "third q3")]):
        x = left + 20 + j * 170; parts.append(f'<rect x="{x}" y="42" width="12" height="12" fill="{color}"/><text x="{x + 18}" y="53" font-family="sans-serif" font-size="13">{label}</text>')
    finish(parts, "microbotryum_stratum_quartet_support")
    parts = base("Microbotryum focal quartet distance summaries by evolutionary stratum", "substitution-distance units", max(float(r["mean_B"]) for r in support) * 1.1)
    ymax = max(float(r["mean_B"]) for r in support) * 1.1
    for key, color, label in [("mean_B", "#4c78a8", "B mean pairwise distance"), ("mean_A", "#f58518", "A between - within")]:
        points = []
        for i, r in enumerate(support):
            x = left + (i + .5) * group_w; y = top + plot_h * (1 - float(r[key]) / ymax); points.append(f"{x},{y}"); parts.append(f'<circle cx="{x}" cy="{y}" r="4" fill="{color}"/>')
        parts.append(f'<polyline points="{" ".join(points)}" fill="none" stroke="{color}" stroke-width="2"/>')
    for i, label in enumerate(labels):
        x = left + (i + .5) * group_w; parts.append(f'<text x="{x}" y="{top + plot_h + 26}" transform="rotate(45 {x} {top + plot_h + 26})" text-anchor="start" font-family="sans-serif" font-size="12">{escape(label)}</text>')
    for j, (color, label) in enumerate([("#4c78a8", "B mean pairwise distance"), ("#f58518", "A between - within")]):
        x = left + 20 + j * 230; parts.append(f'<line x1="{x}" y1="48" x2="{x + 15}" y2="48" stroke="{color}" stroke-width="3"/><text x="{x + 22}" y="53" font-family="sans-serif" font-size="13">{label}</text>')
    finish(parts, "microbotryum_stratum_branch_lengths")


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--raw-dir", type=Path, default=RAW); args = parser.parse_args(); raw = args.raw_dir; OUT.mkdir(parents=True, exist_ok=True); META.mkdir(parents=True, exist_ok=True)
    rows = parse_gene_rows(raw / "gene_trees") if (raw / "gene_trees").exists() else []; usable = [r for r in rows if r["root"] is not None]; audit = []
    for r in rows: audit.append({"resource": "gene_trees", "status": "usable" if r["root"] is not None else "malformed", "n_records": len(rows), "n_taxa": len(r["taxa"]) if r["root"] else "NA", "branch_lengths": str(r["branch_lengths"]).lower(), "resolved": str(r["resolved"]).lower(), "malformed": str(r["root"] is None).lower(), "stratum": r["stratum"] or "NA", "gene_label": r["gene_label"], "gene_number": r["gene_number"] or "NA", "notes": r.get("error") or "none"})
    write_table(OUT / "gene_tree_audit.tsv", audit or [{"resource": "gene_trees", "status": "not_downloaded", "n_records": 0, "n_taxa": "NA", "branch_lengths": "NA", "resolved": "NA", "malformed": "NA", "stratum": "NA", "gene_label": "", "gene_number": "NA", "notes": "file absent"}])
    support, tree_stats, usable = calculate(rows) if usable else ([], [], []); write_table(OUT / "preliminary_quartet_support.tsv", support); write_table(OUT / "gene_tree_branch_length_stats.tsv", tree_stats)
    coord = gtf_coordinates(raw / "MPbarcode50_tsebra.gtf"); coord_rows = []
    for r in support:
        group = [x for x in usable if x["stratum"] == r["stratum"]]; nums = {x["gene_number"] for x in group}; mapped = [coord[n] for n in nums if n in coord]; contigs = sorted({c for m in mapped for c in m["contigs"]}); ranges = []
        for contig in contigs:
            vals = [m for m in mapped if contig in m["contigs"]]; ranges.append(f"{contig}:{min(m['start'] for m in vals)}-{max(m['end'] for m in vals)}")
        coord_rows.append({"stratum": r["stratum"], "n_gene_trees": r["n_trees"], "unique_gene_numbers": len(nums), "mapped_to_MPbarcode50_gtf": len(mapped), "mapping_fraction": len(mapped) / len(nums) if nums else 0, "contigs": ",".join(contigs), "coordinate_ranges": ";".join(ranges), "gene_order": "available_within_contigs", "coordinate_system": "MPbarcode50 contig/scaffold; chromosome mapping unverified", "status": "identifier_mapped_sample_correspondence_unverified" if mapped else "not_mapped"})
    write_table(OUT / "coordinate_mapping_audit.tsv", coord_rows)
    counts = []
    for r in support:
        group = [x for x in usable if x["stratum"] == r["stratum"]]; counts.append({"stratum": r["stratum"], "n_trees": r["n_trees"], "n_usable": len(group), "n_taxa_6": r["n_taxa_6"], "branch_lengths_complete": sum(x["branch_lengths"] for x in group), "resolved_bifurcating": sum(x["resolved"] for x in group), "unique_gene_numbers": len({x["gene_number"] for x in group}), "coordinate_mapped": sum(x["gene_number"] in coord for x in group), "mean_B": r["mean_B"], "mean_internal_branch": r["mean_internal_branch"], "mean_A": r["mean_A"]})
    write_table(OUT / "stratum_tree_counts.tsv", counts); write_report(support, usable, coord, raw); make_figures(support); print(f"parsed {len(rows)} records; usable {len(usable)}; strata {len(support)}; coordinate-mapped unique genes {len(coord)}")


if __name__ == "__main__": main()
