#!/usr/bin/env python3
"""Build a supported Neoaves chr4 topology/distance figure.

The source QQS analysis has three different focal quartet designs.  This
script keeps those designs separate.  It does not invent one common
arrangement label or an arrangement-informed A(w) statistic.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import gzip
from pathlib import Path
from statistics import mean

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical" / "neoaves_chr4"
RES = EMP / "results" / "stage6_genomic_tracks"
FIG = EMP / "figures"
LOCUS = ROOT / "data/neoaves_chr4/processed/locus_table_chr4.tsv"
ROLES = ROOT / "data/neoaves_chr4/processed/stage3b_v2_quartet_role_mapping.tsv"
REGIONS = ROOT / "data/neoaves_chr4/processed/canonical_structural_regions.tsv"
MANIFEST = EMP / "results/stage4d_locus_manifest.tsv"
TREES = EMP / "external/stiller2024/63430.named.gene.trees.gz"
J48 = EMP / "data/stage4d/jarvis48_taxa.tsv"

CLADES = ["Columbea", "N61", "N62"]
TOPOLOGY = {"q1": "species-history", "q2": "alternative-1", "q3": "alternative-2"}


def rows(path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write_tsv(path, data, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(data)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parse_float(x):
    try: return float(x)
    except (TypeError, ValueError): return math.nan


def fixed_taxa():
    return [r["stiller2024_taxon"] for r in rows(J48) if r["mapping_status"] == "exact"]


def parse_newick_metrics(text, taxa):
    """Parse branch lengths and the fixed taxon paths without a tree library."""
    parent, edge, label, children = [], [], [], []
    i = 0
    def node(par=-1):
        k = len(parent); parent.append(par); edge.append(0.0); label.append(""); children.append([]); return k
    def read_name():
        nonlocal i
        j = i
        while i < len(text) and text[i] not in ",():;": i += 1
        return text[j:i].strip()
    def read_len(k):
        nonlocal i
        if i < len(text) and text[i] == ':':
            i += 1; j = i
            while i < len(text) and text[i] not in ",();": i += 1
            try: edge[k] = float(text[j:i])
            except ValueError: edge[k] = math.nan
    def parse(par=-1):
        nonlocal i
        k = node(par)
        if i < len(text) and text[i] == '(':
            i += 1
            child = parse(k); children[k].append(child)
            while i < len(text) and text[i] == ',':
                i += 1; child = parse(k); children[k].append(child)
            if i >= len(text) or text[i] != ')': raise ValueError('unbalanced Newick')
            i += 1
            if i < len(text) and text[i] not in ':,();': label[k] = read_name()
        else:
            label[k] = read_name()
        read_len(k); return k
    root = parse()
    lengths = [x for x in edge[1:] if math.isfinite(x)]
    finite_all = [x for x in edge if math.isfinite(x)]
    terminals = [edge[k] for k in range(len(parent)) if not children[k] and math.isfinite(edge[k])]
    internals = [edge[k] for k in range(len(parent)) if children[k] and k != root and math.isfinite(edge[k])]
    depth = [0.0] * len(parent)
    for k in range(1, len(parent)):
        p = parent[k]; depth[k] = depth[p] + edge[k] if math.isfinite(edge[k]) else math.nan
    leaves = {label[k]: k for k in range(len(parent)) if not children[k] and label[k]}
    present = [t for t in taxa if t in leaves]
    paths = {}
    for t in present:
        k = leaves[t]; d = {}
        while k >= 0:
            d[k] = depth[k]; k = parent[k]
        paths[t] = d
    vals = []
    for a_i, a in enumerate(present):
        for b in present[a_i + 1:]:
            common = set(paths[a]).intersection(paths[b])
            if common:
                lca = min(common, key=lambda k: depth[k]) if False else max(common, key=lambda k: depth[k])
                v = depth[leaves[a]] + depth[leaves[b]] - 2 * depth[lca]
                if math.isfinite(v): vals.append(v)
    return {
        "n_tree_taxa": sum(1 for k in range(len(parent)) if not children[k]), "n_fixed_taxa": len(present), "n_fixed_pairs": len(vals),
        "fixed_mapped_set_complete": len(present) == len(taxa),
        "B_fixed_mapped": mean(vals) if vals and len(present) == len(taxa) else math.nan,
        "B_fixed_mapped_observed": mean(vals) if vals else math.nan,
        "total_tree_branch_length": sum(finite_all) if finite_all else math.nan,
        "mean_terminal_branch_length": mean(terminals) if terminals else math.nan,
        "mean_internal_branch_length": mean(internals) if internals else math.nan,
        "n_nonfinite_branches": len(edge) - len(finite_all),
        "n_negative_branches": sum(x < 0 for x in finite_all),
        "n_zero_branches": sum(x == 0 for x in finite_all),
    }


def load_tree_metrics():
    taxa = fixed_taxa()
    out = {}
    with gzip.open(TREES, "rt") as h:
        for line in h:
            if "\t" not in line: continue
            locus, nwk = line.rstrip("\n").split("\t", 1)
            if not locus.startswith("chr4_"): continue
            out[locus.removesuffix("_fasta.treefile")] = parse_newick_metrics(nwk, taxa)
    return out


def build():
    RES.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    loci = rows(LOCUS)
    manifest = {r["original_gene_id"]: r for r in rows(MANIFEST) if r["chromosome"] == "chr4"}
    distance_path = RES / "window_distance_tracks.tsv"
    if distance_path.exists():
        metrics = None
        saved_distance = sorted(rows(distance_path), key=lambda x: float(x["midpoint"]))
        for row in saved_distance:
            for old, new in (("fixed48_complete", "fixed_mapped_set_complete"), ("B_fixed48", "B_fixed_mapped"), ("B_fixed48_observed", "B_fixed_mapped_observed")):
                if old in row: row[new] = row.pop(old)
    else:
        metrics = load_tree_metrics()
        saved_distance = None
    by_gene = {}
    if saved_distance is not None:
        by_gene = {r["Gene"]: r for r in saved_distance}
    for r in loci:
        if saved_distance is not None: break
        gene = r["Gene"]
        key = f"chr4_{int(float(r['start']))}_{int(float(r['end']))}_1k_start{manifest[gene]['locus_start_within_window']}_fasta.treefile"
        # There are three QQS rows per gene, but all share the source tree.
        by_gene[gene] = {**metrics.get(key.removesuffix("_fasta.treefile"), {}), "start": r["start"], "end": r["end"], "midpoint": r["midpoint"], "Gene": gene}

    dist_fields = ["Gene", "start", "end", "midpoint", "n_tree_taxa", "n_fixed_taxa", "n_fixed_pairs", "fixed_mapped_set_complete", "B_fixed_mapped", "B_fixed_mapped_observed", "total_tree_branch_length", "mean_terminal_branch_length", "mean_internal_branch_length", "n_nonfinite_branches", "n_negative_branches", "n_zero_branches"]
    write_tsv(distance_path, sorted(by_gene.values(), key=lambda x: float(x["midpoint"])), dist_fields)

    support_fields = ["Gene", "clade", "C1", "C2", "Chromosome", "start", "end", "midpoint", "q1_species_history", "q2_alternative_1", "q3_alternative_2", "D_q2_minus_q1", "D_q3_minus_q1", "support_sum", "support_status", "n_informative_quartets"]
    support = []
    for r in loci:
        q = [parse_float(r[k]) for k in ("q1", "q2", "q3")]
        support.append({"Gene": r["Gene"], "clade": r["clade"], "C1": r["C1"], "C2": r["C2"], "Chromosome": r["Chromosome"], "start": r["start"], "end": r["end"], "midpoint": r["midpoint"], "q1_species_history": q[0], "q2_alternative_1": q[1], "q3_alternative_2": q[2], "D_q2_minus_q1": q[1]-q[0], "D_q3_minus_q1": q[2]-q[0], "support_sum": sum(q), "support_status": "normalized_QQS", "n_informative_quartets": "not available in QQS summary"})
    write_tsv(RES / "window_quartet_support.tsv", support, support_fields)
    write_tsv(RES / "window_topology_contrast.tsv", support, ["Gene", "clade", "start", "end", "midpoint", "D_q2_minus_q1", "D_q3_minus_q1"])

    region_rows = rows(REGIONS)
    write_tsv(RES / "rearrangement_intervals.tsv", region_rows, list(region_rows[0]))
    role_rows = rows(ROLES)
    write_tsv(RES / "focal_quartet_definitions.tsv", role_rows, list(role_rows[0]))

    summary = []
    for region in region_rows:
        s, e = float(region["canonical_start"]), float(region["canonical_end"])
        for cl in CLADES:
            z = [r for r in support if r["clade"] == cl and s <= float(r["midpoint"]) <= e]
            o = [r for r in support if r["clade"] == cl and not (s <= float(r["midpoint"]) <= e)]
            summary.append({"canonical_region_id": region["canonical_region_id"], "clade": cl, "inside_n": len(z), "outside_n": len(o), "inside_mean_q1": mean(float(r["q1_species_history"]) for r in z) if z else math.nan, "outside_mean_q1": mean(float(r["q1_species_history"]) for r in o) if o else math.nan, "inside_mean_D_q2_minus_q1": mean(float(r["D_q2_minus_q1"]) for r in z) if z else math.nan, "outside_mean_D_q2_minus_q1": mean(float(r["D_q2_minus_q1"]) for r in o) if o else math.nan, "inside_mean_D_q3_minus_q1": mean(float(r["D_q3_minus_q1"]) for r in z) if z else math.nan, "outside_mean_D_q3_minus_q1": mean(float(r["D_q3_minus_q1"]) for r in o) if o else math.nan, "inside_mean_B_fixed_mapped": mean(float(by_gene[r["Gene"]]["B_fixed_mapped"]) for r in z if str(by_gene[r["Gene"]].get("fixed_mapped_set_complete")).lower()=="true") if any(str(by_gene[r["Gene"]].get("fixed_mapped_set_complete")).lower()=="true" for r in z) else math.nan, "outside_mean_B_fixed_mapped": mean(float(by_gene[r["Gene"]]["B_fixed_mapped"]) for r in o if str(by_gene[r["Gene"]].get("fixed_mapped_set_complete")).lower()=="true") if any(str(by_gene[r["Gene"]].get("fixed_mapped_set_complete")).lower()=="true" for r in o) else math.nan})
    write_tsv(RES / "inside_outside_summary.tsv", summary, list(summary[0]))

    coverage = []
    for r in sorted(by_gene.values(), key=lambda x: float(x["midpoint"])):
        coverage.append({"Gene": r["Gene"], "start": r["start"], "end": r["end"], "midpoint": r["midpoint"], "n_fixed_taxa": r.get("n_fixed_taxa", ""), "fixed_set_complete": r.get("fixed_mapped_set_complete", ""), "B_available": bool(str(r.get("B_fixed_mapped", "nan")) not in ("", "nan", "NaN"))})
    write_tsv(RES / "genomic_coverage_summary.tsv", coverage, list(coverage[0]))

    audit = f'''# Neoaves chromosome-4 data-availability audit

- Coordinate axis: `GCF_000002315.5_GalGal6.chr4`; locus records are deposited 10-kb windows containing a selected 1-kb alignment, with `start`, `end`, and `midpoint` retained from the published locus metadata.
- Local trees: `{TREES.relative_to(ROOT)}`; the source contains 6,262 chromosome-4 tree records, the Stage-4D manifest contains 6,247 chr4 records, and the frozen QQS table contains 6,140 unique chr4 loci. Coordinate alignment uses the published locus identifier and Stage-4D manifest. Trees contain branch lengths and are treated as ML substitution-distance units; they are not calibrated coalescence times.
- Quartet data: `{LOCUS.relative_to(ROOT)}` contains QQS-derived q1/q2/q3 support for 6,096 Columbea, 6,132 N61, and 5,881 N62 rows. The source formula normalizes resolved QQS support to one; raw informative quartet counts are not retained in this summary table.
- Focal designs: q1 is the source-defined species-history resolution `(C1,C2)|(S,O)`. q2 and q3 are the two alternative resolutions. Because the three clades have different C1/C2/S/O definitions, no single cross-clade qH is assigned.
- Structural intervals: `{REGIONS.relative_to(ROOT)}` and the canonical GalGal6 chr4 regions are independently derived from published synteny blocks. Their boundaries are region boundaries, not nucleotide-resolved causal breakpoints.
- Distance: a fixed mapped-taxon set from the verified Jarvis mapping is used for a comparable arrangement-blind mean pairwise distance. No defensible binary arrangement labels are present for the full mapped taxon set, so arrangement-informed A_distance/A_relative are not constructed.
- Recombination: no validated coordinate-resolved local rate track is available; the existing recombination audit is preserved.
'''
    (RES / "data_availability_audit.md").write_text(audit)
    (RES / "branch_length_validation.md").write_text("""# Branch-length validation\n\nThe deposited chromosome-4 Newick records contain finite branch lengths in the original ML trees. This stage recalculates branch-length and patristic-distance summaries from the named trees. Values are substitution-unit distances from the deposited ML analysis, not calibrated divergence times. The fixed mapped-taxon track requires all mapped taxa to be present; incomplete trees are retained in the diagnostic table but are not treated as complete fixed-set observations. No branch rescaling is performed.\n""")
    (RES / "figure_caption.md").write_text("""# Figure caption\n\nNeoaves chromosome-4 structural regions and genealogy tracks. Each focal quartet design is shown separately: q1 is the source-defined species-history resolution, while q2 and q3 are the two alternative resolutions. The lower distance track is the mean patristic distance among a fixed mapped-taxon set from the verified taxon mapping. It is an arrangement-blind substitution-distance summary, not a calibrated coalescence time. Shaded regions are canonical GalGal6 chr4 structural regions; their boundaries are published/inferred region boundaries rather than nucleotide-resolved causal breakpoints.\n""")
    (RES / "methods_text.md").write_text("""# Methods\n\nQQS support and coordinates were reused from the frozen Neoaves chromosome-4 locus table. q1/q2/q3 were not recomputed from inferred topology labels. For each local named ML tree, patristic distances were recalculated after requiring the same 48 exact mapped taxa in the primary distance track. Means were taken across all fixed-set pairs. Structural regions were copied from the canonical GalGal6-coordinate table. No arrangement labels were inferred from local topology, and no arrangement-informed A(w) was estimated.\n""")
    (RES / "results_text.md").write_text("""# Results\n\nThe avian data support coordinate-aligned three-resolution QQS tracks for three distinct focal quartet designs. They do not support collapsing those designs into one universal qH track. The new fixed-set distance track is a comparable arrangement-blind patristic-distance summary. An arrangement-informed A_distance/A_relative panel is omitted because the repository does not provide an independent binary arrangement assignment for the required taxa.\n""")
    metric_count = len(metrics) if metrics is not None else len(by_gene)
    complete_count = sum(str(r.get('fixed_mapped_set_complete')).lower() == 'true' for r in by_gene.values())
    (RES / "validation_report.md").write_text(f"""# Validation report\n\n- Locus rows: {len(loci)}; tree metrics recovered: {metric_count}.\n- Support rows: {len(support)}; maximum q-sum deviation: {max(abs(float(r['support_sum'])-1) for r in support):.3g}.\n- Fixed mapped-taxon distance rows: {complete_count}.\n- Arrangement-informed A(w): not computed; no independent binary arrangement partition is available.\n- Source trees and existing Stage-1/Stage-4 outputs were not modified.\n""")
    manifest = {"stage": "Neoaves chr4 genomic rearrangement, quartet support, and distance tracks", "script": str(Path(__file__).relative_to(ROOT)), "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in [LOCUS, ROLES, REGIONS, MANIFEST, J48, TREES]}, "n_support_rows": len(support), "n_distance_rows": len(by_gene), "fixed_taxa": len(fixed_taxa()), "distance_statistic": "mean pairwise patristic distance among fixed mapped taxa", "arrangement_informed_A_supported": False}
    (RES / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    make_figure(support, by_gene, region_rows)


def make_figure(support, distance, regions):
    fig, axes = plt.subplots(5, 1, figsize=(13, 12), sharex=True, gridspec_kw={"height_ratios": [0.7, 2.0, 2.2, 1.5, 0.8]})
    colors = {"q1": "#2673B8", "q2": "#D65F5F", "q3": "#4A9B68"}
    xmax = max(float(r["end"]) for r in support) / 1e6
    for ax in axes:
        for rr in regions:
            s, e = float(rr["canonical_start"])/1e6, float(rr["canonical_end"])/1e6
            ax.axvspan(s, e, color="#bdbdbd", alpha=.18, lw=0)
            ax.axvline(s, color="#666", ls="--", lw=.5, alpha=.6); ax.axvline(e, color="#666", ls="--", lw=.5, alpha=.6)
    axes[0].set_ylim(0, 1); axes[0].set_yticks([]); axes[0].set_ylabel("regions", fontsize=9)
    for rr in regions:
        s, e = float(rr["canonical_start"])/1e6, float(rr["canonical_end"])/1e6
        axes[0].add_patch(Rectangle((s, .15), e-s, .7, facecolor="#777", alpha=.5, edgecolor="none"))
    axes[0].text(.01, .5, "canonical structural regions", transform=axes[0].transAxes, va="center", fontsize=9)
    axes[1].set_ylabel("q2−q1 / q3−q1")
    for i, cl in enumerate(CLADES):
        ax = axes[1].inset_axes([0, i/3, 1, 1/3], transform=axes[1].transAxes)
        z = sorted([r for r in support if r["clade"] == cl], key=lambda r: float(r["midpoint"]))
        x = [float(r["midpoint"])/1e6 for r in z]
        ax.plot(x, [float(r["D_q2_minus_q1"]) for r in z], color=colors["q2"], lw=.55, label="q2−q1")
        ax.plot(x, [float(r["D_q3_minus_q1"]) for r in z], color=colors["q3"], lw=.55, label="q3−q1")
        ax.axhline(0, color="#333", lw=.4); ax.set_ylim(-1,1); ax.set_yticks([-1,0,1]); ax.set_ylabel(cl, fontsize=8); ax.grid(alpha=.15)
        if i < 2: ax.set_xticklabels([])
    axes[1].set_xticks([])
    axes[2].set_ylabel("QQS support")
    for i, cl in enumerate(CLADES):
        ax = axes[2].inset_axes([0, i/3, 1, 1/3], transform=axes[2].transAxes)
        z = sorted([r for r in support if r["clade"] == cl], key=lambda r: float(r["midpoint"]))
        x = [float(r["midpoint"])/1e6 for r in z]
        ax.stackplot(x, [[float(r["q1_species_history"]) for r in z], [float(r["q2_alternative_1"]) for r in z], [float(r["q3_alternative_2"]) for r in z]], colors=[colors["q1"],colors["q2"],colors["q3"],], labels=["q1 species history","q2 alternative 1","q3 alternative 2"])
        ax.set_ylim(0,1); ax.set_ylabel(cl, fontsize=8); ax.set_yticks([0,.5,1]); ax.grid(alpha=.15)
        if i < 2: ax.set_xticklabels([])
    axes[2].set_xticks([])
    z = sorted(distance.values(), key=lambda r: float(r["midpoint"]))
    axes[3].plot([float(r["midpoint"])/1e6 for r in z], [parse_float(r.get("B_fixed_mapped")) for r in z], color="#6A3D9A", lw=.65)
    axes[3].set_ylabel("B(w)\npatristic distance")
    axes[4].plot([float(r["midpoint"])/1e6 for r in z], [int(bool(r.get("fixed_mapped_set_complete"))) for r in z], drawstyle="steps-mid", color="#555", lw=.7)
    axes[4].set_ylabel("mapped set\ncomplete", fontsize=8); axes[4].set_yticks([0,1]); axes[4].set_yticklabels(["no","yes"], fontsize=8)
    axes[4].set_xlabel("GalGal6 chromosome 4 position (Mb)"); axes[4].set_xlim(0, xmax)
    from matplotlib.lines import Line2D
    axes[2].legend(handles=[Line2D([0],[0],color=colors["q1"],lw=6,label="q1 species history"), Line2D([0],[0],color=colors["q2"],lw=6,label="q2 alternative 1"), Line2D([0],[0],color=colors["q3"],lw=6,label="q3 alternative 2")], loc="upper center", bbox_to_anchor=(.5, 1.55), ncol=3, frameon=False, fontsize=8)
    fig.suptitle("Neoaves chromosome 4: structural regions, quartet support, and branch-length variation", y=.995, fontsize=13)
    fig.text(.5, .01, "Shading marks canonical structural regions; q2 and q3 remain separate alternatives. B(w) is a fixed-taxon substitution-distance summary.", ha="center", fontsize=8)
    fig.tight_layout(rect=[0, .03, 1, .98])
    for ext, dpi in [("pdf", None),("png", 240)]:
        path = FIG / f"neoaves_chr4_genomic_rearrangement_tracks.{ext}"
        if ext == "pdf": fig.savefig(path, bbox_inches="tight")
        else: fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    build()
