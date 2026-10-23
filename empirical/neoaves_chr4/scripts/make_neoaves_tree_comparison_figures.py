#!/usr/bin/env python3
"""Make auditable Neoaves 48/363 tree and branch-length comparison figures.

This script only reads completed ASTRAL outputs.  It never runs ASTRAL or
reconstructs gene trees.  Trees are pruned after inference for display; the
branch comparison tables are calculated from the complete inferred trees.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import dendropy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "stage4e"
CLUSTER = RESULTS / "cluster"
FIGURES = ROOT / "figures"

SCALES = {
    "48": ["J48_T0", "J48_T_PNAS", "J48_T_STRUCT_PRIMARY"],
    "363": ["T0_extracted", "T_PNAS", "T_STRUCT_PRIMARY"],
}
ALL_TREATMENTS = {
    "48": ["J48_T0", "J48_T_PNAS", "J48_T_STRUCT_PRIMARY", "J48_T_BLOCK_500K", "J48_NONCHR4"],
    "363": ["T0_extracted", "T_PNAS", "T_STRUCT_PRIMARY", "T_BLOCK_500K", "NONCHR4"],
}
DISPLAY_TAXA = [
    "Acanthisitta_chloris", "Chlamydotis_macqueenii", "Columba_livia",
    "Cuculus_canorus", "Eurypyga_helias", "Gallus_gallus",
    "Mesitornis_unicolor", "Phoenicopterus_ruber", "Podiceps_cristatus",
    "Pterocles_gutturalis", "Tauraco_erythrolophus", "Tyto_alba",
]
DISPLAY_LABELS = {
    "Acanthisitta_chloris": "Acanthisitta",
    "Chlamydotis_macqueenii": "Chlamydotis",
    "Columba_livia": "Columba",
    "Cuculus_canorus": "Cuculus",
    "Eurypyga_helias": "Eurypyga",
    "Gallus_gallus": "Gallus",
    "Mesitornis_unicolor": "Mesitornis",
    "Phoenicopterus_ruber": "Phoenicopterus",
    "Podiceps_cristatus": "Podiceps",
    "Pterocles_gutturalis": "Pterocles",
    "Tauraco_erythrolophus": "Tauraco",
    "Tyto_alba": "Tyto",
}
COLORS = {"N61": "#2a6fbb", "N62 baseline": "#238b45", "N62 alternative": "#c0397a", "Columbea": "#d97706"}
TREATMENT_LABELS = {
    "J48_T0": "T0", "J48_T_PNAS": "PNAS", "J48_T_STRUCT_PRIMARY": "Structural",
    "J48_T_BLOCK_500K": "500-kb", "J48_NONCHR4": "Non-chr4",
    "T0_extracted": "T0", "T_PNAS": "PNAS", "T_STRUCT_PRIMARY": "Structural",
    "T_BLOCK_500K": "500-kb", "NONCHR4": "Non-chr4",
}


def split_key(tree: dendropy.Tree, node: dendropy.Node) -> str:
    all_taxa = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    side = {leaf.taxon.label for leaf in node.leaf_iter()}
    other = all_taxa - side
    chosen = side if (len(side), sorted(side)) <= (len(other), sorted(other)) else other
    return "|".join(sorted(chosen))


def tree_splits(tree: dendropy.Tree) -> dict[str, dendropy.Node]:
    return {split_key(tree, n): n for n in tree.postorder_node_iter() if n.parent_node and n.taxon is None}


def load_tree(treatment: str, suffix: str) -> dendropy.Tree:
    path = CLUSTER / treatment / f"tree.{suffix}.nwk"
    return dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True)


def read_branch_table(treatment: str) -> dict[str, dict[str, str]]:
    with (CLUSTER / treatment / "branch_lengths.tsv").open() as fh:
        return {row["split_key"]: row for row in csv.DictReader(fh, delimiter="\t")}


def fmt(v: float | None) -> str:
    return "NA" if v is None or not math.isfinite(v) else f"{v:.10g}"


def numeric(v: str | None) -> float | None:
    if not v or v == "NA":
        return None
    try:
        return float(v)
    except ValueError:
        return None


def write_branch_comparisons(out: Path) -> tuple[dict, dict]:
    rows = []
    summary = defaultdict(dict)
    for scale, treatments in ALL_TREATMENTS.items():
        base = treatments[0]
        base_table = read_branch_table(base)
        base_internal = {k: v for k, v in base_table.items() if v["branch_type"] == "internal"}
        for treatment in treatments[1:]:
            table = read_branch_table(treatment)
            shared = sorted(set(base_internal) & {k for k, v in table.items() if v["branch_type"] == "internal"})
            for metric in ("CU_length", "SU_length"):
                vals = []
                for key in shared:
                    x, y = numeric(base_internal[key].get(metric)), numeric(table[key].get(metric))
                    if x is None or y is None:
                        continue
                    delta = y - x
                    vals.append(abs(delta))
                    rows.append({"scale": scale, "baseline": base, "treatment": treatment, "metric": metric,
                                 "split_key": key, "baseline_value": fmt(x), "treatment_value": fmt(y),
                                 "delta": fmt(delta), "shared_internal_split": "1"})
                summary[(scale, treatment, metric)] = {
                    "shared": len(vals), "median_abs": sorted(vals)[len(vals)//2] if vals else None,
                    "mean_abs": sum(vals)/len(vals) if vals else None,
                }
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    return summary, {scale: {t: read_branch_table(t) for t in ts} for scale, ts in ALL_TREATMENTS.items()}


def write_focal_table(out: Path, tables: dict) -> None:
    # These are explicit display roles, defined by the changed full-tree splits.
    focal = {"N61": "N61", "N62": "N62", "Columbea": "Columbea"}
    rows = []
    for scale, treatments in ALL_TREATMENTS.items():
        # N61 is stable in both scales.  The best reproducible focal records
        # are recovered from the existing audited branch table when available.
        existing = RESULTS / ("j48_focal_branch_comparison.tsv" if scale == "48" else "full363_focal_branch_lengths.tsv")
        existing_rows = []
        if existing.exists():
            with existing.open() as fh:
                existing_rows = list(csv.DictReader(fh, delimiter="\t"))
        for treatment in treatments:
            tree = load_tree(treatment, "cu")
            table = tables[scale][treatment]
            for role in focal:
                source_treatment = treatment.replace("_extracted", "")
                matching = [r for r in existing_rows if r.get("treatment") == source_treatment and (r.get("focal") or r.get("clade")) == role]
                # Existing tables are authoritative for focal annotations.
                if scale == "48":
                    rec = matching[0] if matching else {}
                    present = rec.get("branch_present", "false").lower() == "true"
                    observed = rec.get("observed_split_key", "") if present else ""
                    state = rec.get("state", "absent")
                else:
                    matching = [r for r in matching if r.get("present", "False").lower() == "true"]
                    rec = matching[0] if matching else {}
                    present = bool(matching)
                    observed = ""
                    if present:
                        matches = [k for k, v in table.items() if v.get("branch_type") == "internal"
                                   and v.get("CU_length") == rec.get("CU_length")
                                   and v.get("SU_length") == rec.get("SU_length")]
                        if len(matches) == 1: observed = matches[0]
                    state = rec.get("state", "absent")
                row = {"scale": scale, "treatment": treatment, "focal": role,
                       "state": state, "present": "1" if present else "0", "split_key": observed,
                       "CU_length": rec.get("CU_length", "NA"), "SU_length": rec.get("SU_length", "NA"),
                       "localPP": rec.get("localPP", "NA"), "q1": rec.get("q1", "NA"),
                       "q2": rec.get("q2", "NA"), "q3": rec.get("q3", "NA"),
                       "source_table": str(existing) if existing.exists() else "branch_lengths.tsv"}
                if present and observed in table:
                    row["CU_length"] = table[observed].get("CU_length", row["CU_length"])
                    row["SU_length"] = table[observed].get("SU_length", row["SU_length"])
                rows.append(row)
    fields = list(rows[0])
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t"); w.writeheader(); w.writerows(rows)


def prune_tree(tree: dendropy.Tree, labels: list[str]) -> dendropy.Tree:
    pruned = tree.clone(depth=2)
    pruned.retain_taxa_with_labels(labels)
    pruned.suppress_unifurcations(update_bipartitions=False)
    return pruned


def display_mapping(scale: str, treatments: list[str], out: Path) -> None:
    rows = []
    for treatment in treatments:
        topo = load_tree(treatment, "cu")
        cu = load_tree(treatment, "cu")
        su = load_tree(treatment, "su")
        source_table = read_branch_table(treatment)
        labels = [x for x in DISPLAY_TAXA if x in {l.taxon.label for l in topo.leaf_node_iter()}]
        display = prune_tree(topo, labels)
        full_splits = tree_splits(topo)
        for node in display.preorder_node_iter():
            if node.parent_node is None: continue
            if node.taxon is not None:
                source = node.taxon.label
                kind = "terminal_cu_unavailable"
                source_ids = source
                dkey = source
            else:
                dkey = split_key(display, node)
                candidates = []
                dside = {x.taxon.label for x in node.leaf_iter()}
                for key in full_splits:
                    projected = set(key.split("|")) & set(labels)
                    comp = set(labels) - projected
                    canonical = "|".join(sorted(projected if (len(projected), sorted(projected)) <= (len(comp), sorted(comp)) else comp))
                    if canonical == dkey: candidates.append(key)
                kind = "exact_internal" if len(candidates) == 1 else "composite_internal_path"
                source_ids = ";".join(candidates)
                source = dkey
            cu_len = node.edge_length
            su_len = None
            if node.taxon is not None:
                cu_len = None
                su_len = source_table.get(node.taxon.label, {}).get("SU_length")
            else:
                # DendroPy pruning sums the original branch lengths; retain
                # that value in the mapping and separately mark composite paths.
                su_node = next(n for n in prune_tree(su, labels).preorder_node_iter()
                               if (n.taxon.label if n.is_leaf() else split_key(prune_tree(su, labels), n)) == dkey)
                su_len = su_node.edge_length
            rows.append({"scale": scale, "treatment": treatment, "display_edge": dkey,
                         "source_edge_kind": kind, "source_split_or_taxon": source,
                         "source_edge_ids": source_ids, "display_CU_length": fmt(cu_len),
                         "display_SU_length": fmt(numeric(su_len) if isinstance(su_len, str) else su_len), "taxa_retained": len(labels)})
    mode = "a" if out.exists() else "w"
    with out.open(mode, newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)


def draw_tree(ax, tree: dendropy.Tree, metric: str, focal_keys: dict[str, set[str]], title: str) -> None:
    leaves = list(tree.leaf_node_iter())
    y = {leaf.taxon.label: i for i, leaf in enumerate(reversed(leaves))}
    internal_lengths = [n.edge_length for n in tree.postorder_node_iter() if n.parent_node and not n.is_leaf() and n.edge_length]
    schematic_terminal = max(internal_lengths or [1.0]) * 0.08
    def rec(node, x0):
        if node.taxon is not None:
            yy = y[node.taxon.label]; x1 = x0 + (node.edge_length or 0.0)
            if metric == "CU": x1 = x0 + schematic_terminal
            ax.plot([x0, x1], [yy, yy], color="#777", lw=0.8, ls="--" if metric == "CU" else "-")
            ax.text(x1 + 0.01, yy, DISPLAY_LABELS.get(node.taxon.label, node.taxon.label), va="center", fontsize=6)
            return yy, x1
        children = [rec(c, x0 + (node.edge_length or 0.0)) for c in node.child_node_iter()]
        ys = [v[0] for v in children]; yy = sum(ys)/len(ys); x1 = x0 + (node.edge_length or 0.0)
        ax.plot([x1, x1], [min(ys), max(ys)], color="#333", lw=0.85)
        for c, (cy, cx) in zip(node.child_node_iter(), children):
            color = "#333"
            k = split_key(tree, c)
            for role, keys in focal_keys.items():
                if k in keys: color = COLORS.get(role, "#333")
            ax.plot([x1, cx], [cy, cy], color=color, lw=2.0 if color != "#333" else 0.85)
            if color != "#333" and c.edge_length is not None:
                ax.text((x1 + cx) / 2, cy + 0.18, f"{c.edge_length:.3g}", color=color,
                        fontsize=5.5, ha="center", va="bottom")
        return yy, x1
    rec(tree.seed_node, 0.0)
    ax.set_title(title, fontsize=8, loc="left", pad=3)
    ax.set_yticks([]); ax.set_xticks([]); ax.spines[:].set_visible(False)


def make_tree_figure(metric: str, outstem: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), squeeze=False)
    for i, (scale, treatments) in enumerate(SCALES.items()):
        # Highlight all changed splits seen between T0 and PNAS after projection.
        base = load_tree(treatments[0], "cu"); alt = load_tree(treatments[1], "cu")
        changed = set(tree_splits(base)) ^ set(tree_splits(alt))
        focal = {"N62 baseline": set(), "N62 alternative": set()}
        display_set = set(DISPLAY_TAXA)
        for key in tree_splits(base):
            if key in changed:
                side = set(key.split("|")) & display_set; other = display_set - side
                focal["N62 baseline"].add("|".join(sorted(side if (len(side), sorted(side)) <= (len(other), sorted(other)) else other)))
        for key in tree_splits(alt):
            if key in changed:
                side = set(key.split("|")) & display_set; other = display_set - side
                focal["N62 alternative"].add("|".join(sorted(side if (len(side), sorted(side)) <= (len(other), sorted(other)) else other)))
        for j, treatment in enumerate(treatments):
            tree = prune_tree(load_tree(treatment, "cu" if metric == "CU" else "su"), DISPLAY_TAXA)
            draw_tree(axes[i][j], tree, metric, focal, f"{'48 taxa' if scale == '48' else '363 taxa'}  ·  {TREATMENT_LABELS[treatment]}")
            if j == 0: axes[i][j].set_ylabel("display subset", fontsize=7)
    axes[0][0].legend(handles=[Line2D([0],[0],color=COLORS['N62 baseline'],lw=2,label='baseline focal split'),
                                Line2D([0],[0],color=COLORS['N62 alternative'],lw=2,label='alternative focal split'),
                                Line2D([0],[0],color='#777',lw=0.8,ls='--',label='CU terminal schematic' if metric=='CU' else 'terminal branch')],
                       loc='lower left', bbox_to_anchor=(0,-0.22), fontsize=7, frameon=False, ncol=3)
    fig.suptitle(f"Neoaves species-tree comparison ({metric} branch lengths)", fontsize=12, y=.995)
    fig.tight_layout(rect=(0, .03, 1, .98))
    for ext in ("pdf", "png", "svg"):
        fig.savefig(outstem.with_suffix(f".{ext}"), dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def make_branch_figure(summary: dict, outstem: Path, tables: dict) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    colors = {
        "J48_T_PNAS": "#c0397a", "J48_T_STRUCT_PRIMARY": "#238b45", "J48_T_BLOCK_500K": "#777777", "J48_NONCHR4": "#d97706",
        "T_PNAS": "#c0397a", "T_STRUCT_PRIMARY": "#238b45", "T_BLOCK_500K": "#777777", "NONCHR4": "#d97706",
    }
    for r, metric in enumerate(("CU_length", "SU_length")):
        for c, scale in enumerate(("48", "363")):
            ax = axes[r][c]; base = ALL_TREATMENTS[scale][0]; base_tab=tables[scale][base]
            plotted=0
            for treatment in ALL_TREATMENTS[scale][1:]:
                tab=tables[scale][treatment]
                xs=[]; ys=[]
                for key,row in base_tab.items():
                    if row['branch_type']!='internal' or key not in tab: continue
                    x=numeric(row.get(metric)); y=numeric(tab[key].get(metric))
                    if x is not None and y is not None: xs.append(x); ys.append(y)
                if xs:
                    ax.scatter(xs,ys,s=7,alpha=.45,label=TREATMENT_LABELS[treatment],color=colors[treatment]); plotted+=len(xs)
            lim=max([v for coll in ax.collections for v in coll.get_offsets().data.ravel()] or [1])
            ax.plot([0,lim],[0,lim],color='#333',lw=.7,ls=':')
            ax.set_title(f"{'48' if scale=='48' else '363'} taxa · {metric.replace('_length','')}",fontsize=9)
            ax.set_xlabel('T0'); ax.set_ylabel('treatment')
            ax.text(.04,.94,f"shared observations: {plotted}",transform=ax.transAxes,fontsize=7,va='top')
            if r==0 and c==1: ax.legend(fontsize=7,frameon=False,loc='lower right')
    fig.suptitle('Homologous internal branch-length comparisons',fontsize=12)
    fig.tight_layout(rect=(0,0,1,.97))
    for ext in ('pdf','png','svg'):
        fig.savefig(outstem.with_suffix(f'.{ext}'),dpi=300 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)


def main() -> None:
    parser=argparse.ArgumentParser(); parser.add_argument('--outdir',type=Path,default=RESULTS); parser.add_argument('--figdir',type=Path,default=FIGURES)
    args=parser.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True); args.figdir.mkdir(parents=True,exist_ok=True)
    summary,tables=write_branch_comparisons(args.outdir/'neoaves_cross_taxon_branch_lengths.tsv')
    write_focal_table(args.outdir/'neoaves_cross_taxon_focal_splits.tsv',tables)
    mapping_path = args.outdir/'neoaves_display_tree_edge_mapping.tsv'
    if mapping_path.exists(): mapping_path.unlink()
    for scale,ts in ALL_TREATMENTS.items(): display_mapping(scale,ts,args.outdir/'neoaves_display_tree_edge_mapping.tsv')
    lines=['# Neoaves cross-taxon branch-length summary','', 'Branch lengths are compared only for exact full-tree unrooted split keys shared by T0 and a treatment. Missing splits are topology changes and are not assigned zero length.','']
    for (scale,treatment,metric),v in summary.items(): lines.append(f"- {scale} taxa, {TREATMENT_LABELS[treatment]}, {metric}: {v['shared']} shared; median absolute difference {fmt(v['median_abs'])}; mean absolute difference {fmt(v['mean_abs'])}.")
    (args.outdir/'neoaves_branch_length_summary.md').write_text('\n'.join(lines)+'\n')
    make_tree_figure('CU',args.figdir/'neoaves_48_363_trees_cu'); make_tree_figure('SU',args.figdir/'neoaves_48_363_trees_su'); make_branch_figure(summary,args.figdir/'neoaves_48_363_branch_length_comparison',tables)


if __name__=='__main__': main()
