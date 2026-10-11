#!/usr/bin/env python3
"""Create the preliminary four-panel Mai30 scaled-WF figure."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT = Path(__file__).resolve()
STUDY_ROOT = SCRIPT.parents[2]
RESULTS = STUDY_ROOT / "analysis/results/mai30_scaled_wf_astral.tsv"
OUT = STUDY_ROOT / "analysis/figures/mai30_scaled_wf_pilot.png"


@dataclass
class Node:
    name: str | None = None
    length: float = 0.0
    children: list["Node"] = field(default_factory=list)
    y: float = 0.0
    x: float = 0.0


def parse_tree(text: str) -> Node:
    tokens = re.findall(r"[(),:;]|[^(),:;\s]+", text.strip())
    pos = 0

    def rec() -> Node:
        nonlocal pos
        if tokens[pos] == "(":
            pos += 1
            children = [rec()]
            while tokens[pos] == ",":
                pos += 1
                children.append(rec())
            pos += 1
            if pos < len(tokens) and tokens[pos] not in ",);:":
                pos += 1
            node = Node(children=children)
        else:
            node = Node(name=tokens[pos])
            pos += 1
        if pos < len(tokens) and tokens[pos] == ":":
            pos += 1
            node.length = float(tokens[pos])
            pos += 1
        return node

    return rec()


def layout(root: Node) -> list[Node]:
    tips = []

    def rec(node: Node, parent_x: float) -> None:
        node.x = parent_x + node.length
        if not node.children:
            node.y = len(tips)
            tips.append(node)
            return
        for child in node.children:
            rec(child, node.x)
        node.y = float(np.mean([child.y for child in node.children]))

    # Newick lengths are backward from tips; use a transformed coordinate so
    # the displayed tree grows left-to-right from the root.
    def max_height(node: Node) -> float:
        return node.length + (max((max_height(c) for c in node.children), default=0.0) if node.children else 0.0)
    height = max_height(root)

    def rec2(node: Node, age: float) -> None:
        node.x = age
        if not node.children:
            node.y = len(tips)
            tips.append(node)
            return
        for child in node.children:
            rec2(child, age - child.length)
        node.y = float(np.mean([child.y for child in node.children]))

    tips.clear()
    rec2(root, height)
    return tips


def draw_tree(ax) -> None:
    source = next((STUDY_ROOT / "datasets/processed/mai30_scaled_wf_runs").glob("*replicate_042*/history_000/species_tree.nwk"))
    root = parse_tree(source.read_text())
    tips = layout(root)
    def draw(node: Node):
        if not node.children:
            return
        ys = [child.y for child in node.children]
        ax.plot([node.x, node.x], [min(ys), max(ys)], color="#8a8a8a", lw=0.45)
        for child in node.children:
            ax.plot([node.x, child.x], [child.y, child.y], color="#8a8a8a", lw=0.45)
            draw(child)
    draw(root)
    for tip in tips:
        ax.text(tip.x + 0.15, tip.y, tip.name, fontsize=4.5, va="center")
    ax.axvline(root.x - 4.0, color="#c43c39", lw=2.0, alpha=0.8)
    ax.text(root.x - 3.9, 29.3, "representative ancestral origin\n(red guide; branch in config: N8)", color="#a22", fontsize=7)
    ax.set_title("A  Mai 30-taxon backbone and ancestral origin", fontsize=10, loc="left")
    ax.set_yticks([])
    ax.set_xlabel("scaled generations (c = 0.001)", fontsize=8)


def main() -> None:
    with RESULTS.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        row["rf"] = float(row["normalized_rf"])
        row["cu"] = float(row["mean_absolute_cu_error_common_splits"]) if row["mean_absolute_cu_error_common_splits"] else np.nan
        group = row["group"]
        if "msc_control" in group:
            row["treatment"] = "MSC"
        elif "no_suppression" in group:
            row["treatment"] = "MSRC 1.0"
        elif "moderate" in group:
            row["treatment"] = "MSRC 0.5"
        else:
            row["treatment"] = "MSRC 0.1"
        row["rep"] = group.split("__replicate_")[1].split("__")[0]
    order = {"042": 0.004, "005": 0.303, "074": 0.579, "007": 0.979}
    fig = plt.figure(figsize=(12, 8.2))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], hspace=0.36, wspace=0.28)
    ax_a = fig.add_subplot(grid[0, 0])
    draw_tree(ax_a)
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])
    treatments = ["MSC", "MSRC 1.0", "MSRC 0.5", "MSRC 0.1"]
    colors = {"MSC": "#444444", "MSRC 1.0": "#377eb8", "MSRC 0.5": "#ff8c42", "MSRC 0.1": "#d73027"}
    offsets = np.linspace(-0.24, 0.24, len(treatments))
    for ax, metric, ylabel, title in [
        (ax_b, "rf", "normalized RF error", "B  Species-tree topology accuracy"),
        (ax_c, "cu", "mean absolute CULength error", "C  Branch-length accuracy on common splits"),
    ]:
        for rep_i, rep in enumerate(["042", "005", "074", "007"]):
            base = rep_i
            for off, treatment in zip(offsets, treatments):
                vals = [float(r[metric]) for r in rows if r["rep"] == rep and r["treatment"] == treatment and np.isfinite(float(r[metric]))]
                if vals:
                    ax.scatter(base + off, np.mean(vals), color=colors[treatment], s=35, zorder=3)
        ax.set_xticks(range(4), ["low\n0.004", "moderate\n0.303", "high\n0.579", "extreme\n0.979"])
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=10, loc="left")
        ax.grid(axis="y", alpha=0.25)
    for treatment in treatments:
        ax_b.scatter([], [], color=colors[treatment], label=treatment)
    ax_b.legend(frameon=False, fontsize=7, ncol=2, loc="upper left")
    ax_d.set_title("D  MSRC effect relative to MSC", fontsize=10, loc="left")
    for rep_i, rep in enumerate(["042", "005", "074", "007"]):
        control = next((float(r["rf"]) for r in rows if r["rep"] == rep and r["treatment"] == "MSC"), np.nan)
        for treatment, color in [("MSRC 1.0", colors["MSRC 1.0"]), ("MSRC 0.5", colors["MSRC 0.5"]), ("MSRC 0.1", colors["MSRC 0.1"])]:
            vals = [float(r["rf"]) - control for r in rows if r["rep"] == rep and r["treatment"] == treatment]
            if vals:
                ax_d.scatter(rep_i + (treatments.index(treatment)-1) * 0.18, np.mean(vals), color=color, s=35)
    ax_d.axhline(0, color="black", lw=0.7)
    ax_d.set_xticks(range(4), ["low", "moderate", "high", "extreme"])
    ax_d.set_ylabel("Δ normalized RF (MSRC − MSC)")
    ax_d.grid(axis="y", alpha=0.25)
    fig.suptitle("Mai30 scaled Wright–Fisher MSRC pilot", fontsize=14)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=220, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
