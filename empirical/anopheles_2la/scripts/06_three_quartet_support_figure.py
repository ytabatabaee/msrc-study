#!/usr/bin/env python3
"""Create the Anopheles 2La five-panel three-topology support figure."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc-paper-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical/anopheles_2la"
RES = EMP / "results"
FIG = EMP / "figures"
OUT = RES / "three_quartet_support"
TOPO = RES / "stage2r_quartet_topology_signal.tsv"
DIST = RES / "stage12_distance/stage12_window_distance.tsv"
CROSS = RES / "stage2r_crossed_distance_signal.tsv"
MARGIN = RES / "stage2r_quartet_margin_signal.tsv"
RECOMB = ROOT / "data/anopheles_2la/processed/2la_recombination_summary.tsv"
INV_START, INV_END = 20_524_058, 42_165_532
SEED = 20261020
COLORS = {"q_species": "#4c78a8", "q_arrangement": "#f58518", "q_third": "#54a24b"}


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        out = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        out.writeheader()
        out.writerows({key: row.get(key, "") for key in fields} for row in rows)


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def finite(value):
    return math.isfinite(number(value))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def segments(x):
    if len(x) < 2:
        return [np.arange(len(x))]
    gap = max(float(np.median(np.diff(x))) * 2.5, 0.01)
    return [part for part in np.split(np.arange(len(x)), np.flatnonzero(np.diff(x) > gap) + 1) if len(part)]


def shade(axis):
    axis.axvspan(INV_START / 1e6, INV_END / 1e6, color="#888888", alpha=.14, zorder=0)
    axis.axvline(INV_START / 1e6, color="#555555", ls="--", lw=.7)
    axis.axvline(INV_END / 1e6, color="#555555", ls="--", lw=.7)


def track(axis, x, values, color, label, lw=.8, marker="."):
    y = np.array([number(value) for value in values])
    for part in segments(x):
        valid = np.isfinite(y[part])
        if not valid.any():
            continue
        axis.plot(x[part][valid], y[part][valid], color=color, lw=lw,
                  marker=marker, ms=1.5, label=label if part[0] == 0 else None)


def save(fig, stem):
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def main_figure(rows):
    x = np.array([number(row["midpoint"]) / 1e6 for row in rows])
    fig, axes = plt.subplots(5, 1, figsize=(13, 11), sharex=True,
                             gridspec_kw={"height_ratios": [0.9, 1.35, 0.9, 1, 1]},
                             constrained_layout=True)
    for axis in axes:
        shade(axis); axis.grid(axis="y", alpha=.16); axis.set_xlim(0, 49.364325)
    # The source is a published regional-summary schematic, not a digitized map.
    axes[0].plot([0, INV_START / 1e6], [2, 2], color="#555555", lw=2, label="~2.0 cM/Mb homokaryotype background")
    axes[0].plot([INV_END / 1e6, 49.364325], [2, 2], color="#555555", lw=2)
    axes[0].plot([INV_START / 1e6, INV_END / 1e6], [.5, .5], color="#8c6bb1", lw=3, label="<0.5 cM/Mb 2La/+a")
    axes[0].text((INV_START + INV_END) / 2e6, .66, "published regional estimate", ha="center", fontsize=8)
    axes[0].set_ylim(0, 2.35); axes[0].set_ylabel("cM/Mb"); axes[0].set_title("(a) published recombination-rate reference", loc="left", weight="bold"); axes[0].legend(frameon=False, fontsize=8, loc="upper right")
    base = np.zeros(len(rows))
    for key, label in (("q_species", "q_S species history"), ("q_arrangement", "q_H arrangement-associated"), ("q_third", "q_3 third")):
        y = np.array([number(row[key]) for row in rows])
        for part in segments(x):
            valid = np.isfinite(y[part])
            if valid.any():
                axes[1].fill_between(x[part][valid], base[part][valid], base[part][valid] + y[part][valid], color=COLORS[key], alpha=.84, label=label if part[0] == 0 else None)
        base += np.nan_to_num(y)
    axes[1].set_ylim(0, 1); axes[1].set_ylabel("support"); axes[1].set_title("(b) three-topology quartet support", loc="left", weight="bold"); axes[1].legend(frameon=False, fontsize=8, ncol=3, loc="upper right")
    track(axes[2], x, [row["D"] for row in rows], "#333333", "D(w) = q_H - q_S")
    axes[2].axhline(0, color="#777", lw=.6); axes[2].set_ylabel("D(w)"); axes[2].set_title("(c) quartet topology contrast", loc="left", weight="bold")
    track(axes[3], x, [row["A_distance"] for row in rows], "#c44e52", "A_distance")
    axes[3].axhline(0, color="#777", lw=.6); axes[3].set_ylabel("NJ units"); axes[3].set_title("(d) absolute divergence-distance contrast", loc="left", weight="bold")
    track(axes[4], x, [row["A_relative"] for row in rows], "#55a868", "A_relative")
    axes[4].axhline(0, color="#777", lw=.6); axes[4].set_ylabel("relative contrast"); axes[4].set_xlabel("2L position (Mb)"); axes[4].set_title("(e) relative divergence-distance contrast", loc="left", weight="bold")
    fig.suptitle("Anopheles 2La: recombination, quartet support, and divergence-distance tracks", fontsize=15, weight="bold")
    save(fig, "anopheles_2la_three_quartet_support")


def full_figure(rows):
    x = np.array([number(row["midpoint"]) / 1e6 for row in rows])
    fig, axes = plt.subplots(5, 1, figsize=(13, 11), sharex=True, constrained_layout=True)
    for axis in axes: shade(axis); axis.grid(axis="y", alpha=.16); axis.set_xlim(0, 49.364325)
    base = np.zeros(len(rows))
    for key, label in (("q_species", "q_S"), ("q_arrangement", "q_H"), ("q_third", "q_3")):
        y = np.array([number(row[key]) for row in rows])
        for part in segments(x):
            valid = np.isfinite(y[part])
            if valid.any(): axes[0].fill_between(x[part][valid], base[part][valid], base[part][valid] + y[part][valid], color=COLORS[key], alpha=.84, label=label if part[0] == 0 else None)
        base += np.nan_to_num(y)
    axes[0].set_ylim(0, 1); axes[0].set_ylabel("support"); axes[0].set_title("(a) three-topology support", loc="left", weight="bold"); axes[0].legend(frameon=False, ncol=3, fontsize=8)
    for axis, key, color, label, ylabel, title in ((axes[1], "D", "#333333", "D(w)", "D(w)", "(b) quartet topology contrast"), (axes[2], "A_distance", "#c44e52", "A_distance", "NJ units", "(c) absolute divergence-distance contrast"), (axes[3], "A_relative", "#55a868", "A_relative", "relative contrast", "(d) relative divergence-distance contrast")):
        track(axis, x, [row[key] for row in rows], color, label); axis.axhline(0, color="#777", lw=.6); axis.set_ylabel(ylabel); axis.set_title(title, loc="left", weight="bold")
    axes[4].set_visible(False)
    axes[3].set_xlabel("2L position (Mb)")
    fig.suptitle("Anopheles 2La topology and divergence-distance tracks", fontsize=15, weight="bold")
    save(fig, "anopheles_2la_full_topology_distance")


def supplementary(rows, cross, margin):
    x = np.array([number(row["midpoint"]) / 1e6 for row in rows])
    fig, axes = plt.subplots(3, 1, figsize=(13, 7), sharex=True, constrained_layout=True)
    for axis in axes: shade(axis); axis.grid(axis="y", alpha=.16); axis.set_xlim(0, 49.364325)
    for axis, values, color, label, ylabel in ((axes[0], [row.get("C", "nan") for row in rows], "#4c78a8", "C(w)", "C(w)"), (axes[1], [row.get("M", "nan") for row in rows], "#f58518", "M(w)", "M(w)"), (axes[2], [row.get("n_post_filter_snps", "nan") for row in rows], "#54a24b", "retained SNPs", "SNP count")):
        track(axis, x, values, color, label); axis.set_ylabel(ylabel); axis.legend(frameon=False, fontsize=8)
    axes[2].set_xlabel("2L position (Mb)")
    fig.suptitle("Anopheles 2La supplementary genealogy controls", fontsize=14, weight="bold")
    save(fig, "anopheles_2la_supplementary_genealogy_controls")


def main():
    topo = {row["window_id"]: row for row in read(TOPO)}
    dist = {row["window_id"]: row for row in read(DIST)}
    cross = {row["window_id"]: row for row in read(CROSS)}
    margin = {row["window_id"]: row for row in read(MARGIN)}
    ids = sorted(set(topo) | set(dist), key=lambda key: int((topo[key] if key in topo else dist[key])["start"]))
    rows = []
    validation = []
    for wid in ids:
        t = topo.get(wid, {}); d = dist.get(wid, {}); c = cross.get(wid, {}); m = margin.get(wid, {})
        qsum = sum(number(t.get(key, "nan")) for key in ("q_species", "q_arrangement", "q_third"))
        row = {"window_id": wid, "start": t.get("start", d.get("start")), "end": t.get("end", d.get("end")), "midpoint": d.get("midpoint", (number(t.get("start")) + number(t.get("end"))) / 2), "region_class": t.get("region_class", d.get("region_class", "")), "n_informative_quartets": t.get("n_informative_quartets", "0"), "q_species": t.get("q_species", "nan"), "q_arrangement": t.get("q_arrangement", "nan"), "q_third": t.get("q_third", "nan"), "q_unresolved": 0.0 if math.isfinite(qsum) else math.nan, "D_original": t.get("D", "nan"), "D_reconstructed": number(t.get("q_arrangement", "nan")) - number(t.get("q_species", "nan")), "D": t.get("D", "nan"), "A_distance": d.get("A_distance", "nan"), "A_relative": d.get("A_relative", "nan"), "C": c.get("C", "nan"), "M": m.get("mean_M", "nan"), "n_post_filter_snps": d.get("n_post_filter_snps", "nan")}
        rows.append(row)
        validation.append({"window_id": wid, "n_informative_quartets": row["n_informative_quartets"], "q_sum": qsum, "normalization_ok": math.isfinite(qsum) and abs(qsum - 1) < 1e-8, "D_original": row["D_original"], "D_reconstructed": row["D_reconstructed"], "D_difference": number(row["D_original"]) - number(row["D_reconstructed"]), "stage12_distance_available": finite(row["A_distance"])})
    rows.sort(key=lambda row: int(row["start"]))
    write(OUT / "window_three_quartet_support.tsv", rows, list(rows[0]))
    write(OUT / "D_consistency_validation.tsv", validation, list(validation[0]))
    valid = [row for row in rows if finite(row["q_species"]) and finite(row["A_distance"])]
    summary = []
    for label, selector in (("inside_2La", lambda row: INV_START <= int(row["start"]) and int(row["end"]) <= INV_END), ("outside_2La", lambda row: int(row["end"]) < INV_START or int(row["start"]) > INV_END)):
        group = [row for row in valid if selector(row)]
        for key in ("q_species", "q_arrangement", "q_third", "D_original", "A_distance", "A_relative"):
            vals = [number(row[key]) for row in group if finite(row[key])]
            summary.append({"region": label, "metric": key, "n": len(vals), "mean": float(np.mean(vals)) if vals else math.nan, "median": float(np.median(vals)) if vals else math.nan})
    write(OUT / "inside_outside_support_summary.tsv", summary, list(summary[0]))
    main_figure(rows); full_figure(rows); supplementary(rows, cross, margin)
    (OUT / "quartet_definition_audit.md").write_text("""# Quartet definition audit

The primary support values are the frozen Stage 2R local-quartet fractions from `stage2r_quartet_topology_signal.tsv`. The three mutually exclusive resolved resolutions are `q_species` (species-history topology), `q_arrangement` (the previously defined arrangement-associated alternative), and `q_third` (the remaining quartet resolution). They use the same 2,304 informative quartet design in valid windows; windows with zero informative quartets remain missing.

The preferred contrast is `D = q_arrangement - q_species`. No arrangement labels were inferred from local topology; the arrangement designation is inherited from the frozen quartet design. The Stage 12 distance tracks are separate local NJ patristic-distance summaries and are not calibrated coalescence times.

The recombination panel is the existing conservative Stump et al. regional-summary schematic. Marker-level physical positions were not recovered, so the flanking estimates are not drawn as interpolated genomic curves.
""")
    max_d = max(abs(number(row["D_difference"])) for row in validation if math.isfinite(number(row["D_difference"])))
    normalization = sum(1 for row in validation if row["normalization_ok"])
    (OUT / "validation_report.md").write_text(f"""# Validation report

- Union windows: {len(rows)}; Stage 2R topology windows: {len(topo)}; Stage 12 distance windows: {len(dist)}.
- Valid normalized three-support windows: {normalization}; maximum finite q-sum deviation from one: {max(abs(number(row['q_sum'])-1) for row in validation if math.isfinite(number(row['q_sum']))):.3g}.
- Maximum finite difference between saved D and q_arrangement-q_species: {max_d:.3g}.
- Stage 12 A_distance/A_relative values were copied without recomputation.
- Inversion coordinates: 2L:{INV_START}-{INV_END}.
- Missing windows and zero-informative-quartet rows remain missing; no interpolation was introduced.
""")
    (OUT / "figure_caption.md").write_text("""**Figure. Anopheles 2La recombination, quartet support, and divergence-distance tracks.** Panel A retains the published Stump et al. regional recombination estimates as the existing schematic reference. Panel B shows the three resolved quartet supports from the frozen Stage 2R design: species history, arrangement-associated alternative, and third resolution. Panel C is D(w)=q_H-q_S. Panels D and E are the frozen Stage 12 local-NJ patristic-distance contrasts in NJ substitution units. Gray shading marks 2La at 20.524058–42.165532 Mb; dashed lines mark its boundaries. Missing windows remain disconnected.
""")
    (OUT / "methods_text.md").write_text("""# Methods

The figure joins the frozen Stage 2R quartet support table with the frozen Stage 12 distance table by window ID. The support fractions are not recomputed from incompatible summaries: q_S, q_H, and q_3 are the saved q_species, q_arrangement, and q_third values from the same local quartet design. The preferred contrast is q_H-q_S. Absolute and relative distance tracks are copied from Stage 12 and retain their NJ substitution-distance interpretation. The published recombination reference is the existing Stump et al. regional-summary schematic.
""")
    manifest = {"stage":"three_quartet_support","seed":SEED,"command":"python empirical/anopheles_2la/scripts/06_three_quartet_support_figure.py","inputs":{str(path.relative_to(ROOT)):sha(path) for path in (TOPO,DIST,CROSS,MARGIN,RECOMB)},"inversion_start_bp":INV_START,"inversion_end_bp":INV_END,"existing_outputs_preserved":True}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"n_union_windows":len(rows),"n_valid_support_windows":normalization,"max_D_difference":max_d}, indent=2))


if __name__ == "__main__":
    main()
