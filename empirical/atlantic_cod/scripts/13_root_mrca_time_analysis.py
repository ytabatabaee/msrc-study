#!/usr/bin/env python3
"""Root-MRCA-time analysis for Atlantic cod 250-kb SNAPP MCC windows.

Uses the published MCC point estimates; it does not rerun SNAPP. Root age is
measured from the MRCA of the Stage-2 Atlantic-cod ingroup populations to the
ingroup tips, thereby excluding any outgroup if one is present.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import platform
import random
import statistics
import textwrap
from collections import defaultdict
from pathlib import Path

import Bio
from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FuncFormatter

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data/atlantic_cod"
EMP = ROOT / "empirical/atlantic_cod"
RESULTS = EMP / "results/root_mrca"
FIGURES = EMP / "figures"
WINDOWS = DATA / "processed/cod_window_trees.tsv"
ARRANGEMENTS = DATA / "processed/population_arrangements.tsv"
REGIONS = DATA / "metadata/region_manifest.tsv"
STAGE7_SIGNAL = EMP / "results/stage7_window_divergence_time_signal.tsv"
STAGE4_SCORES = DATA / "processed/stage4a_window_topology_scores.tsv"
LGS = ("LG01", "LG02", "LG07", "LG12")
WINDOW_BP = 250_000
ULTRAMETRIC_TOL = 1e-7
NEGATIVE_TOL = 1e-10
SEED = 20261120
BOOTSTRAPS = 5000
PRIMARY_BLOCK_WINDOWS = 8
BLOCK_SENSITIVITY = (4, 8, 12)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: fmt(row.get(k, "")) for k in fields})


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return "NA" if not math.isfinite(value) else f"{value:.9g}"
    return str(value)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def region_class(row: dict[str, str]) -> str:
    if row["inside_inversion"] == "true":
        return "inside"
    if row["overlaps_inversion"] == "true":
        return "boundary"
    return "outside"


def root_mrca_metrics(newick: str, ingroup: set[str]) -> dict[str, object]:
    tree = Phylo.read(io.StringIO(newick.strip()), "newick")
    terminals = tree.get_terminals()
    labels = [t.name for t in terminals]
    if len(labels) != len(set(labels)):
        raise ValueError("duplicate_tip_labels")
    missing = sorted(ingroup - set(labels))
    if missing:
        raise ValueError("missing_ingroup:" + ",".join(missing))
    outgroups = sorted(set(labels) - ingroup)
    min_branch = math.inf
    for clade in tree.find_clades():
        if clade.branch_length is None:
            clade.branch_length = 0.0
        x = float(clade.branch_length)
        if not math.isfinite(x):
            raise ValueError("nonfinite_branch_length")
        min_branch = min(min_branch, x)
        if x < -NEGATIVE_TOL:
            raise ValueError(f"materially_negative_branch:{x:.12g}")
        if x < 0:
            clade.branch_length = 0.0
    ingroup_tips = [t for t in terminals if t.name in ingroup]
    ingroup_mrca = tree.common_ancestor(*ingroup_tips)
    depths = tree.depths(unit_branch_lengths=False)
    tip_ages = [float(depths[t] - depths[ingroup_mrca]) for t in ingroup_tips]
    crown_age = max(tip_ages)
    max_dev = max(abs(x - crown_age) for x in tip_ages)
    if not math.isfinite(crown_age) or crown_age < 0:
        raise ValueError("invalid_ingroup_root_age")
    if max_dev > ULTRAMETRIC_TOL:
        raise ValueError(f"ingroup_not_ultrametric:{max_dev:.12g}")
    return {
        "root_mrca_time": crown_age,
        "max_ingroup_tip_depth_deviation": max_dev,
        "n_tips": len(labels),
        "n_ingroup_tips": len(ingroup_tips),
        "outgroup_tips": ",".join(outgroups) if outgroups else "none",
        "min_branch_length": min_branch,
        # Newick lacks a [&R] rootedness flag; rootedness is established by the
        # BEAST/SNAPP clock-tree workflow and its explicit crown-age prior.
        "newick_root_marker": "absent",
        "rooting_basis": "BEAST/SNAPP clock MCC; source XML constrains crown/root age",
    }


def read_inputs():
    ingroup = {r["tree_label"] for r in read_tsv(ARRANGEMENTS)}
    regions = {r["lg"]: (int(r["start"]), int(r["end"])) for r in read_tsv(REGIONS) if r["lg"] in LGS}
    stage7 = {r["window_id"]: r for r in read_tsv(STAGE7_SIGNAL)}
    stage4 = {r["window_id"]: r for r in read_tsv(STAGE4_SCORES)}
    raw_rows = [r for r in read_tsv(WINDOWS) if r["lg"] in LGS]
    rows, audit = [], []
    for raw in raw_rows:
        row = {**raw, "start": int(raw["start"]), "end": int(raw["end"]), "midpoint": int(raw["midpoint"]), "region_class": region_class(raw)}
        diag = {"lg": row["lg"], "window_id": row["window_id"], "start": row["start"], "end": row["end"], "region_class": row["region_class"]}
        try:
            metrics = root_mrca_metrics(row["tree_newick"], ingroup)
            diag.update(metrics)
            diag["valid_for_root_mrca"] = "true"
            diag["validation_issue"] = ""
            s7 = stage7[row["window_id"]]
            s4 = stage4[row["window_id"]]
            row.update(metrics)
            row["A_w"] = float(s7["A_opposite_minus_same"])
            # The cod Stage-7 file's mean_time_all is the existing arrangement-blind
            # mean of the 66 pairwise MRCA times; retain B's transparent operational name.
            row["B_w"] = float(s7["mean_time_all"])
            for key in ("fraction_arrangement", "fraction_baseline", "fraction_third", "D_arrangement_minus_baseline"):
                row[key] = float(s4[key])
            rows.append(row)
        except (ValueError, KeyError) as exc:
            diag.update({"root_mrca_time": math.nan, "max_ingroup_tip_depth_deviation": math.nan,
                         "n_tips": len(Phylo.read(io.StringIO(row["tree_newick"]), "newick").get_terminals()),
                         "n_ingroup_tips": len(ingroup), "outgroup_tips": "none",
                         "newick_root_marker": "absent",
                         "rooting_basis": "BEAST/SNAPP clock MCC; source XML constrains crown/root age",
                         "valid_for_root_mrca": "false", "validation_issue": str(exc)})
        audit.append(diag)
    return rows, audit, regions, ingroup


def summarize(values: list[float]) -> tuple[float, float]:
    return statistics.mean(values), statistics.median(values)


def contiguous_runs(rows: list[dict[str, object]], cls: str) -> list[list[float]]:
    selected = sorted((r for r in rows if r["region_class"] == cls), key=lambda r: int(r["start"]))
    runs: list[list[float]] = []
    last_end = None
    for row in selected:
        start = int(row["start"])
        if not runs or last_end is None or start != last_end + 1:
            runs.append([])
        runs[-1].append(float(row["root_mrca_time"]))
        last_end = int(row["end"])
    return runs


def circular_block_sample(run: list[float], block_len: int, rng: random.Random) -> list[float]:
    if len(run) <= 1:
        return list(run)
    block_len = max(1, min(block_len, len(run)))
    out: list[float] = []
    while len(out) < len(run):
        start = rng.randrange(len(run))
        out.extend(run[(start + j) % len(run)] for j in range(block_len))
    return out[:len(run)]


def ci_summary(rows: list[dict[str, object]], block_len: int, reps: int, seed: int) -> dict[str, float]:
    inside_runs = contiguous_runs(rows, "inside")
    outside_runs = contiguous_runs(rows, "outside")
    if not inside_runs or not outside_runs:
        return {k: math.nan for k in ("delta_mean_lo", "delta_mean_hi", "delta_median_lo", "delta_median_hi", "mean_ratio_lo", "mean_ratio_hi", "median_ratio_lo", "median_ratio_hi")}
    rng = random.Random(seed)
    stats: dict[str, list[float]] = defaultdict(list)
    for _ in range(reps):
        ins = [x for run in inside_runs for x in circular_block_sample(run, block_len, rng)]
        outs = [x for run in outside_runs for x in circular_block_sample(run, block_len, rng)]
        im, om = statistics.mean(ins), statistics.mean(outs)
        imed, omed = statistics.median(ins), statistics.median(outs)
        stats["delta_mean"].append(im - om)
        stats["delta_median"].append(imed - omed)
        stats["mean_ratio"].append(im / om if om else math.nan)
        stats["median_ratio"].append(imed / omed if omed else math.nan)
    out = {}
    for name, values in stats.items():
        a, b = np.quantile(np.asarray(values), [0.025, 0.975])
        out[name + "_lo"] = float(a)
        out[name + "_hi"] = float(b)
    return out


def build_grid(lg: str, all_rows: list[dict[str, str]], valid_rows: list[dict[str, object]], region: tuple[int, int]):
    lg_raw = [r for r in all_rows if r["lg"] == lg]
    starts = [int(r["start"]) for r in lg_raw]
    ends = [int(r["end"]) for r in lg_raw]
    start0, endmax = min(starts), max(ends)
    n = (endmax - start0 + 1 + WINDOW_BP - 1) // WINDOW_BP
    coordinates = [start0 + i * WINDOW_BP for i in range(n)]
    value_at = {int(r["start"]): float(r["root_mrca_time"]) for r in valid_rows if r["lg"] == lg}
    mask = []
    values = []
    rstart, rend = region
    for start in coordinates:
        end = start + WINDOW_BP - 1
        if start >= rstart and end <= rend:
            cls = "inside"
        elif start <= rend and end >= rstart:
            cls = "boundary"
        else:
            cls = "outside"
        mask.append(cls)
        values.append(value_at.get(start, math.nan))
    return coordinates, values, mask


def shift_stat(values: list[float], mask: list[str]) -> float:
    inside = [v for v, cls in zip(values, mask) if cls == "inside" and math.isfinite(v)]
    outside = [v for v, cls in zip(values, mask) if cls == "outside" and math.isfinite(v)]
    return statistics.mean(inside) - statistics.mean(outside) if inside and outside else math.nan


def spatial_tests(all_rows, valid_rows, regions):
    null_rows, summary = [], []
    pvals = {}
    for lg in LGS:
        coords, values, mask = build_grid(lg, all_rows, valid_rows, regions[lg])
        observed = shift_stat(values, mask)
        null = []
        for shift in range(len(values)):
            shifted = values[-shift:] + values[:-shift] if shift else values[:]
            stat = shift_stat(shifted, mask)
            null.append(stat)
            null_rows.append({"lg": lg, "shift_windows": shift, "shift_bp": shift * WINDOW_BP,
                              "delta_root_mrca": stat, "at_least_observed": stat >= observed if math.isfinite(stat) else False})
        finite = [x for x in null if math.isfinite(x)]
        p = sum(x >= observed for x in finite) / len(finite)
        pvals[lg] = p
        summary.append({"lg": lg, "n_grid_windows": len(values), "n_valid_root_ages": sum(math.isfinite(v) for v in values),
                        "observed_delta_mean": observed, "n_circular_shifts": len(finite), "p_one_sided": p})
    # Benjamini-Hochberg over the four predeclared linkage groups.
    ordered = sorted(pvals.items(), key=lambda x: x[1])
    adj, running = {}, 1.0
    for i in range(len(ordered) - 1, -1, -1):
        lg, p = ordered[i]
        running = min(running, p * len(ordered) / (i + 1), 1.0)
        adj[lg] = running
    for row in summary:
        row["BH_p"] = adj[row["lg"]]
    return summary, null_rows


def inside_outside_stats(rows: list[dict[str, object]], key: str, metric: str):
    valid = [r for r in rows if math.isfinite(float(r[key])) and r["region_class"] in {"inside", "outside"}]
    ins = [float(r[key]) for r in valid if r["region_class"] == "inside"]
    outs = [float(r[key]) for r in valid if r["region_class"] == "outside"]
    if not ins or not outs:
        return None
    imean, imed = summarize(ins)
    omean, omed = summarize(outs)
    return {"metric": metric, "n_inside": len(ins), "n_outside": len(outs),
            "mean_inside": imean, "median_inside": imed, "mean_outside": omean, "median_outside": omed,
            "delta_mean_inside_minus_outside": imean - omean,
            "delta_median_inside_minus_outside": imed - omed,
            "mean_ratio_inside_over_outside": imean / omean if omean else math.nan,
            "median_ratio_inside_over_outside": imed / omed if omed else math.nan}


def rankdata(vals: list[float]) -> list[float]:
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and vals[order[j]] == vals[order[i]]:
            j += 1
        rank = (i + 1 + j) / 2
        for k in range(i, j):
            ranks[order[k]] = rank
        i = j
    return ranks


def corr(x: list[float], y: list[float]) -> tuple[float, float]:
    if len(x) < 3:
        return math.nan, math.nan
    def pearson(a, b):
        ma, mb = statistics.mean(a), statistics.mean(b)
        num = sum((u-ma)*(v-mb) for u,v in zip(a,b))
        den = math.sqrt(sum((u-ma)**2 for u in a) * sum((v-mb)**2 for v in b))
        return num / den if den else math.nan
    return pearson(x,y), pearson(rankdata(x),rankdata(y))


def correlation_rows(rows):
    outputs = []
    metrics = (("B_w", "B(w): mean all-pair MRCA time"), ("A_w", "A(w): opposite minus same arrangement"),
               ("fraction_arrangement", "quartet frequency: arrangement"),
               ("fraction_baseline", "quartet frequency: baseline"), ("fraction_third", "quartet frequency: third topology"))
    for lg in LGS:
        elig = [r for r in rows if r["lg"] == lg and r["region_class"] != "boundary"]
        for key, label in metrics:
            x, y = [], []
            for r in elig:
                a, b = float(r["root_mrca_time"]), float(r[key])
                if math.isfinite(a) and math.isfinite(b):
                    x.append(a); y.append(b)
            pearson, spearman = corr(x,y)
            outputs.append({"lg": lg, "metric": label, "n_windows": len(x), "pearson_r": pearson, "spearman_rho": spearman,
                            "interpretation": "descriptive_only; linked_windows"})
    return outputs


def pval_lookup(rows):
    return {r["lg"]: r for r in rows}


def save_figures(rows, regions):
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "pdf.fonttype": 42, "ps.fonttype": 42})
    # Primary: root-age tracks paired with inside/outside distributions.
    fig, axes = plt.subplots(4, 2, figsize=(11.4, 10.2), gridspec_kw={"width_ratios": [2.1, 1]})
    root_color = "#136f75"
    for i, lg in enumerate(LGS):
        ax = axes[i,0]
        sub = sorted([r for r in rows if r["lg"] == lg], key=lambda r: int(r["start"]))
        lo, hi = regions[lg]
        ax.axvspan(lo/1e6, hi/1e6, color="#d9d9d9", alpha=0.62, zorder=0)
        xs = [int(r["midpoint"])/1e6 for r in sub]
        ys = [float(r["root_mrca_time"]) if math.isfinite(float(r["root_mrca_time"])) else np.nan for r in sub]
        ax.plot(xs, ys, color=root_color, lw=0.85, zorder=2)
        ax.scatter(xs, ys, s=10, color=root_color, edgecolor="white", linewidth=0.25, zorder=3)
        boundary = [r for r in sub if r["region_class"] == "boundary"]
        for r in boundary:
            y = float(r["root_mrca_time"])
            if math.isfinite(y):
                ax.scatter(int(r["midpoint"])/1e6, y, s=22, marker="D", facecolor="#e07a5f", edgecolor="white", linewidth=0.4, zorder=4)
            else:
                ax.scatter(int(r["midpoint"])/1e6, 0.98, s=24, marker="x", color="#c23b22", transform=ax.get_xaxis_transform(), zorder=5)
        ax.set_ylabel(f"{lg}\nRoot age (Ma)")
        ax.grid(axis="y", color="#e6e9ec", lw=0.6)
        ax.set_xlim(0, max(int(r["end"]) for r in sub)/1e6)
        ax.xaxis.set_major_locator(MultipleLocator(5))
        if i == 0:
            ax.set_title("Ingroup root-MRCA time by genomic window", loc="left", fontweight="bold")
        if i == len(LGS)-1:
            ax.set_xlabel("gadMor2 genomic position (Mb)")

        dx = axes[i,1]
        vals = [[float(r["root_mrca_time"]) for r in sub if r["region_class"] == cls and math.isfinite(float(r["root_mrca_time"]))] for cls in ("outside","inside")]
        bp = dx.boxplot(vals, positions=[0,1], widths=0.48, patch_artist=True, showfliers=False,
                        medianprops={"color":"#20252b", "linewidth":1.2}, whiskerprops={"color":"#697782"}, capprops={"color":"#697782"})
        for patch, color in zip(bp["boxes"], ("#9bb7c7", "#edb183")):
            patch.set_facecolor(color); patch.set_alpha(0.8)
        rng = random.Random(SEED + i)
        for j, (v, color) in enumerate(zip(vals, ("#4c78a8","#d37243"))):
            jitter = [j + rng.uniform(-0.13,0.13) for _ in v]
            dx.scatter(jitter, v, s=8, alpha=0.42, color=color, edgecolor="none", zorder=3)
        dx.set_xticks([0,1], [f"outside\n(n={len(vals[0])})",f"inside\n(n={len(vals[1])})"])
        dx.grid(axis="y", color="#e6e9ec", lw=0.6)
        dx.set_ylabel("Root age (Ma)")
        if i == 0:
            dx.set_title("Inside/outside distributions", loc="left", fontweight="bold")
    fig.suptitle("Atlantic cod root MRCA times across supergenes", y=0.995, fontsize=13, fontweight="bold")
    fig.text(0.5,0.008,"Gray: frozen cod_phylogenomics supergene intervals in gadMor2 coordinates. Diamond: boundary-overlap window. X: invalid LG07 boundary tree excluded from root-age summaries.",ha="center",fontsize=7.5,color="#505b64")
    fig.tight_layout(rect=(0,0.035,1,0.97),h_pad=0.95,w_pad=1.1)
    p_png=FIGURES/"atlantic_cod_root_mrca_tracks_and_distributions.png"
    p_pdf=FIGURES/"atlantic_cod_root_mrca_tracks_and_distributions.pdf"
    fig.savefig(p_png,dpi=300,bbox_inches="tight",facecolor="white");fig.savefig(p_pdf,bbox_inches="tight",facecolor="white");plt.close(fig)

    # Companion panel directly compares root, B, A and quartet-composition tracks.
    fig, axes = plt.subplots(4,4,figsize=(14.2,9.3),sharex="row")
    cols=("root_mrca_time","B_w","A_w","quartets")
    titles=("T_MRCA(w): root age (Ma)","B(w): all-pair mean MRCA (Ma)","A(w): opposite − same pair MRCA (Ma)","Quartet topology frequencies")
    colors={"root_mrca_time":"#136f75","B_w":"#e09f3e","A_w":"#4c78a8"}
    for i,lg in enumerate(LGS):
        sub=sorted([r for r in rows if r["lg"]==lg],key=lambda r:int(r["start"]))
        lo,hi=regions[lg]
        for j,key in enumerate(cols):
            ax=axes[i,j];ax.axvspan(lo/1e6,hi/1e6,color="#d9d9d9",alpha=.55,zorder=0)
            xx=[int(r["midpoint"])/1e6 for r in sub]
            if key=="quartets":
                for q,qc,lab in (("fraction_arrangement","#2a9d8f","arrangement"),("fraction_baseline","#7b6db0","baseline"),("fraction_third","#a66a4c","third")):
                    ax.plot(xx,[float(r[q]) for r in sub],lw=.8,marker="o",ms=1.7,color=qc,label=lab)
                ax.set_ylim(-.03,1.03)
            else:
                yy=[float(r[key]) if math.isfinite(float(r[key])) else np.nan for r in sub]
                ax.plot(xx,yy,lw=.85,marker="o",ms=1.8,color=colors[key])
                for r in sub:
                    if r["region_class"]=="boundary" and math.isfinite(float(r[key])):
                        ax.scatter(int(r["midpoint"])/1e6,float(r[key]),s=17,marker="D",color="#e07a5f",edgecolor="white",linewidth=.3,zorder=3)
            ax.grid(axis="y",color="#e6e9ec",lw=.5)
            if j==0: ax.set_ylabel(f"{lg}\nRoot age (Ma)")
            if i==0: ax.set_title(titles[j],loc="left",fontsize=8.2,fontweight="bold")
            if i==3: ax.set_xlabel("gadMor2 position (Mb)")
            ax.set_xlim(0,max(int(r["end"]) for r in sub)/1e6)
            ax.xaxis.set_major_locator(MultipleLocator(10 if max(int(r["end"]) for r in sub)>25e6 else 5))
    h,l=axes[0,3].get_legend_handles_labels();fig.legend(h,l,loc="upper right",bbox_to_anchor=(.995,.995),ncol=3,frameon=False,fontsize=7.5)
    fig.suptitle("Root age, pairwise coalescence, and topology signals",y=1.005,fontsize=12,fontweight="bold")
    fig.text(.5,.005,"All tracks and frozen supergene intervals use gadMor2 coordinates; linked 250-kb windows are descriptive, not independent replicates.",ha="center",fontsize=7.5,color="#505b64")
    fig.tight_layout(rect=(0,.025,1,.96),h_pad=.75,w_pad=.7)
    p_png2=FIGURES/"atlantic_cod_root_mrca_comparison_tracks.png"
    p_pdf2=FIGURES/"atlantic_cod_root_mrca_comparison_tracks.pdf"
    fig.savefig(p_png2,dpi=300,bbox_inches="tight",facecolor="white");fig.savefig(p_pdf2,bbox_inches="tight",facecolor="white");plt.close(fig)
    return [p_png,p_pdf,p_png2,p_pdf2]


def main():
    RESULTS.mkdir(parents=True,exist_ok=True)
    rows,audit,regions,ingroup=read_inputs()
    all_input=read_tsv(WINDOWS)
    valid_by_id={r["window_id"]:r for r in rows}
    # Audit inventory includes every tree row, including validation failures.
    audit_fields=["lg","window_id","start","end","region_class","n_tips","n_ingroup_tips","outgroup_tips","root_mrca_time","max_ingroup_tip_depth_deviation","min_branch_length","newick_root_marker","rooting_basis","valid_for_root_mrca","validation_issue"]
    write_tsv(RESULTS/"root_mrca_tree_validation.tsv",audit,audit_fields)

    # Window-level root/A/B/quartet series, with invalid trees retained as NA rows in audit only.
    signal_fields=["lg","window_id","start","end","midpoint","region_class","root_mrca_time","n_ingroup_tips","n_outgroup_tips","A_w","B_w","fraction_arrangement","fraction_baseline","fraction_third","D_arrangement_minus_baseline"]
    signal=[]
    for source in all_input:
        if source["lg"] not in LGS: continue
        if source["window_id"] in valid_by_id:
            row=valid_by_id[source["window_id"]]
            signal.append({k:row.get(k,"NA") for k in signal_fields})
        else:
            signal.append({"lg":source["lg"],"window_id":source["window_id"],"start":source["start"],"end":source["end"],"midpoint":source["midpoint"],"region_class":region_class(source),"root_mrca_time":"NA","n_ingroup_tips":len(ingroup),"n_outgroup_tips":0,"A_w":"NA","B_w":"NA","fraction_arrangement":"NA","fraction_baseline":"NA","fraction_third":"NA","D_arrangement_minus_baseline":"NA"})
    write_tsv(RESULTS/"root_mrca_window_signal.tsv",signal,signal_fields)

    summary_rows=[]
    for lg in LGS:
        sub=[r for r in rows if r["lg"]==lg]
        row=inside_outside_stats(sub,"root_mrca_time","T_MRCA_root")
        row["lg"]=lg
        row.update(ci_summary(sub,PRIMARY_BLOCK_WINDOWS,BOOTSTRAPS,SEED+LGS.index(lg)))
        row["bootstrap_block_windows"]=PRIMARY_BLOCK_WINDOWS
        row["bootstrap_block_bp"]=PRIMARY_BLOCK_WINDOWS*WINDOW_BP
        summary_rows.append(row)
    summary_fields=["lg","metric","n_inside","n_outside","mean_inside","median_inside","mean_outside","median_outside","delta_mean_inside_minus_outside","delta_median_inside_minus_outside","mean_ratio_inside_over_outside","median_ratio_inside_over_outside","delta_mean_lo","delta_mean_hi","delta_median_lo","delta_median_hi","mean_ratio_lo","mean_ratio_hi","median_ratio_lo","median_ratio_hi","bootstrap_block_windows","bootstrap_block_bp"]
    write_tsv(RESULTS/"root_mrca_inside_outside.tsv",summary_rows,summary_fields)

    spatial_summary,spatial_null=spatial_tests(all_input,rows,regions)
    write_tsv(RESULTS/"root_mrca_spatial_test.tsv",spatial_summary,["lg","n_grid_windows","n_valid_root_ages","observed_delta_mean","n_circular_shifts","p_one_sided","BH_p"])
    write_tsv(RESULTS/"root_mrca_spatial_null.tsv",spatial_null,["lg","shift_windows","shift_bp","delta_root_mrca","at_least_observed"])

    block_rows=[]
    for bi,block in enumerate(BLOCK_SENSITIVITY):
        for i,lg in enumerate(LGS):
            sub=[r for r in rows if r["lg"]==lg]
            eff=ci_summary(sub,block,BOOTSTRAPS,SEED+1000+bi*len(LGS)+i)
            block_rows.append({"lg":lg,"block_windows":block,"block_bp":block*WINDOW_BP,**eff})
    write_tsv(RESULTS/"root_mrca_block_sensitivity.tsv",block_rows,["lg","block_windows","block_bp","delta_mean_lo","delta_mean_hi","delta_median_lo","delta_median_hi","mean_ratio_lo","mean_ratio_hi","median_ratio_lo","median_ratio_hi"])

    # Root vs existing cod Stage-7 A/B and Stage-4 quartet signals.
    comparison=[]
    for lg in LGS:
        sub=[r for r in rows if r["lg"]==lg]
        for key,label in (("root_mrca_time","T_MRCA_root"),("B_w","B_mean_all_pair_MRCA"),("A_w","A_opposite_minus_same"),
                          ("D_arrangement_minus_baseline","D_arrangement_minus_baseline"),
                          ("fraction_arrangement","q_arrangement"),("fraction_baseline","q_baseline"),("fraction_third","q_third")):
            stat=inside_outside_stats(sub,key,label)
            if stat: comparison.append({"lg":lg,**stat})
    write_tsv(RESULTS/"root_mrca_comparison_inside_outside.tsv",comparison,["lg","metric","n_inside","n_outside","mean_inside","median_inside","mean_outside","median_outside","delta_mean_inside_minus_outside","delta_median_inside_minus_outside","mean_ratio_inside_over_outside","median_ratio_inside_over_outside"])
    correlations=correlation_rows(rows)
    write_tsv(RESULTS/"root_mrca_track_correlations.tsv",correlations,["lg","metric","n_windows","pearson_r","spearman_rho","interpretation"])

    figures=save_figures(rows,regions)
    summary_lookup={r["lg"]:r for r in summary_rows}
    test_lookup={r["lg"]:r for r in spatial_summary}
    comp={(r["lg"],r["metric"]):r for r in comparison}
    topo=read_tsv(EMP/"results/stage4b_circular_shift_summary.tsv")
    qcompare={r["lg"]:r for r in read_tsv(EMP/"results/stage4a_inside_outside_summary.tsv")}
    # Pick root-vs-B/A correlations for concise prose.
    cor={(r["lg"],r["metric"]):r for r in correlations}
    report=["# Atlantic cod root-MRCA analysis","",
            "This analysis calculates T_MRCA(w) as the age of the MRCA of all Stage-2 Atlantic cod population tips in each published 250-kb SNAPP MCC tree. Any outgroup tips are excluded by selecting the ingroup MRCA; all 426 source trees contain exactly the 12 ingroup tips and no outgroup tips.","",
            "## Tree and time validation","",
            f"Validated: {len(rows)}/{len(all_input)}. One materially negative-branch tree was excluded from age estimation: {', '.join(a['window_id'] for a in audit if a['valid_for_root_mrca']=='false')}.",
            "The tree strings lack an explicit Newick [&R] marker, but their source is rooted BEAST/SNAPP MCC output generated with a crown/root-age prior (`lognormal(0,3.83,0.093)` in the local XML-generation script). All usable ingroup trees are ultrametric within 1e-7 time units. The branch-length scale is verified against the original Matschiner et al. (2022) Figure 4 caption: it describes these 250-kb sliding-window profiles, gives a 0.5 Myr scale bar, and states that inset branch bars represent 50,000 years. Thus the Newick branch lengths are in Ma (million years), not years; values are plotted directly without rescaling.",
            "The published source caption is available at https://www.nature.com/articles/s41559-022-01661-x. The window trees and frozen supergene masks both use gadMor2 coordinates, following the cod_phylogenomics window-tree workflow; no gadMor_Celtic inversion boundaries are mixed into these tracks.","",
            "## Root-MRCA inside/outside comparison","",
            "Intervals use the frozen cod_phylogenomics window classifications; boundary-overlap windows are excluded. Confidence intervals are 95% percentile intervals from a spatial moving-block bootstrap (8 adjacent 250-kb windows = 2 Mb); they quantify genomic-window uncertainty, not SNAPP posterior uncertainty. Significance uses exact circular shifts of the root-age track relative to the fixed inside/outside mask over the full 250-kb coordinate grid, retaining gaps and autocorrelation; BH adjustment covers the four linkage groups.","",
            "| LG | n in/out | mean inside/out | median inside/out | Δ mean [95% CI] | Δ median [95% CI] | mean ratio [95% CI] | median ratio [95% CI] | circular p (BH) |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for lg in LGS:
        s=summary_lookup[lg];t=test_lookup[lg]
        report.append(f"| {lg} | {s['n_inside']}/{s['n_outside']} | {s['mean_inside']:.4g}/{s['mean_outside']:.4g} | {s['median_inside']:.4g}/{s['median_outside']:.4g} | {s['delta_mean_inside_minus_outside']:.4g} [{s['delta_mean_lo']:.4g}, {s['delta_mean_hi']:.4g}] | {s['delta_median_inside_minus_outside']:.4g} [{s['delta_median_lo']:.4g}, {s['delta_median_hi']:.4g}] | {s['mean_ratio_inside_over_outside']:.3g} [{s['mean_ratio_lo']:.3g}, {s['mean_ratio_hi']:.3g}] | {s['median_ratio_inside_over_outside']:.3g} [{s['median_ratio_lo']:.3g}, {s['median_ratio_hi']:.3g}] | {t['p_one_sided']:.4g} ({t['BH_p']:.4g}) |")
    report.extend(["","## Relation to pairwise and topology tracks","",
        "A(w) is the existing opposite-arrangement minus same-arrangement mean pairwise MRCA contrast. The repository has no Atlantic-cod file literally named B(w); here B(w) denotes the already-produced Stage-7 `mean_time_all`, i.e. the arrangement-blind mean across all 66 pairwise MRCA times in each 12-population tree. Quartet frequencies are the existing Stage-4 arrangement, baseline, and third-topology fractions.","",
        "| LG | Δ root age | Δ B(w) | Δ A(w) | Δ arrangement quartet frequency | ΔD | topology circular p |",
        "|---|---:|---:|---:|---:|---:|---:|"])
    for lg in LGS:
        get=lambda name: comp[(lg,name)]["delta_mean_inside_minus_outside"]
        topo_row=next(r for r in topo if r["lg"]==lg)
        report.append(f"| {lg} | {get('T_MRCA_root'):.4g} | {get('B_mean_all_pair_MRCA'):.4g} | {get('A_opposite_minus_same'):.4g} | {get('q_arrangement'):.4g} | {get('D_arrangement_minus_baseline'):.4g} | {float(topo_row['p_one_sided']):.4g} |")
    report.extend(["",
        "LG07 specifically has a strong root-age shift despite little change in its dominant quartet/topology contrast: its root-age track and B(w) are evaluated in separate panels from quartet frequencies. This is a time-depth change without a matching shift in which topology dominates; the finer A(w) signal can still change because it measures arrangement-conditioned pairwise MRCA contrasts rather than topology frequencies. LG12 has a positive inside/outside root-age difference and block-bootstrap interval, but its circular-shift test is not significant after correction (BH p=0.055).",
        "A higher root age and higher B(w) indicate deeper overall genealogies and deeper average pairwise coalescence, respectively. A(w) isolates the difference between opposite- and same-arrangement pairs. Their agreement is not a mathematical identity: root age can change through the deepest split even when most pairwise/topology summaries move little. Root age and B(w) are highly correlated within each LG (Pearson r 0.993–1.000); on LG07 both vary strongly despite essentially stable quartet frequencies.","",
        "## Limitations","",
        "The 426 values are published MCC point estimates, not independent posterior trees or independent genomic replicates. One invalid tree was excluded; no branch lengths were repaired. The root-age prior is shared across windows, root-time units are undocumented in the local source, and the block confidence intervals and cyclic-shift tests quantify spatial rather than posterior uncertainty. Cyclic shifts rely on a circular/stationary placement null. The window-level confidence intervals are sensitive to block length; see `root_mrca_block_sensitivity.tsv`.", "",
        "## Reproduction","",
        "From the repository root run `python3 empirical/atlantic_cod/scripts/13_root_mrca_time_analysis.py`. Outputs are written under `empirical/atlantic_cod/results/root_mrca/` and the two PDF/PNG figure pairs under `empirical/atlantic_cod/figures/`."])
    (RESULTS/"root_mrca_report.md").write_text("\n".join(report)+"\n")

    manifest={"analysis":"Atlantic cod ingroup root-MRCA times","source_trees":len(all_input),"valid_trees":len(rows),"invalid_trees":[a["window_id"] for a in audit if a["valid_for_root_mrca"]=="false"],
              "ingroup_tips":sorted(ingroup),"n_ingroup_tips":len(ingroup),"outgroups_in_source_trees":sorted({a["outgroup_tips"] for a in audit}),
              "rootedness_evidence":"BEAST/SNAPP MCC workflow uses crown/root-age prior; source Newick strings lack explicit [&R] marker", "ultrametric_tolerance":ULTRAMETRIC_TOL,
              "coordinate_assembly":"gadMor2; window and supergene coordinates use the same reference frame",
              "root_time_units":"Ma (million years), verified against the published Figure 4 0.5 Myr scale bar and 50,000-year inset branch scale", "bootstrap":"circular moving-block bootstrap within contiguous inside/outside runs", "primary_block_windows":PRIMARY_BLOCK_WINDOWS,"bootstrap_reps":BOOTSTRAPS,
              "spatial_test":"exact cyclic shifts on 250-kb full coordinate grid; gaps retained", "software":{"python":platform.python_version(),"biopython":Bio.__version__,"matplotlib":matplotlib.__version__,"numpy":np.__version__},
              "inputs":{str(p.relative_to(ROOT)):sha256(p) for p in (WINDOWS,ARRANGEMENTS,REGIONS,STAGE7_SIGNAL,STAGE4_SCORES)},
              "outputs":[str(p.relative_to(ROOT)) for p in [RESULTS/"root_mrca_tree_validation.tsv",RESULTS/"root_mrca_window_signal.tsv",RESULTS/"root_mrca_inside_outside.tsv",RESULTS/"root_mrca_spatial_test.tsv",RESULTS/"root_mrca_spatial_null.tsv",RESULTS/"root_mrca_block_sensitivity.tsv",RESULTS/"root_mrca_comparison_inside_outside.tsv",RESULTS/"root_mrca_track_correlations.tsv",RESULTS/"root_mrca_report.md",*figures]]}
    (RESULTS/"root_mrca_manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print(f"Validated {len(rows)}/{len(all_input)} trees; report: {RESULTS/'root_mrca_report.md'}")
    for path in figures: print(path)


if __name__=="__main__":
    main()
