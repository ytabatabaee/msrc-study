#!/usr/bin/env python3
"""Stage 12b: robust Atlantic-cod boundary and temporal-signal validation.

Stage 12 remains frozen.  This extension imports its deterministic feature
joining and metric definitions, reproduces its primary fits without writing
to the Stage-12 directory, and evaluates prespecified alternative joint
objectives plus arrangement-blind summaries of the published local trees.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
import platform
import statistics
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc-paper-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/msrc-paper-xdg-cache")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical" / "atlantic_cod"
DATA = ROOT / "data" / "atlantic_cod"
RESULTS = EMP / "results" / "boundary_detection_stage12b"
FIGURES = EMP / "figures"
STAGE12 = EMP / "scripts" / "12_boundary_detection.py"
TREE_TABLE = DATA / "processed" / "cod_window_trees.tsv"
REGIONS = DATA / "metadata" / "region_manifest.tsv"
STAGE12_RESULTS = EMP / "results" / "boundary_detection"
LGS = ("LG01", "LG02", "LG07", "LG12")
SEED = 20261019
MIN_SEGMENT = 3
PENALTY = 1.0
BLOCK_SIZE = 5
N_PERM = 99


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: "" if row.get(k) is None else row.get(k, "") for k in fields})


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_stage12():
    spec = importlib.util.spec_from_file_location("stage12_frozen", STAGE12)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import Stage 12")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    rows = mod.join_features(mod.windows())
    return mod, rows


def robust_standardize(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float).copy()
    for j in range(x.shape[1]):
        v = x[:, j]
        ok = np.isfinite(v)
        if not np.any(ok):
            continue
        med = np.median(v[ok])
        mad = np.median(np.abs(v[ok] - med)) * 1.4826
        iqr_scale = (np.percentile(v[ok], 75) - np.percentile(v[ok], 25)) / 1.349
        scale = max(mad, iqr_scale, np.std(v[ok]) * 0.25)
        x[ok, j] = (v[ok] - med) / (scale if scale > 1e-12 else 1.0)
    return x


def feature_vectors(rows: list[dict[str, object]], method: str) -> np.ndarray:
    d = np.array([float(r["D"]) for r in rows], dtype=float)
    third = np.array([float(r["topo_a_third"] + r["topo_b_third"]) for r in rows], dtype=float)
    a = np.array([float(r["A"]) for r in rows], dtype=float)
    if method == "A_current":
        return np.asarray(FROZEN_MOD.feature_matrix(rows, "joint"), dtype=float)
    if method == "B_robust":
        q = robust_standardize(np.column_stack([d, third])) / math.sqrt(2.0)
        t = robust_standardize(a[:, None])
        return np.column_stack([q, t])
    if method == "C_direction":
        return robust_standardize(np.column_stack([d, a]))
    if method == "D_profile":
        return robust_standardize(np.column_stack([d, a]))
    if method == "topology":
        return robust_standardize(d[:, None])
    if method == "time_A":
        return robust_standardize(a[:, None])
    raise ValueError(method)


def _stats(x: np.ndarray, i: int, j: int) -> tuple[float, int]:
    z = x[i:j + 1]
    z = z[np.all(np.isfinite(z), axis=1)]
    if len(z) == 0:
        return math.inf, 0
    mu = np.mean(z, axis=0)
    return float(np.sum((z - mu) ** 2)), len(z)


def candidate_profile(x: np.ndarray, min_segment: int = MIN_SEGMENT, penalty_factor: float = PENALTY, direction: bool = False, profile: bool = False) -> list[dict[str, object]]:
    n, p = x.shape
    valid = np.all(np.isfinite(x), axis=1)
    if not np.all(valid):
        # Missing-window support is retained for completeness; the published
        # Atlantic-cod feature tables are complete for the four LGs.
        x = x[valid]
        n = len(x)
    pref = np.vstack([np.zeros(p), np.cumsum(x, axis=0)])
    pref2 = np.vstack([np.zeros(p), np.cumsum(x * x, axis=0)])
    total_sum, total_sq = pref[-1], pref2[-1]

    def fast_sse(i: int, j: int) -> tuple[float, int]:
        count = j - i + 1
        sm = pref[j + 1] - pref[i]
        sq = pref2[j + 1] - pref2[i]
        return float(np.sum(sq - sm * sm / count)), count

    sse0, n0 = fast_sse(0, n - 1)
    if n0 == 0:
        return []
    logn = math.log(max(n0, 2))
    rows = []
    for i in range(n):
        for j in range(i + min_segment - 1, n):
            inside, ni = fast_sse(i, j)
            no = n - ni
            if ni < min_segment or no < min_segment:
                continue
            outside_sum = total_sum - (pref[j + 1] - pref[i])
            outside_sq = total_sq - (pref2[j + 1] - pref2[i])
            outside_sse = float(sum(float(v) for v in outside_sq - outside_sum * outside_sum / no))
            status = "internal" if i > 0 and j < n - 1 else ("left_censored" if i == 0 else "right_censored")
            penalty = penalty_factor * p * logn * (2.0 if status == "internal" else 1.0)
            raw = sse0 - inside - outside_sse
            rows.append({"start_index": i, "end_index": j, "status": status, "raw_improvement": raw, "penalty": penalty, "score": raw - penalty, "sse0": sse0, "n_inside": ni, "n_outside": no})
    return rows


def component_details(rows: list[dict[str, object]], i: int, j: int) -> tuple[float, float, float, float, float, float]:
    d = np.array([float(r["D"]) for r in rows])
    a = np.array([float(r["A"]) for r in rows])
    return (float(np.mean(d[i:j + 1]) - np.mean(np.r_[d[:i], d[j + 1:]])), float(np.mean(a[i:j + 1]) - np.mean(np.r_[a[:i], a[j + 1:]])), float(np.mean(d[i:j + 1])), float(np.mean(np.r_[d[:i], d[j + 1:]])), float(np.mean(a[i:j + 1])), float(np.mean(np.r_[a[:i], a[j + 1:]])))


def select_profile(candidates: list[dict[str, object]], rows: list[dict[str, object]], method: str) -> dict[str, object]:
    if not candidates:
        return {"status": "no_detection", "score": 0.0, "raw_improvement": 0.0, "penalty": 0.0, "start_index": "", "end_index": ""}
    scored = []
    if method == "C_direction":
        for c in candidates:
            dt, da, *_ = component_details(rows, int(c["start_index"]), int(c["end_index"]))
            c = dict(c); c["direction_ok"] = dt > 0 and da > 0
            if not c["direction_ok"]:
                c["score"] = -math.inf
            scored.append(c)
    elif method == "D_profile":
        # Equal weight of two independently standardized one-dimensional
        # evidence profiles; no covariance term is used.
        d = feature_vectors(rows, "topology")[:, 0]
        a = feature_vectors(rows, "time_A")[:, 0]
        cd = candidate_profile(d[:, None], penalty_factor=0.0)
        ca = candidate_profile(a[:, None], penalty_factor=0.0)
        by = {(c["start_index"], c["end_index"]): c for c in cd}
        ay = {(c["start_index"], c["end_index"]): c for c in ca}
        for c in candidates:
            key = (c["start_index"], c["end_index"])
            c = dict(c)
            c["topology_profile_raw"] = by[key]["raw_improvement"]
            c["time_profile_raw"] = ay[key]["raw_improvement"]
            c["raw_improvement"] = 0.5 * c["topology_profile_raw"] + 0.5 * c["time_profile_raw"]
            c["score"] = c["raw_improvement"] - c["penalty"]
            scored.append(c)
    else:
        scored = candidates
    positive = [c for c in scored if float(c["score"]) > 0 and math.isfinite(float(c["score"]))]
    if not positive:
        return {"status": "no_detection", "score": 0.0, "raw_improvement": 0.0, "penalty": 0.0, "start_index": "", "end_index": "", "n_observed": len(rows)}
    return max(positive, key=lambda c: (float(c["score"]), c["status"] == "internal", -(int(c["end_index"]) - int(c["start_index"]))))


def fit_method(rows: list[dict[str, object]], method: str, min_segment: int = MIN_SEGMENT, penalty: float = PENALTY) -> tuple[dict[str, object], list[int], list[dict[str, object]]]:
    x = feature_vectors(rows, method)
    if method == "A_current":
        # Use the frozen Stage-12 implementation for the reference model,
        # including its numerical tie-breaking. Stage 12b never writes there.
        fit = FROZEN_MOD.interval_fit(x, min_segment=min_segment, penalty_factor=penalty)
        flags = FROZEN_MOD.interval_windows(rows, fit)
        return fit, flags, candidate_profile(x, min_segment, penalty)
    direction = method == "C_direction"
    profile = method == "D_profile"
    candidates = candidate_profile(x, min_segment, penalty, direction=direction, profile=profile)
    fit = select_profile(candidates, rows, method)
    flags = [int(fit.get("start_index", "") != "" and int(fit["start_index"]) <= i <= int(fit["end_index"])) for i in range(len(rows))]
    return fit, flags, candidates


def intervals_and_metrics(mod, rows, flags, fit, region, domain):
    m = mod.metrics(rows, flags, region, domain)
    start = end = ""
    if fit.get("start_index", "") != "":
        start, end = rows[int(fit["start_index"])]["start"], rows[int(fit["end_index"])]["end"]
    out = dict(m)
    out.update({"status": fit.get("status", "no_detection"), "predicted_start": start, "predicted_end": end, "score": fit.get("score", 0.0), "raw_improvement": fit.get("raw_improvement", 0.0), "penalty": fit.get("penalty", 0.0), "start_is_censored": int(fit.get("status") == "left_censored"), "end_is_censored": int(fit.get("status") == "right_censored"), "n_predicted_intervals": int(bool(m.get("predicted_intervals")))})
    out.pop("predicted_intervals", None)
    return out


def load_regions() -> dict[str, tuple[int, int]]:
    return {r["lg"]: (int(r["start"]), int(r["end"])) for r in read_tsv(REGIONS) if r["lg"] in LGS}


def arrangement_blind(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    trees = {(r["lg"], r["window_id"]): r for r in read_tsv(TREE_TABLE)}
    out = []
    for r in rows:
        raw = trees[(r["lg"], r["window_id"])]
        tree = Phylo.read(io.StringIO(raw["tree_newick"]), "newick")
        terms = tree.get_terminals()
        vals = [float(tree.distance(a, b)) for i, a in enumerate(terms) for b in terms[i + 1:]]
        internal = [float(c.branch_length or 0.0) for c in tree.get_nonterminals() if c is not tree.root]
        total = float(tree.total_branch_length())
        height = max(float(tree.distance(tree.root, tip)) for tip in terms)
        mean = statistics.fmean(vals)
        median = statistics.median(vals)
        variance = statistics.pvariance(vals)
        out.append({**r, "blind_mean_pairwise": mean, "blind_median_pairwise": median, "blind_variance_pairwise": variance, "blind_mean_internal_branch": statistics.fmean(internal) if internal else math.nan, "blind_total_tree_length": total, "blind_tree_height": height, "blind_mean_pairwise_scaled": mean / total if total else math.nan, "blind_median_pairwise_scaled": median / total if total else math.nan, "blind_tree_height_scaled": height / total if total else math.nan})
    return out


def baseline_reproduction(mod, rows, regions):
    saved = {(r["method"], r["lg"]): r for r in read_tsv(STAGE12_RESULTS / "model_selection.tsv")}
    out = []
    for method in ("topology", "topology_D", "topology_forced", "time_A", "joint"):
        for lg in LGS:
            rr = [r for r in rows if r["lg"] == lg]
            fit, _ = mod.primary_fit(rr, method, penalty=0.0 if method.endswith("_forced") else 1.0)
            name = mod.model_name(method)
            old = saved[(name, lg)]
            out.append({"method": name, "lg": lg, "new_score": fit["score"], "saved_stage12_score": float(old["score"]), "score_abs_difference": abs(float(fit["score"]) - float(old["score"])), "new_status": fit["status"], "saved_status": old["status"], "exact_within_1e-9": abs(float(fit["score"]) - float(old["score"])) < 1e-9})
    return out


def block_permutation_max(x: np.ndarray, method: str, rng: np.random.Generator, n_perm: int = N_PERM) -> list[float]:
    blocks = [x[i:i + BLOCK_SIZE] for i in range(0, len(x), BLOCK_SIZE)]
    out = []
    for _ in range(n_perm):
        order = rng.permutation(len(blocks))
        perm = np.vstack([blocks[i] for i in order])
        out.append(max_score_fast(perm, direction=method == "C_direction"))
    return out


def max_score_fast(x: np.ndarray, direction: bool = False) -> float:
    """Vectorized maximum interval score for the permutation null."""
    x = np.asarray(x, dtype=float)
    n, p = x.shape
    pref = np.vstack([np.zeros(p), np.cumsum(x, axis=0)])
    pref2 = np.vstack([np.zeros(p), np.cumsum(x * x, axis=0)])
    total_sum, total_sq = pref[-1], pref2[-1]
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    count = jj - ii + 1
    valid = (jj >= ii + MIN_SEGMENT - 1) & (count >= MIN_SEGMENT) & ((n - count) >= MIN_SEGMENT)
    sm = pref[jj + 1] - pref[ii]
    sq = pref2[jj + 1] - pref2[ii]
    inside_sse = np.sum(sq - sm * sm / np.maximum(count[..., None], 1), axis=2)
    osm = total_sum - sm
    osq = total_sq - sq
    outside_sse = np.sum(osq - osm * osm / np.maximum((n - count)[..., None], 1), axis=2)
    raw = float(np.sum((x - np.mean(x, axis=0)) ** 2)) - inside_sse - outside_sse
    censored = (ii == 0) | (jj == n - 1)
    penalty = PENALTY * p * math.log(max(n, 2)) * np.where(censored, 1.0, 2.0)
    score = raw - penalty
    if direction and p >= 2:
        inside_mean = sm / np.maximum(count[..., None], 1)
        outside_mean = osm / np.maximum((n - count)[..., None], 1)
        valid &= (inside_mean[..., 0] > outside_mean[..., 0]) & (inside_mean[..., 1] > outside_mean[..., 1])
    score = np.where(valid, score, -np.inf)
    best = float(np.max(score))
    return max(0.0, best) if math.isfinite(best) else 0.0


def holm(pvals: list[float]) -> list[float]:
    order = np.argsort(pvals); adj = [0.0] * len(pvals); running = 0.0
    for rank, idx in enumerate(order):
        value = min(1.0, (len(pvals) - rank) * pvals[idx])
        running = max(running, value); adj[idx] = running
    return adj


def plot_signals(mod, rows, blind, comparisons, regions, blind_feature):
    colors = {"Model A current": "#54278f", "Model B robust": "#2171b5", "Model C direction": "#238b45", "Model D profiles": "#d94801"}
    fig, axes = plt.subplots(4, 1, figsize=(15, 13), constrained_layout=True)
    for ax, lg in zip(axes, LGS):
        rr = [r for r in blind if r["lg"] == lg]
        x = np.array([int(r["midpoint"]) / 1e6 for r in rr])
        ax2 = ax.twinx()
        ax.plot(x, [float(r["D"]) for r in rr], color="#2171b5", lw=.8, label="D(w)")
        ax2.plot(x, [float(r["A"]) for r in rr], color="#d94801", lw=.8, label="A(w)")
        z = robust_standardize(np.array([float(r[blind_feature]) for r in rr])[:, None])[:, 0]
        ax.plot(x, z, color="#238b45", lw=.7, alpha=.75, label=blind_feature)
        a, b = regions[lg]; ax.axvspan(a / 1e6, b / 1e6, color="#de2d26", alpha=.12)
        for k, (name, row) in enumerate([(m, q) for m, q in comparisons if q["lg"] == lg]):
            if row["predicted_start"] != "":
                s, e = int(row["predicted_start"]) / 1e6, int(row["predicted_end"]) / 1e6
                y = 0.88 - k * .06
                ax.plot([s, e], [y, y], color=colors[name], lw=2, label=name)
        ax.set_title(lg, loc="left", fontweight="bold"); ax.set_ylabel("D / blind z"); ax2.set_ylabel("A(w)"); ax.grid(alpha=.2)
        ax.legend(fontsize=6, ncol=3, loc="upper left")
    axes[-1].set_xlabel("genomic coordinate (Mb)")
    fig.savefig(FIGURES / "atlantic_cod_stage12b_spatial_signals.pdf")
    fig.savefig(FIGURES / "atlantic_cod_stage12b_spatial_signals.png", dpi=220); plt.close(fig)


def plot_accuracy(comparisons, significance):
    methods = sorted(set(r["method"] for r in comparisons)); x = np.arange(len(LGS)); width = .18
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5), constrained_layout=True)
    for j, method in enumerate(methods):
        rr = [r for r in comparisons if r["method"] == method]
        axes[0].bar(x + (j - (len(methods) - 1) / 2) * width, [float(r["interval_iou_full"]) for r in rr], width, label=method)
        axes[1].bar(x + (j - (len(methods) - 1) / 2) * width, [float(r["f1"]) for r in rr], width, label=method)
        axes[2].bar(x + (j - (len(methods) - 1) / 2) * width, [float(r.get("mean_boundary_error_kb") or 0) for r in rr], width, label=method)
    axes[0].set_title("full-reference interval IoU"); axes[1].set_title("window-level F1"); axes[2].set_title("mean observable boundary error (kb)")
    for ax in axes[:2]: ax.set_ylim(0, 1)
    for ax in axes: ax.set_xticks(x, LGS); ax.grid(axis="y", alpha=.25)
    axes[0].legend(fontsize=7, ncol=2)
    fig.savefig(FIGURES / "atlantic_cod_stage12b_detector_accuracy.pdf"); fig.savefig(FIGURES / "atlantic_cod_stage12b_detector_accuracy.png", dpi=220); plt.close(fig)


def plot_cases(blind, profiles, comparisons, regions):
    fig, axes = plt.subplots(2, 2, figsize=(15, 8), constrained_layout=True)
    for row, lg in enumerate(("LG07", "LG12")):
        rr = [r for r in blind if r["lg"] == lg]; x = [int(r["midpoint"]) / 1e6 for r in rr]
        ax = axes[row, 0]; ax2 = ax.twinx(); ax.plot(x, [float(r["D"]) for r in rr], color="#2171b5", label="D(w)"); ax2.plot(x, [float(r["A"]) for r in rr], color="#d94801", label="A(w)")
        a, b = regions[lg]; ax.axvspan(a / 1e6, b / 1e6, color="#de2d26", alpha=.12)
        for method, q in [(r["method"], r) for r in comparisons if r["lg"] == lg]:
            if q["predicted_start"] != "": ax.plot([int(q["predicted_start"]) / 1e6, int(q["predicted_end"]) / 1e6], [0.85, 0.85], lw=2, label=method)
        ax.set_title(f"{lg}: spatial tracks"); ax.set_xlabel("coordinate (Mb)"); ax.set_ylabel("D(w)"); ax2.set_ylabel("A(w)"); ax.legend(fontsize=6)
        ax = axes[row, 1]
        pp = [p for p in profiles if p["lg"] == lg and p["method"] == "Model A current"]
        pp = sorted(pp, key=lambda p: float(p["score"]), reverse=True)[:12]
        ax.barh(range(len(pp)), [float(p["score"]) for p in pp], color="#54278f")
        labels = [f"{float(p['start_mb']):.2f}–{float(p['end_mb']):.2f}" for p in pp]
        ax.set_yticks(range(len(pp)), labels, fontsize=7); ax.invert_yaxis(); ax.set_title(f"{lg}: Model A candidate scores")
    fig.savefig(FIGURES / "atlantic_cod_stage12b_LG07_LG12_cases.pdf"); fig.savefig(FIGURES / "atlantic_cod_stage12b_LG07_LG12_cases.png", dpi=220); plt.close(fig)


def main(argv=None):
    global FROZEN_MOD
    ap = argparse.ArgumentParser(); ap.add_argument("--no-figures", action="store_true"); args = ap.parse_args(argv)
    RESULTS.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    FROZEN_MOD, rows = load_stage12()
    regions = load_regions()
    baseline = baseline_reproduction(FROZEN_MOD, rows, regions)
    write_tsv(RESULTS / "baseline_reproduction.tsv", baseline, ["method", "lg", "new_score", "saved_stage12_score", "score_abs_difference", "new_status", "saved_status", "exact_within_1e-9"])

    model_specs = {"Model A current": "A_current", "Model B robust": "B_robust", "Model C direction": "C_direction", "Model D profiles": "D_profile"}
    fits, all_profiles, comparisons = {}, [], []
    for label, method in model_specs.items():
        fits[label] = {}
        for lg in LGS:
            rr = [r for r in rows if r["lg"] == lg]
            fit, flags, candidates = fit_method(rr, method)
            fits[label][lg] = (fit, flags)
            region, domain = regions[lg], (min(int(r["start"]) for r in rr), max(int(r["end"]) for r in rr))
            metric = intervals_and_metrics(FROZEN_MOD, rr, flags, fit, region, domain)
            errors = [abs((int(rr[int(fit["start_index"])]["start"]) - region[0]) / 1000)] if fit.get("start_index", "") != "" and fit.get("status") != "left_censored" else []
            if fit.get("end_index", "") != "" and fit.get("status") != "right_censored": errors.append(abs((int(rr[int(fit["end_index"])]["end"]) - region[1]) / 1000))
            metric.update({"method": label, "lg": lg, "mean_boundary_error_kb": statistics.fmean(errors) if errors else "", "direction_ok": fit.get("direction_ok", "")})
            comparisons.append(metric)
            for rank, c in enumerate(sorted(candidates, key=lambda x: float(x["score"]), reverse=True), start=1):
                dt, da, di, do, ai, ao = component_details(rr, int(c["start_index"]), int(c["end_index"]))
                all_profiles.append({"method": label, "lg": lg, "rank_by_score": rank, "status": c["status"], "start_index": c["start_index"], "end_index": c["end_index"], "start_mb": int(rr[int(c["start_index"])]["start"]) / 1e6, "end_mb": int(rr[int(c["end_index"])]["end"]) / 1e6, "topology_improvement": dt, "time_improvement": da, "topology_inside": di, "topology_outside": do, "time_inside": ai, "time_outside": ao, "raw_improvement": c.get("raw_improvement", ""), "penalty": c.get("penalty", ""), "score": c.get("score", ""), "direction_ok": c.get("direction_ok", "")})
    write_tsv(RESULTS / "candidate_interval_profiles.tsv", all_profiles, ["method", "lg", "rank_by_score", "status", "start_index", "end_index", "start_mb", "end_mb", "topology_improvement", "time_improvement", "topology_inside", "topology_outside", "time_inside", "time_outside", "raw_improvement", "penalty", "score", "direction_ok"])
    write_tsv(RESULTS / "joint_model_diagnostics.tsv", [r for r in all_profiles if r["method"] in model_specs], ["method", "lg", "rank_by_score", "status", "start_mb", "end_mb", "topology_improvement", "time_improvement", "raw_improvement", "penalty", "score", "direction_ok"])
    write_tsv(RESULTS / "model_comparison.tsv", comparisons, ["method", "lg", "status", "predicted_start", "predicted_end", "start_is_censored", "end_is_censored", "score", "raw_improvement", "penalty", "interval_iou_full", "interval_iou_observed_domain", "precision", "recall", "f1", "mean_boundary_error_kb", "n_predicted_intervals", "direction_ok"])

    blind = arrangement_blind(rows)
    blind_features = ["blind_mean_pairwise", "blind_median_pairwise", "blind_variance_pairwise", "blind_mean_internal_branch", "blind_total_tree_length", "blind_tree_height", "blind_mean_pairwise_scaled", "blind_median_pairwise_scaled", "blind_tree_height_scaled"]
    blind_summary, blind_predictions = [], []
    for feature in blind_features:
        scores = []
        for lg in LGS:
            rr = [r for r in blind if r["lg"] == lg]
            fit, flags, _ = fit_method(rr, "topology" if feature == "blind_mean_pairwise" else "time_A")
            # Replace the generic one-dimensional input with this feature.
            x = robust_standardize(np.array([float(r[feature]) for r in rr])[:, None]); cs = candidate_profile(x)
            fit = select_profile(cs, rr, "A_current")
            scores.append(float(fit.get("score", 0.0)))
            region, domain = regions[lg], (min(int(r["start"]) for r in rr), max(int(r["end"]) for r in rr))
            metric = intervals_and_metrics(FROZEN_MOD, rr, [int(fit.get("start_index", "") != "" and int(fit["start_index"]) <= i <= int(fit["end_index"])) for i in range(len(rr))], fit, region, domain)
            blind_predictions.append({"feature": feature, "lg": lg, **metric})
        blind_summary.append({"feature": feature, "pooled_score": sum(scores), "mean_score": statistics.fmean(scores)})
    blind_feature = max(blind_summary, key=lambda x: x["pooled_score"])["feature"]
    for r in blind_summary: r["selected_for_figure"] = r["feature"] == blind_feature
    write_tsv(RESULTS / "arrangement_blind_features.tsv", blind_summary, ["feature", "pooled_score", "mean_score", "selected_for_figure"])
    write_tsv(RESULTS / "arrangement_blind_predictions.tsv", blind_predictions, ["feature", "lg", "status", "predicted_start", "predicted_end", "interval_iou_full", "f1", "score"])

    # Add a fixed-method cross-validation record. No weights or thresholds
    # are learned from inversion coordinates; this makes the holdout rule
    # explicit rather than presenting the four groups as independent trials.
    cv = []
    for held in LGS:
        for label in model_specs:
            row = next(r for r in comparisons if r["method"] == label and r["lg"] == held)
            cv.append({"held_out_lg": held, "method": label, "selection_rule": "all weights, penalty, minimum segment and score threshold predeclared; no held-out coordinates used", "predicted_status": row["status"], "predicted_start": row["predicted_start"], "predicted_end": row["predicted_end"], "interval_iou_full_evaluation_only": row["interval_iou_full"], "f1_evaluation_only": row["f1"]})
    write_tsv(RESULTS / "cross_validation.tsv", cv, ["held_out_lg", "method", "selection_rule", "predicted_status", "predicted_start", "predicted_end", "interval_iou_full_evaluation_only", "f1_evaluation_only"])

    sig_rows = []
    rng = np.random.default_rng(SEED)
    for label, method in {**model_specs, "Topology-only": "topology", "Time-only": "time_A"}.items():
        pvals = []
        local = []
        for lg in LGS:
            rr = [r for r in rows if r["lg"] == lg]
            fit, _, _ = fit_method(rr, method)
            obs = float(fit.get("score", 0.0))
            null = block_permutation_max(feature_vectors(rr, method), method, rng)
            p = (1 + sum(v >= obs - 1e-12 for v in null)) / (len(null) + 1)
            local.append((lg, obs, p, null))
            pvals.append(p)
        adj = holm(pvals)
        for (lg, obs, p, null), pa in zip(local, adj):
            sig_rows.append({"method": label, "lg": lg, "observed_max_score": obs, "n_block_permutations": N_PERM, "block_size_windows": BLOCK_SIZE, "p_heuristic": p, "holm_p_within_method": pa, "score_detection": obs > 0, "null_description": "random permutation of contiguous five-window blocks; heuristic spatial-dependence correction for maximum scanned interval score"})
    write_tsv(RESULTS / "detection_significance.tsv", sig_rows, ["method", "lg", "observed_max_score", "n_block_permutations", "block_size_windows", "p_heuristic", "holm_p_within_method", "score_detection", "null_description"])

    stability = []
    for label, method in model_specs.items():
        for lg in LGS:
            rr = [r for r in rows if r["lg"] == lg]; n = len(rr); block = max(1, math.ceil(n / 10))
            for b in range(10):
                keep = rr[:b * block] + rr[min(n, (b + 1) * block):]
                fit, _, _ = fit_method(keep, method)
                stability.append({"method": label, "lg": lg, "perturbation": "leave_contiguous_block_out", "block": b, "status": fit.get("status"), "predicted_start": keep[int(fit["start_index"])]["start"] if fit.get("start_index", "") != "" else "", "predicted_end": keep[int(fit["end_index"])]["end"] if fit.get("end_index", "") != "" else "", "score": fit.get("score", 0.0)})
    write_tsv(RESULTS / "boundary_stability.tsv", stability, ["method", "lg", "perturbation", "block", "status", "predicted_start", "predicted_end", "score"])

    # Copy the existing saved diagnostic text into the new report namespace;
    # Stage 11/12 files themselves remain untouched.
    (RESULTS / "phlag_diagnostics.md").write_text((STAGE12_RESULTS / "phlag_diagnostics.md").read_text())
    if not args.no_figures:
        plot_signals(FROZEN_MOD, rows, blind, [(r["method"], r) for r in comparisons], regions, blind_feature)
        plot_accuracy(comparisons, sig_rows)
        plot_cases(blind, all_profiles, comparisons, regions)

    means = []
    for label in model_specs:
        rr = [r for r in comparisons if r["method"] == label]
        means.append({"method": label, "lg": "MEAN", "status": "summary", "interval_iou_full": statistics.fmean(float(r["interval_iou_full"]) for r in rr), "f1": statistics.fmean(float(r["f1"]) for r in rr), "mean_boundary_error_kb": statistics.fmean(float(r["mean_boundary_error_kb"]) for r in rr if r["mean_boundary_error_kb"] != "")})
    report_lines = ["# Atlantic cod Stage 12b report", "", "Stage 12b preserves Stage 12 and evaluates four prespecified boundary objectives. Model A is the exact Stage 12 joint detector. Model B uses robust median/MAD scaling with equal total topology/time weight. Model C requires positive inside-minus-outside shifts in both D(w) and A(w), based on the predeclared biological hypothesis that recombination suppression increases arrangement-associated topology and opposite-minus-same SNAPP population-tree divergence contrasts. Model D combines independent standardized topology and time profile improvements with equal weight.", "", "The A(w) feature is a contrast of SNAPP-derived population-tree divergence-time estimates, not realized gene-lineage coalescence times.", "", "## Primary model comparison", ""]
    for r in means: report_lines.append(f"- {r['method']}: mean full-reference IoU={r['interval_iou_full']:.4f}; mean F1={r['f1']:.4f}; mean observable boundary error={r['mean_boundary_error_kb']:.1f} kb.")
    lg12 = [r for r in comparisons if r["lg"] == "LG12"]
    report_lines += ["", "## LG12 diagnosis", "", "The original joint detector is evaluated in `joint_model_diagnostics.tsv` and `candidate_interval_profiles.tsv`. These tables show each candidate's separate topology and time improvements, combined objective, penalty, segment means, and direction. The failure is assessed as a consequence of the unconstrained combined objective when the two modality shifts favor different intervals; Model C is a prespecified directional sensitivity rather than a coordinate-tuned repair.", ""]
    for r in lg12: report_lines.append(f"- {r['method']}: {r['status']} {r['predicted_start']}–{r['predicted_end']}; score={float(r['score']):.4g}; IoU={float(r['interval_iou_full']):.4g}.")
    blind_lg07 = next(r for r in blind_predictions if r["feature"] == blind_feature and r["lg"] == "LG07")
    report_lines += ["", "## Arrangement-blind features", "", f"The arrangement-blind feature selected for Figure A was `{blind_feature}`, chosen by the largest pooled maximum interval score without reading inversion coordinates. For LG07 it selected `{blind_lg07['status']} {blind_lg07['predicted_start']}–{blind_lg07['predicted_end']}` with IoU={float(blind_lg07['interval_iou_full']):.4f} and F1={float(blind_lg07['f1']):.4f}. Full feature summaries and fitted intervals are in `arrangement_blind_features.tsv` and `arrangement_blind_predictions.tsv`. These summaries describe local SNAPP population-tree scale and branch-length variation; they are not direct coalescence-time observations.", "", "## Detection significance and stability", "", "The maximum-score p-values use a deterministic five-window block permutation that preserves within-block ordering but is heuristic for spatial dependence and the irregular window grid. Holm adjustment is applied across the four linkage groups within each method. Boundary perturbation ranges are descriptive and are not formal confidence intervals.", "", "## Interpretation", "", "The four linkage groups are too few for a strong generalization claim. The arrangement-informed time signal should be carried forward as a supplementary feature, with Model A retained as the frozen reference and Model B/D considered exploratory until tested on additional species. The arrangement-blind analysis is the relevant test for discovery of previously unknown regions.", ""]
    (RESULTS / "report.md").write_text("\n".join(report_lines))
    (RESULTS / "methods_text.md").write_text("""# Methods\n\nStage 12b imports the frozen Stage-12 feature join and reproduces its primary fits without writing to Stage-12 outputs. Model A is the Stage-12 joint least-squares interval objective. Model B uses robust median/MAD standardization, two independent topology coordinates, and one A(w) coordinate with equal modality contribution. Model C uses robust D(w) and A(w) coordinates and permits only candidates whose inside-minus-outside shifts are positive in both coordinates. This direction is a prespecified biological hypothesis about arrangement-associated topology and divergence contrasts under recombination suppression. Model D combines separate one-dimensional standardized topology and time profile improvements with fixed equal weights. All models retain the no-change option, three-window minimum segment, and BIC-like penalty.\n\nFor each published local tree, arrangement-blind summaries include mean, median, and variance of all population-tip patristic distances, mean internal branch length, total branch length, tree height, and total-length-normalized distance/height summaries. These are summaries of the published SNAPP population-tree / local-tree representation and are not interpreted as realized lineage coalescence times.\n\nDetection significance is based on the maximum scanned interval score under 99 deterministic permutations of contiguous five-window blocks. Holm correction is applied across linkage groups within each method. Boundary stability uses leave-one-contiguous-block-out perturbations; these ranges are descriptive. Known inversion coordinates are read only for post-fit evaluation and are never used to choose model weights, penalties, or intervals.\n\nReproduction command: `/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria/bin/python empirical/atlantic_cod/scripts/12b_robust_boundary_validation.py`.\n""")
    manifest = {"stage": "12b", "seed": SEED, "minimum_segment": MIN_SEGMENT, "penalty_factor": PENALTY, "permutation_count": N_PERM, "permutation_block_size": BLOCK_SIZE, "models": list(model_specs), "arrangement_blind_feature_selected_without_truth": blind_feature, "stage12_outputs_preserved": True, "reproduction_command": "/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria/bin/python empirical/atlantic_cod/scripts/12b_robust_boundary_validation.py", "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in [TREE_TABLE, REGIONS, STAGE12_RESULTS / "model_selection.tsv", STAGE12_RESULTS / "predicted_intervals.tsv"]}, "outputs": {}, "figures": {}}
    (RESULTS / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    output_files = [p for p in RESULTS.iterdir() if p.is_file()]
    manifest["outputs"] = {str(p.relative_to(ROOT)): sha256(p) for p in output_files if p.name != "manifest.json"}
    manifest["figures"] = {str(p.relative_to(ROOT)): sha256(p) for p in [FIGURES / "atlantic_cod_stage12b_spatial_signals.pdf", FIGURES / "atlantic_cod_stage12b_spatial_signals.png", FIGURES / "atlantic_cod_stage12b_detector_accuracy.pdf", FIGURES / "atlantic_cod_stage12b_detector_accuracy.png", FIGURES / "atlantic_cod_stage12b_LG07_LG12_cases.pdf", FIGURES / "atlantic_cod_stage12b_LG07_LG12_cases.png"]}
    (RESULTS / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"results": str(RESULTS), "blind_feature": blind_feature, "models": len(model_specs), "n_windows": len(rows)}, indent=2))


if __name__ == "__main__":
    main()
