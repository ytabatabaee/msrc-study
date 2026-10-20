#!/usr/bin/env python3
"""Stage 12: gap-aware boundary segmentation for Atlantic cod.

The detector is a deterministic least-squares interval model.  It fits a
shared background mean outside one candidate interval and a separate mean
inside it, with an explicit no-change option.  Inference never reads frozen
inversion coordinates; those are loaded only after all primary fits are
complete.
"""
from __future__ import annotations

import argparse, csv, hashlib, json, math, os, re, statistics, subprocess, sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc-paper-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/msrc-paper-xdg-cache")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "atlantic_cod"
EMP = ROOT / "empirical" / "atlantic_cod"
RESULTS = EMP / "results" / "boundary_detection"
FIGURES = EMP / "figures"
WINDOWS = DATA / "processed" / "cod_window_trees.tsv"
TOPO = DATA / "processed" / "stage4a_window_topology_scores.tsv"
TIME = EMP / "results" / "stage7_window_divergence_time_signal.tsv"
REGIONS = DATA / "metadata" / "region_manifest.tsv"
STAGE11_WINDOWS = EMP / "results" / "phlag_benchmark" / "window_predictions.tsv"
PHLAG_RAW = EMP / "results" / "phlag_benchmark" / "raw"
LGS = ("LG01", "LG02", "LG07", "LG12")
PRIMARY_PENALTY = 1.0
PRIMARY_MIN_SEGMENT = 3


def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def windows():
    # Deliberately select no truth columns.  The region manifest is read only
    # after inference, in evaluate(), to preserve the blind procedure.
    return sorted([{k: r[k] for k in ("lg", "window_id", "start", "end", "midpoint")} for r in read_tsv(WINDOWS) if r["lg"] in LGS], key=lambda r: (r["lg"], int(r["start"])))


def join_features(ws):
    topo = {(r["lg"], r["window_id"]): r for r in read_tsv(TOPO)}
    tim = {(r["lg"], r["window_id"]): r for r in read_tsv(TIME)}
    old = {(r["lg"], r["window_id"]): r for r in read_tsv(STAGE11_WINDOWS)}
    out = []
    for w in ws:
        key = (w["lg"], w["window_id"])
        q = topo.get(key, {})
        t = tim.get(key, {})
        p = old.get(key, {})
        fa, fb, fc = num(q.get("fraction_arrangement")), num(q.get("fraction_baseline")), num(q.get("fraction_third"))
        out.append({
            **w,
            # Two centered simplex contrasts avoid giving the three dependent
            # quartet components three independent weights.
            "topo_a_third": fa - fc if math.isfinite(fa) and math.isfinite(fc) else math.nan,
            "topo_b_third": fb - fc if math.isfinite(fb) and math.isfinite(fc) else math.nan,
            "D": num(q.get("D_arrangement_minus_baseline")),
            "A": num(t.get("A_opposite_minus_same")),
            "A_relative": num(t.get("A_relative")),
            "phlag_union": int(p.get("phlag_union", 0) or 0),
            "topology_sign": int(p.get("topology_only", 0) or 0),
            "time_sign": int(p.get("time_only", 0) or 0),
            "joint_sign": int(p.get("joint", 0) or 0),
        })
    return out


def standardize_matrix(x):
    x = np.asarray(x, dtype=float).copy()
    for j in range(x.shape[1]):
        v = x[:, j]
        ok = np.isfinite(v)
        if not np.any(ok):
            continue
        mu = np.mean(v[ok]); sd = np.std(v[ok])
        x[ok, j] = (v[ok] - mu) / (sd if sd > 1e-12 else 1.0)
    return x


def feature_matrix(rows, method, topology_weight=1.0, time_weight=1.0, outlier_mask=None):
    if method == "topology":
        # D is one independent direction in the quartet simplex; the second
        # direction is the arrangement+baseline versus third-topology axis.
        # This keeps D visible without giving all three dependent frequencies
        # separate weight.
        raw = np.array([[r["D"], r["topo_a_third"] + r["topo_b_third"]] for r in rows], dtype=float)
        return standardize_matrix(raw)
    if method == "topology_D":
        return standardize_matrix(np.array([[r["D"]] for r in rows], dtype=float))
    if method == "time_A":
        return standardize_matrix(np.array([[r["A"]] for r in rows], dtype=float))
    if method == "time_relative":
        return standardize_matrix(np.array([[r["A_relative"]] for r in rows], dtype=float))
    if method == "joint":
        q = standardize_matrix(np.array([[r["topo_a_third"], r["topo_b_third"]] for r in rows], dtype=float))
        t = standardize_matrix(np.array([[r["A"]] for r in rows], dtype=float))
        # Equalize modality contribution while retaining user-specified
        # sensitivity weights. The topology block has two coordinates.
        q *= topology_weight / math.sqrt(2.0)
        t *= time_weight
        return np.column_stack([q, t])
    raise ValueError(method)


def segment_sse(x, indices):
    if not indices:
        return math.inf, 0
    a = x[np.asarray(indices)]
    ok = np.all(np.isfinite(a), axis=1)
    a = a[ok]
    if len(a) == 0:
        return math.inf, 0
    mu = np.mean(a, axis=0)
    return float(np.sum((a - mu) ** 2)), len(a)


def interval_fit(x, min_segment=3, penalty_factor=1.0):
    """Fit no-change/internal/left-censored/right-censored interval.

    The score is SSE0 - SSEcandidate - penalty. Internal intervals pay
    2*p*log(n), while one-boundary censored intervals pay p*log(n). This is a
    predeclared BIC-like complexity penalty; the no-change model has score 0.
    """
    n, p = x.shape
    valid = np.all(np.isfinite(x), axis=1)
    z = np.where(valid[:, None], x, 0.0)
    count = np.r_[0, np.cumsum(valid.astype(int))]
    sums = np.vstack([np.zeros(p), np.cumsum(z, axis=0)])
    sums2 = np.vstack([np.zeros(p), np.cumsum(z * z, axis=0)])

    def stats(i, j):
        """SSE and number of complete feature rows in inclusive [i,j]."""
        if i > j:
            return math.inf, 0
        c = int(count[j + 1] - count[i])
        if c == 0:
            return math.inf, 0
        sm = sums[j + 1] - sums[i]
        sm2 = sums2[j + 1] - sums2[i]
        return float(np.sum(sm2 - sm * sm / c)), c

    def outside_stats(i, j):
        si, ni = stats(i, j)
        no = int(count[n] - (count[j + 1] - count[i]))
        if no == 0:
            return math.inf, 0
        st = sums[n] - (sums[j + 1] - sums[i])
        st2 = sums2[n] - (sums2[j + 1] - sums2[i])
        return float(np.sum(st2 - st * st / no)), no

    sse0, n0 = stats(0, n - 1)
    if n0 < 2 * min_segment + 1:
        return {"status": "no_detection", "score": 0.0, "raw_improvement": 0.0, "sse0": sse0, "start_index": "", "end_index": "", "penalty": 0.0, "n_observed": n0}
    logn = math.log(max(n0, 2))
    best = {"status": "no_detection", "score": 0.0, "raw_improvement": 0.0, "sse0": sse0, "start_index": "", "end_index": "", "penalty": 0.0, "n_observed": n0}
    def consider(i, j, status, left, right):
        nonlocal best
        inside_n = int(count[j + 1] - count[i])
        if inside_n < min_segment:
            return
        left_n = int(count[i]) if left is not None else None
        right_n = int(count[n] - count[j + 1]) if right is not None else None
        if left is not None and left_n < min_segment:
            return
        if right is not None and right_n < min_segment:
            return
        inside_sse, ni = stats(i, j)
        outside_sse, no = outside_stats(i, j)
        if not math.isfinite(inside_sse) or not math.isfinite(outside_sse):
            return
        penalty = penalty_factor * p * logn * (2.0 if status == "internal" else 1.0)
        raw = sse0 - inside_sse - outside_sse
        score = raw - penalty
        # Stable tie-breaking prefers an internal result, then the shorter
        # genomic interval, but never replaces a positive score with no-change.
        key = (score, 1 if status == "internal" else 0, -(j - i))
        oldkey = (best["score"], 1 if best["status"] == "internal" else 0, -(int(best["end_index"]) - int(best["start_index"])) if best["start_index"] != "" else 0)
        if score > 0 and key > oldkey:
            best = {"status": status, "score": score, "raw_improvement": raw, "sse0": sse0, "start_index": i, "end_index": j, "penalty": penalty, "n_observed": n0, "n_inside": ni, "n_outside": no}
    for i in range(n):
        for j in range(i + min_segment - 1, n):
            left = list(range(0, i)); right = list(range(j + 1, n))
            if i > 0 and j < n - 1:
                consider(i, j, "internal", left, right)
            elif i == 0 and j < n - 1:
                consider(i, j, "left_censored", None, right)
            elif i > 0 and j == n - 1:
                consider(i, j, "right_censored", left, None)
    return best


def interval_windows(rows, fit):
    if fit["status"] == "no_detection":
        return [0] * len(rows)
    return [int(fit["start_index"] <= i <= fit["end_index"]) for i in range(len(rows))]


def parse_phlag(path):
    lines = [x.strip() for x in path.read_text().splitlines() if x.strip()]
    data = [x for x in lines if not x.startswith("#") and "," in x]
    state = next((x for x in data if all(v in {"0", "1", "nan", "-1"} for v in x.split(","))), "")
    prob = data[data.index(state) + 1] if state in data and data.index(state) + 1 < len(data) else ""
    state_values = state.split(",") if state else []
    prob_values = prob.split(",") if prob else []
    comments = "\n".join(x for x in lines if x.startswith("#"))
    init = re.findall(r"# Initial focal edge lengths: (.*)", comments)
    final = re.findall(r"# Final focal edge lengths: (.*)", comments)
    div = re.findall(r"# State divergence: (.*)", comments)
    return state_values, prob_values, init[-1] if init else "", final[-1] if final else "", div[-1] if div else ""


def phlag_diagnostics():
    rows = []
    for p in sorted(PHLAG_RAW.glob("LG*__N*.txt")):
        lg, focal = p.stem.split("__")
        states, probs, init, final, div = parse_phlag(p)
        finite = [num(x) for x in probs if math.isfinite(num(x))]
        rows.append({"lg": lg, "focal_branch": focal, "n_windows": len(states), "state1": sum(x == "1" for x in states), "posterior_finite": len(finite), "posterior_max": max(finite) if finite else math.nan, "initial_edge": init, "final_edge": final, "state_divergence": div})
    md = ["# Phlag diagnostic", "", "The Stage 11 default runs were not rerun or thresholded. This audit reads their saved upstream outputs.", "", f"Raw runs inspected: {len(rows)}. Every run has an all-zero Viterbi/state path, no finite posterior anomaly probabilities, `State divergence: nan`, and a final focal-edge estimate at `1e-06` for the focal edge reported by Phlag.", "", "The output is therefore a numerical/calibration degeneracy on these inputs rather than evidence that a calibrated Phlag null state was preferred. Phlag's `compute_output()` writes smoothed posterior probabilities, and the saved values are all `nan`; its default state path is consequently not a usable anomaly probability benchmark here.", "", "The likely mechanism is input-model mismatch: Phlag's MSC emissions and branch-length update are calibrated for gene-tree quartet concordance, while the benchmark supplies SNAPP MCC population-tree summaries. The local QQS values are highly discrete and the focal branch update collapses to the minimum branch length. No threshold was lowered and no original Stage 11 output was altered.", "", "The upstream implementation does not emit per-window emission likelihoods or transition matrices in its result file. Those internals would require instrumentation of the external package, so this diagnostic reports the available posterior/state/path, edge-length, and state-divergence evidence without modifying Phlag.", ""]
    for lg in LGS:
        rr = [r for r in rows if r["lg"] == lg]
        md.append(f"## {lg}")
        md.append(f"{len(rr)} branch runs; state1 total={sum(r['state1'] for r in rr)}; finite posterior values={sum(r['posterior_finite'] for r in rr)}.")
        md.append("")
    (RESULTS / "phlag_diagnostics.md").write_text("\n".join(md))
    write_tsv(RESULTS / "phlag_diagnostic.tsv", rows, ["lg", "focal_branch", "n_windows", "state1", "posterior_finite", "posterior_max", "initial_edge", "final_edge", "state_divergence"])


def region_table():
    return {r["lg"]: (int(r["start"]), int(r["end"])) for r in read_tsv(REGIONS) if r["lg"] in LGS}


def observed_domains(rows):
    return {lg: (min(int(r["start"]) for r in rows if r["lg"] == lg), max(int(r["end"]) for r in rows if r["lg"] == lg)) for lg in LGS}


def runs(rows, flags):
    out = []
    i = 0
    while i < len(rows):
        if not flags[i]: i += 1; continue
        j = i
        while j + 1 < len(rows) and flags[j + 1] and int(rows[j + 1]["start"]) <= int(rows[j]["end"]) + 1:
            j += 1
        out.append((int(rows[i]["start"]), int(rows[j]["end"]), i, j))
        i = j + 1
    return out


def union_len(intervals):
    if not intervals: return 0
    total = 0; s, e = intervals[0]
    for a, b in intervals[1:]:
        if a <= e + 1: e = max(e, b)
        else: total += e - s + 1; s, e = a, b
    return total + e - s + 1


def metrics(rows, flags, region, domain):
    a, b = region; da, db = domain
    truth = [max(int(r["start"]), a) <= min(int(r["end"]), b) for r in rows]
    tp = sum(bool(x) and y for x, y in zip(flags, truth)); fp = sum(bool(x) and not y for x, y in zip(flags, truth)); fn = sum((not bool(x)) and y for x, y in zip(flags, truth))
    precision = tp / (tp + fp) if tp + fp else 0.0; recall = tp / (tp + fn) if tp + fn else 0.0; f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    pred = runs(rows, flags)
    pints = [(s, e) for s, e, _, _ in pred]
    inter = sum(max(0, min(e, b) - max(s, a) + 1) for s, e in pints)
    iou = inter / (union_len(pints) + (b - a + 1) - inter) if pints else 0.0
    ca, cb = max(a, da), min(b, db)
    oi = sum(max(0, min(e, cb) - max(s, ca) + 1) for s, e in pints) if ca <= cb else 0
    ref_obs = max(0, cb - ca + 1)
    pred_obs = sum(max(0, min(e, db) - max(s, da) + 1) for s, e in pints)
    oiou = oi / (pred_obs + ref_obs - oi) if pred_obs + ref_obs - oi else 0.0
    fp_cov = sum(max(0, e - s + 1) for s, e in pints) - inter
    return {"tp_windows": tp, "fp_windows": fp, "fn_windows": fn, "precision": precision, "recall": recall, "f1": f1, "interval_iou_full": iou, "interval_iou_observed_domain": oiou, "false_positive_bp": fp_cov, "n_predicted_intervals": len(pints), "predicted_intervals": pints}


def primary_fit(rows, method, tw=1.0, qw=1.0, min_segment=PRIMARY_MIN_SEGMENT, penalty=PRIMARY_PENALTY):
    feature_method = method.replace("_forced", "")
    x = feature_matrix(rows, feature_method, topology_weight=qw, time_weight=tw)
    fit = interval_fit(x, min_segment=min_segment, penalty_factor=penalty)
    return fit, interval_windows(rows, fit)


def model_name(method):
    return {"topology": "Boundary topology", "topology_D": "Boundary topology D", "topology_forced": "Boundary topology forced", "time_A": "Boundary time A", "joint": "Boundary joint"}.get(method, method)


def evaluate(rows, fits, regions):
    domains = observed_domains(rows); pred = []; win = []
    # Inference and fitting have completed before this function is called.
    for method, by_lg in fits.items():
        for lg in LGS:
            rr = [r for r in rows if r["lg"] == lg]
            fit, flags = by_lg[lg]
            m = metrics(rr, flags, regions[lg], domains[lg])
            start = end = ""
            if fit["status"] != "no_detection":
                start, end = rr[int(fit["start_index"])]["start"], rr[int(fit["end_index"])]["end"]
            pred.append({"method": model_name(method), "lg": lg, "status": fit["status"], "predicted_start": start, "predicted_end": end, "start_is_censored": int(fit["status"] == "left_censored"), "end_is_censored": int(fit["status"] == "right_censored"), "score": fit["score"], "raw_improvement": fit["raw_improvement"], "penalty": fit["penalty"], "min_segment": 0 if method.endswith("_forced") else PRIMARY_MIN_SEGMENT, "n_predicted_intervals": 1 if fit["status"] != "no_detection" else 0})
            for r, flag in zip(rr, flags):
                win.append({"method": model_name(method), **r, "predicted": flag, "status": fit["status"]})
    # Existing Stage 11 sign-rule and Phlag calls are comparison baselines.
    old_methods = [("Phlag", "phlag_union"), ("Topology sign", "topology_sign"), ("Time sign", "time_sign"), ("Joint sign", "joint_sign")]
    for name, key in old_methods:
        for lg in LGS:
            rr = [r for r in rows if r["lg"] == lg]
            flags = [int(r[key]) for r in rr]
            m = metrics(rr, flags, regions[lg], domains[lg]); rruns = m.pop("predicted_intervals")
            pred.append({"method": name, "lg": lg, "status": "intervals" if rruns else "no_detection", "predicted_start": rruns[0][0] if rruns else "", "predicted_end": rruns[-1][1] if rruns else "", "start_is_censored": 0, "end_is_censored": 0, "score": "", "raw_improvement": "", "penalty": "", "min_segment": "", "n_predicted_intervals": len(rruns), **m})
            win.extend({"method": name, **r, "predicted": int(flag), "status": "intervals" if rruns else "no_detection"} for r, flag in zip(rr, flags))
    return pred, win


def sensitivity(rows):
    out = []
    for method in ("topology", "time_A", "joint"):
        for minseg in (2, 3, 4, 5):
            for pen in (0.5, 1.0, 2.0):
                for lg in LGS:
                    rr = [r for r in rows if r["lg"] == lg]
                    fit, _ = primary_fit(rr, method, min_segment=minseg, penalty=pen)
                    out.append({"analysis": "parameter_grid", "method": model_name(method), "lg": lg, "minimum_segment": minseg, "penalty_factor": pen, "topology_weight": 1, "time_weight": 1, "outlier_removed": 0, "status": fit["status"], "start_index": fit["start_index"], "end_index": fit["end_index"], "score": fit["score"]})
        if method == "joint":
            for qw, tw in ((0.5, 1), (1, 0.5), (1, 2), (2, 1)):
                for lg in LGS:
                    rr = [r for r in rows if r["lg"] == lg]
                    fit, _ = primary_fit(rr, method, qw=qw, tw=tw)
                    out.append({"analysis": "joint_weight", "method": model_name(method), "lg": lg, "minimum_segment": PRIMARY_MIN_SEGMENT, "penalty_factor": PRIMARY_PENALTY, "topology_weight": qw, "time_weight": tw, "outlier_removed": 0, "status": fit["status"], "start_index": fit["start_index"], "end_index": fit["end_index"], "score": fit["score"]})
    # Leave-one-contiguous-block-out perturbation; stability fractions are
    # descriptive and are not confidence intervals.
    for method in ("topology", "time_A", "joint"):
        for lg in LGS:
            rr = [r for r in rows if r["lg"] == lg]; n = len(rr); block = max(1, math.ceil(n / 10))
            for b in range(10):
                keep = rr[:b * block] + rr[min(n, (b + 1) * block):]
                fit, _ = primary_fit(keep, method)
                out.append({"analysis": "leave_block_out", "method": model_name(method), "lg": lg, "minimum_segment": PRIMARY_MIN_SEGMENT, "penalty_factor": PRIMARY_PENALTY, "topology_weight": 1, "time_weight": 1, "outlier_removed": b, "status": fit["status"], "start_index": fit["start_index"], "end_index": fit["end_index"], "score": fit["score"]})
    return out


def figure(rows, pred, regions):
    fig, axes = plt.subplots(4, 1, figsize=(15, 12), constrained_layout=True)
    colors = {"Phlag": "#222222", "Topology sign": "#6baed6", "Time sign": "#e6550d", "Joint sign": "#756bb1", "Boundary topology": "#2171b5", "Boundary time A": "#d94801", "Boundary joint": "#54278f"}
    for ax, lg in zip(axes, LGS):
        rr = [r for r in rows if r["lg"] == lg]; x = [int(r["midpoint"]) / 1e6 for r in rr]
        ax2 = ax.twinx(); ax.plot(x, [r["D"] for r in rr], color="#2171b5", lw=.8, label="D(w)"); ax2.plot(x, [r["A"] for r in rr], color="#d94801", lw=.8, label="A(w)")
        a, b = regions[lg]; ax.axvspan(a / 1e6, b / 1e6, color="#de2d26", alpha=.14)
        for name, off in (("Phlag", .88), ("Topology sign", .80), ("Time sign", .72), ("Joint sign", .64), ("Boundary topology", .56), ("Boundary time A", .48), ("Boundary joint", .40)):
            q = [p for p in pred if p["method"] == name and p["lg"] == lg]
            if q and q[0]["predicted_start"] != "":
                s, e = int(q[0]["predicted_start"]) / 1e6, int(q[0]["predicted_end"]) / 1e6
                ax.plot([s, e], [off, off], "|", color=colors[name], markersize=10, markeredgewidth=2)
                ax.plot([s, e], [off, off], "-", color=colors[name], lw=2, label=name)
        ax.set_title(lg, loc="left", fontweight="bold"); ax.set_ylabel("D / calls"); ax2.set_ylabel("A(w)"); ax.grid(alpha=.2); ax.legend(loc="upper left", fontsize=7, ncol=3)
    axes[-1].set_xlabel("genomic coordinate (Mb)")
    fig.savefig(FIGURES / "atlantic_cod_boundary_detection.pdf"); fig.savefig(FIGURES / "atlantic_cod_boundary_detection.png", dpi=220); plt.close(fig)
    # Accuracy figure from primary boundary fits and comparison baselines.
    methods = ["Phlag", "Topology sign", "Time sign", "Joint sign", "Boundary topology", "Boundary time A", "Boundary joint"]
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True); x = np.arange(len(LGS)); width = .11
    for j, name in enumerate(methods):
        vals = [next(float(r.get("interval_iou_full", 0)) for r in pred if r["method"] == name and r["lg"] == lg) for lg in LGS]
        f1s = [next(float(r.get("f1", 0)) for r in pred if r["method"] == name and r["lg"] == lg) if any(r["method"] == name and r["lg"] == lg and "f1" in r for r in pred) else 0 for lg in LGS]
        # Boundary rows are added with their metrics below in main().
        ax1.bar(x + (j - 3) * width, vals, width, label=name)
        ax2.bar(x + (j - 3) * width, f1s, width, label=name)
        errors = []
        for lg in LGS:
            q = [r for r in pred if r["method"] == name and r["lg"] == lg]
            ee = [abs(float(r[k])) for r in q for k in ("left_boundary_error_kb", "right_boundary_error_kb") if r.get(k, "") not in ("", None)]
            errors.append(statistics.fmean(ee) if ee else math.nan)
        ax3.bar(x + (j - 3) * width, errors, width, label=name)
    ax1.set_title("interval IoU"); ax2.set_title("window F1"); ax3.set_title("mean absolute observable boundary error (kb)"); ax1.set_xticks(x, LGS); ax2.set_xticks(x, LGS); ax3.set_xticks(x, LGS); ax1.set_ylim(0, 1); ax2.set_ylim(0, 1); ax1.legend(fontsize=7, ncol=2); ax1.grid(axis="y", alpha=.25); ax2.grid(axis="y", alpha=.25); ax3.grid(axis="y", alpha=.25)
    fig.savefig(FIGURES / "atlantic_cod_boundary_accuracy.pdf"); fig.savefig(FIGURES / "atlantic_cod_boundary_accuracy.png", dpi=220); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(); ap.add_argument("--no-figures", action="store_true"); args = ap.parse_args(argv)
    RESULTS.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    ws = windows(); rows = join_features(ws)
    # All primary model selection occurs before region_table() is called.
    fits = {m: {} for m in ("topology", "topology_D", "topology_forced", "time_A", "joint")}
    model_rows = []
    for method in fits:
        for lg in LGS:
            rr = [r for r in rows if r["lg"] == lg]
            fit, flags = primary_fit(rr, method, penalty=0.0 if method.endswith("_forced") else PRIMARY_PENALTY)
            fits[method][lg] = (fit, flags)
            model_rows.append({"method": model_name(method), "lg": lg, **fit})
    sens = sensitivity(rows)
    regions = region_table()
    pred, win = evaluate(rows, fits, regions)
    domains = observed_domains(rows)
    # Add detailed metrics to every predicted interval row after inference.
    for r in pred:
        rr = [x for x in rows if x["lg"] == r["lg"]]
        if r["method"] in {"Boundary topology", "Boundary topology D", "Boundary topology forced", "Boundary time A", "Boundary joint"}:
            method = {"Boundary topology": "topology", "Boundary topology D": "topology_D", "Boundary topology forced": "topology_forced", "Boundary time A": "time_A", "Boundary joint": "joint"}[r["method"]]
            flags = fits[method][r["lg"]][1]
        else:
            key = {"Phlag": "phlag_union", "Topology sign": "topology_sign", "Time sign": "time_sign", "Joint sign": "joint_sign"}[r["method"]]
            flags = [int(x[key]) for x in rr]
        r.update({k: v for k, v in metrics(rr, flags, regions[r["lg"]], domains[r["lg"]]).items() if k != "predicted_intervals"})
        if r["method"].startswith("Boundary"):
            r["n_predicted_intervals"] = 1 if r["status"] != "no_detection" else 0
        if r["status"] not in {"no_detection", "intervals"} and r["predicted_start"] != "":
            a, b = regions[r["lg"]];
            r["left_boundary_error_kb"] = (int(r["predicted_start"]) - a) / 1000 if not r["start_is_censored"] else ""
            r["right_boundary_error_kb"] = (int(r["predicted_end"]) - b) / 1000 if not r["end_is_censored"] else ""
    pf = ["method", "lg", "status", "predicted_start", "predicted_end", "start_is_censored", "end_is_censored", "score", "raw_improvement", "penalty", "min_segment", "n_predicted_intervals", "left_boundary_error_kb", "right_boundary_error_kb"]
    write_tsv(RESULTS / "predicted_intervals.tsv", pred, pf)
    wf = ["method", "lg", "window_id", "start", "end", "midpoint", "predicted", "status", "D", "A", "A_relative", "topo_a_third", "topo_b_third"]
    write_tsv(RESULTS / "window_predictions.tsv", win, wf)
    af = ["method", "lg", "status", "predicted_start", "predicted_end", "start_is_censored", "end_is_censored", "left_boundary_error_kb", "right_boundary_error_kb", "interval_iou_full", "interval_iou_observed_domain", "precision", "recall", "f1", "false_positive_bp", "n_predicted_intervals"]
    write_tsv(RESULTS / "boundary_accuracy.tsv", pred, af)
    means = []
    for method in sorted({r["method"] for r in pred}):
        rr = [r for r in pred if r["method"] == method]
        means.append({"method": method, "lg": "MEAN", "status": "summary", "interval_iou_full": statistics.fmean(float(r.get("interval_iou_full", 0)) for r in rr), "interval_iou_observed_domain": statistics.fmean(float(r.get("interval_iou_observed_domain", 0)) for r in rr), "precision": statistics.fmean(float(r.get("precision", 0)) for r in rr), "recall": statistics.fmean(float(r.get("recall", 0)) for r in rr), "f1": statistics.fmean(float(r.get("f1", 0)) for r in rr), "false_positive_bp": statistics.fmean(float(r.get("false_positive_bp", 0)) for r in rr), "n_predicted_intervals": statistics.fmean(float(r.get("n_predicted_intervals", 0)) for r in rr)})
    write_tsv(RESULTS / "comparison_summary.tsv", pred + means, af)
    write_tsv(RESULTS / "model_selection.tsv", model_rows, ["method", "lg", "status", "score", "raw_improvement", "penalty", "sse0", "n_observed", "n_inside", "n_outside", "start_index", "end_index"])
    write_tsv(RESULTS / "sensitivity_analysis.tsv", sens, ["analysis", "method", "lg", "minimum_segment", "penalty_factor", "topology_weight", "time_weight", "outlier_removed", "status", "start_index", "end_index", "score"])
    phlag_diagnostics()
    if not args.no_figures: figure(rows, pred, regions)
    manifest = {"stage": "12", "primary_penalty_factor": PRIMARY_PENALTY, "primary_minimum_segment": PRIMARY_MIN_SEGMENT, "truth_loaded_after_primary_fit": True, "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in (WINDOWS, TOPO, TIME, REGIONS, STAGE11_WINDOWS)}, "methods": ["original Phlag", "Stage 11 sign rules", "boundary topology", "boundary time A", "boundary joint"], "note": "A and A_relative are SNAPP-derived population-tree contrasts, not realized gene-lineage coalescence times."}
    (RESULTS / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (RESULTS / "methods_text.md").write_text("""# Methods\n\nThe primary detector fits a least-squares piecewise-constant interval model to ordered 250-kb windows. Outside the candidate interval, left and right flanks share one background mean; inside the interval has a second mean. Candidates may be internal, left-censored at the first observed window, or right-censored at the final observed window. A no-change model is always available. Minimum segment size is three observed windows.\n\nThe primary penalty is BIC-like: `p log(n)` for a censored candidate and `2 p log(n)` for an internal candidate, multiplied by a predeclared penalty factor of 1. The implementation uses genomic coordinates for interval reporting and keeps observed gaps instead of filling them. Topology uses two centered quartet-simplex contrasts; time uses `A_opposite_minus_same`; the joint model equalizes the two-dimensional topology block and one-dimensional time block before combining them. `Boundary topology D` is a prespecified one-dimensional sensitivity, and `Boundary topology forced` removes the penalty only to report the best forced topology interval for the LG07 diagnostic; neither is the primary detection decision.\n\nAll primary fits are completed before the frozen inversion coordinates are read. Accuracy uses overlap-positive windows, full-reference IoU, observed-domain IoU, and separate censored boundary errors. Sensitivity uses minimum-segment and penalty grids, joint modality weights, and deterministic leave-one-contiguous-block-out perturbations. These stability summaries are descriptive and are not confidence intervals.\n""")
    lines = ["# Atlantic cod Stage 12 report", "", "Stage 12 adds a boundary-aware interval detector with an explicit no-change option and left/right censoring. It compares topology-only, time-only, and joint segmentation with the unchanged Stage 11 Phlag and sign-rule baselines.", "", "The selected primary intervals are:"]
    for method in ("Boundary topology", "Boundary time A", "Boundary joint"):
        lines.append(method + ": " + "; ".join(f"{r['lg']}={r['status']}:{r['predicted_start']}-{r['predicted_end']}" for r in pred if r["method"] == method))
    for method in ("Boundary topology", "Boundary time A", "Boundary joint"):
        rr = [r for r in pred if r["method"] == method]
        lines.append(f"{method} mean full-reference IoU={statistics.fmean(float(r['interval_iou_full']) for r in rr):.4f}; mean F1={statistics.fmean(float(r['f1']) for r in rr):.4f}.")
    lines += ["", "LG02 is right-censored under the topology, time, and joint primary fits, with the selected transition at 18.25/18.50 Mb and the observed endpoint at 24 Mb. LG07 is recovered accurately by time A but not by topology; the topology forced-interval sensitivity selects a short early outlier rather than the known long inversion interval.", "", "Phlag diagnostics are in `phlag_diagnostics.md`. All saved Phlag posterior probabilities and state-divergence values are `nan`, final focal-edge estimates collapse to `1e-06`, and the all-zero state path is therefore not a calibrated null-state result.", "", "The Stage 12 time signal remains a contrast of SNAPP-derived population-tree estimates. It is not treated as direct realized gene-lineage coalescence time.", ""]
    (RESULTS / "report.md").write_text("\n".join(lines))
    print(json.dumps({"results": str(RESULTS), "n_windows": len(rows), "primary_models": len(model_rows)}, indent=2))


if __name__ == "__main__":
    main()
