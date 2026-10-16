#!/usr/bin/env python3
"""Stage 10: fit a two-background MSRC to MRCA times and predict quartets.

The primary fit is a pairwise composite likelihood.  It uses the frozen
collinear topology and outside-window SNAPP times to project a new baseline
time tree, then fits one shared same-background coalescence rate and one
effective switching rate per linkage group.  Quartet observations are used
only after each cross-validation fit has been completed.

The quartet prediction is deliberately restricted to the exact symmetric 2:2
special case documented in the repository manuscript bundle.  For that case
the available theory gives A_m(t) = P(arrangement) - P(third), and symmetry
implies P(third) = P(baseline).  The code does not manufacture predictions
for other background configurations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import random
import sys
import tempfile
import textwrap
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "msrc-paper-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "msrc-paper-xdg-cache"))

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "atlantic_cod"
EMP = ROOT / "empirical" / "atlantic_cod"
RESULTS = EMP / "results"
FIGURES = EMP / "figures"

BASELINE_TOPOLOGY = DATA / "processed" / "baseline_population_tree_12taxa.nwk"
MRCA_PATH = DATA / "processed" / "stage7_pairwise_mrca_times.tsv"
QUARTET_PATH = DATA / "processed" / "stage4a_local_quartets.tsv"
SCORE_PATH = DATA / "processed" / "stage4a_window_topology_scores.tsv"
ARRANGEMENT_PATH = DATA / "processed" / "population_arrangements.tsv"
REGION_PATH = DATA / "metadata" / "region_manifest.tsv"
S8_PATH = EMP / "scripts" / "08_summary_tree_time_projection.py"

LGS = ("LG01", "LG02", "LG07", "LG12")
SEED = 20261016

BASELINE_TREE_OUT = RESULTS / "stage10_baseline_time_tree.nwk"
PROJECTION_OUT = RESULTS / "stage10_baseline_time_projection.tsv"
PARAM_OUT = RESULTS / "stage10_msrc_parameter_estimates.tsv"
PROFILE_OUT = RESULTS / "stage10_msrc_parameter_profiles.tsv"
CV_ASSIGN_OUT = RESULTS / "stage10_cv_window_assignment.tsv"
CV_PRED_OUT = RESULTS / "stage10_cv_predictions.tsv"
PRED_SUMMARY_OUT = RESULTS / "stage10_prediction_summary.tsv"
INCOMPAT_OUT = RESULTS / "stage10_incompatible_mrca_audit.tsv"
SIM_OUT = RESULTS / "stage10_simulation_validation.tsv"
METHODS_OUT = RESULTS / "stage10_methods_text.md"
RESULTS_TEXT_OUT = RESULTS / "stage10_results_text.md"
REPORT_OUT = RESULTS / "stage10_report.md"
MANIFEST_OUT = RESULTS / "stage10_manifest.json"
FIG_PDF = FIGURES / "atlantic_cod_stage10_msrc_parameter_prediction.pdf"
FIG_PNG = FIGURES / "atlantic_cod_stage10_msrc_parameter_prediction.png"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return "NA" if not math.isfinite(value) else f"{value:.10g}"
    return str(value)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_stage8():
    spec = importlib.util.spec_from_file_location("stage8_projection_for_stage10", S8_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {S8_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


S8 = load_stage8()


def pair_key(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a, b)))


def parse_pair_rows() -> list[dict[str, object]]:
    rows = []
    for raw in read_tsv(MRCA_PATH):
        rows.append({
            **raw,
            "start": int(raw["start"]),
            "end": int(raw["end"]),
            "mrca_time": float(raw["mrca_time"]),
            "pair": pair_key(raw["population1"], raw["population2"]),
            "state_sign": 1.0 if raw["arrangement1"] == raw["arrangement2"] else -1.0,
        })
    return rows


def outside_matrix(rows: list[dict[str, object]]) -> dict[tuple[str, str], float]:
    values: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        if row["region_class"] == "outside":
            values[row["pair"]].append(float(row["mrca_time"]))
    if len(values) != 66 or any(not values[pair] for pair in values):
        raise ValueError(f"Expected 66 outside population pairs, found {len(values)}")
    return {pair: float(np.median(times)) for pair, times in values.items()}


def fit_baseline_time_tree(rows: list[dict[str, object]]):
    matrix = outside_matrix(rows)
    nodes = S8.parse_topology(BASELINE_TOPOLOGY)
    edge, root_audit = S8.choose_root("baseline_collinear", nodes, matrix)
    root = S8.rooted_from_edge(nodes, edge)
    fit = S8.fit_node_ages(root, matrix, "baseline_collinear", "outside", S8.root_edge_label(nodes, edge))
    S8.validate_projected(fit.rooted_root, sorted(S8.taxa_of_graph(nodes)))
    RESULTS.mkdir(parents=True, exist_ok=True)
    S8.write_projected_tree(BASELINE_TREE_OUT, fit)

    projection_rows = []
    for pair in sorted(matrix):
        node = S8.mrca_node(fit.rooted_root, pair[0], pair[1])
        projected = float(fit.ages[node.id])
        projection_rows.append({
            "population1": pair[0],
            "population2": pair[1],
            "n_outside": sum(1 for r in rows if r["region_class"] == "outside" and r["pair"] == pair),
            "outside_median_mrca": matrix[pair],
            "projected_divergence_time": projected,
            "residual": matrix[pair] - projected,
            "abs_residual": abs(matrix[pair] - projected),
            "mrca_node_id": node.id,
        })
    projection_rows.extend({
        "population1": "ROOT_AUDIT",
        "population2": row["candidate_root_edge"],
        "n_outside": "",
        "outside_median_mrca": "",
        "projected_divergence_time": "",
        "residual": row["loss_outside"],
        "abs_residual": "selected" if row["selected"] else "",
        "mrca_node_id": "",
    } for row in root_audit)
    write_tsv(PROJECTION_OUT, projection_rows, [
        "population1", "population2", "n_outside", "outside_median_mrca",
        "projected_divergence_time", "residual", "abs_residual", "mrca_node_id",
    ])
    divergence = {
        (row["population1"], row["population2"]): float(row["projected_divergence_time"])
        for row in projection_rows if row["population1"] != "ROOT_AUDIT"
    }
    return fit, divergence, matrix, root_audit


def make_observation_rows(rows: list[dict[str, object]], divergence: dict[tuple[str, str], float], clipped: bool = False) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    observations = []
    audit = []
    for row in rows:
        if row["region_class"] == "boundary":
            continue
        d = divergence[row["pair"]]
        t_obs = float(row["mrca_time"])
        incompatible = t_obs < d - 1e-12
        audit.append({
            "lg": row["lg"], "window_id": row["window_id"], "start": row["start"], "end": row["end"],
            "population1": row["population1"], "population2": row["population2"],
            "observed_mrca_time": t_obs, "baseline_divergence_time": d,
            "incompatibility_amount": max(0.0, d - t_obs), "incompatible": incompatible,
            "primary_status": "excluded_incompatible" if incompatible else "included",
        })
        if incompatible and not clipped:
            continue
        observations.append({
            **row,
            "d": d,
            "T": max(t_obs, d + 1e-9) if clipped else t_obs,
            "incompatible": incompatible,
        })
    return observations, audit


def observations_to_arrays(rows: list[dict[str, object]], lgs: tuple[str, ...]):
    if not rows:
        raise ValueError("No observations available for fitting")
    return {
        "T": np.asarray([float(r["T"]) for r in rows], dtype=float),
        "d": np.asarray([float(r["d"]) for r in rows], dtype=float),
        "sign": np.asarray([float(r["state_sign"]) for r in rows], dtype=float),
        "lg_index": np.asarray([lgs.index(str(r["lg"])) for r in rows], dtype=int),
    }


def matrix_exponential_density(T: np.ndarray, d: np.ndarray, sign: np.ndarray, lam: float, m: float) -> np.ndarray:
    """Evaluate alpha exp(Q (T-d)) c using the closed 2x2 matrix exponential.

    The expression is the analytic matrix exponential for the exact generator
    Q in the methods text; it avoids numerical time discretization and keeps
    optimization fast over tens of thousands of observations.
    """
    u = np.maximum(T - d, 0.0)
    alpha_s = 0.5 * (1.0 + sign * np.exp(np.clip(-4.0 * m * d, -745.0, 0.0)))
    q00 = -(lam + 2.0 * m)
    q01 = 2.0 * m
    q10 = 2.0 * m
    q11 = -2.0 * m
    half_trace = -(lam + 4.0 * m) / 2.0
    delta = math.sqrt((lam / 2.0) ** 2 + 4.0 * m * m)
    z = delta * u
    exp_part = np.exp(np.clip(half_trace * u, -745.0, 700.0))
    cosh_z = np.cosh(np.clip(z, -350.0, 350.0))
    sinh_over_delta = np.sinh(np.clip(z, -350.0, 350.0)) / delta
    e00 = exp_part * (cosh_z + sinh_over_delta * (q00 - half_trace))
    e10 = exp_part * (sinh_over_delta * q10)
    first_column = alpha_s * e00 + (1.0 - alpha_s) * e10
    return np.maximum(lam * first_column, np.finfo(float).tiny)


def log_likelihood(log_params: np.ndarray, arrays: dict[str, np.ndarray], lgs: tuple[str, ...], model: str) -> float:
    values = np.exp(np.clip(log_params, -20.0, 20.0))
    lam = float(values[0])
    if model == "shared_lambda_lg_m":
        m_by_lg = values[1:]
        m = m_by_lg[arrays["lg_index"]]
        density = np.empty_like(arrays["T"])
        for index, value in enumerate(m_by_lg):
            mask = arrays["lg_index"] == index
            density[mask] = matrix_exponential_density(arrays["T"][mask], arrays["d"][mask], arrays["sign"][mask], lam, float(value))
    elif model == "lg_lambda_lg_m":
        lam_by_lg = values[: len(lgs)]
        m_by_lg = values[len(lgs):]
        density = np.empty_like(arrays["T"])
        for index in range(len(lgs)):
            mask = arrays["lg_index"] == index
            density[mask] = matrix_exponential_density(arrays["T"][mask], arrays["d"][mask], arrays["sign"][mask], float(lam_by_lg[index]), float(m_by_lg[index]))
    elif model == "one_lambda_one_m":
        density = matrix_exponential_density(arrays["T"], arrays["d"], arrays["sign"], lam, float(values[1]))
    else:
        raise ValueError(model)
    return float(np.log(density).sum())


def nelder_mead(objective, start: np.ndarray, max_iter: int = 700) -> tuple[np.ndarray, float]:
    n = len(start)
    simplex = [np.asarray(start, dtype=float)]
    for i in range(n):
        point = np.asarray(start, dtype=float).copy()
        point[i] += 0.55
        simplex.append(point)
    values = [objective(point) for point in simplex]
    for _ in range(max_iter):
        order = np.argsort(values)
        simplex = [simplex[i] for i in order]
        values = [values[i] for i in order]
        if np.max(np.abs(np.asarray(simplex[1:]) - simplex[0])) < 1e-6 and abs(values[-1] - values[0]) < 1e-5:
            break
        centroid = np.mean(simplex[:-1], axis=0)
        reflected = centroid + (centroid - simplex[-1])
        reflected_value = objective(reflected)
        if values[0] <= reflected_value < values[-2]:
            simplex[-1], values[-1] = reflected, reflected_value
            continue
        if reflected_value < values[0]:
            expanded = centroid + 2.0 * (reflected - centroid)
            expanded_value = objective(expanded)
            simplex[-1], values[-1] = (expanded, expanded_value) if expanded_value < reflected_value else (reflected, reflected_value)
            continue
        contracted = centroid + 0.5 * (simplex[-1] - centroid)
        contracted_value = objective(contracted)
        if contracted_value < values[-1]:
            simplex[-1], values[-1] = contracted, contracted_value
            continue
        best = simplex[0]
        simplex = [best + 0.5 * (point - best) for point in simplex]
        values = [objective(point) for point in simplex]
    order = np.argsort(values)
    return simplex[int(order[0])], float(values[int(order[0])])


def initial_log_params(rows: list[dict[str, object]], lgs: tuple[str, ...], model: str) -> np.ndarray:
    positive = np.asarray([float(r["T"]) for r in rows if float(r["T"]) > 0], dtype=float)
    scale = 1.0 / max(float(np.median(positive)) if positive.size else 1.0, 1e-4)
    if model == "shared_lambda_lg_m":
        return np.log(np.repeat(scale, 1 + len(lgs)))
    if model == "lg_lambda_lg_m":
        return np.log(np.repeat(scale, 2 * len(lgs)))
    return np.log([scale, scale])


def fit_model(rows: list[dict[str, object]], lgs: tuple[str, ...], model: str, start: np.ndarray | None = None, max_iter: int = 700) -> dict[str, object]:
    arrays = observations_to_arrays(rows, lgs)
    if start is None:
        start = initial_log_params(rows, lgs, model)
    starts = [np.asarray(start, dtype=float)]
    for shift in (-1.0, 1.0):
        starts.append(np.asarray(start, dtype=float) + shift)
    best_x, best_obj = None, math.inf
    for candidate in starts:
        x, obj = nelder_mead(lambda p: -log_likelihood(p, arrays, lgs, model), candidate, max_iter=max_iter)
        if obj < best_obj:
            best_x, best_obj = x, obj
    assert best_x is not None
    return {
        "model": model, "log_params": best_x, "params": np.exp(np.clip(best_x, -20.0, 20.0)),
        "log_likelihood": -best_obj, "n_observations": len(rows), "arrays": arrays,
    }


def parameter_rows(fit: dict[str, object], label: str, lgs: tuple[str, ...], extra: dict[str, object] | None = None) -> list[dict[str, object]]:
    params = np.asarray(fit["params"], dtype=float)
    out = []
    model = str(fit["model"])
    if model == "shared_lambda_lg_m":
        names = ["lambda", *[f"m_{lg}" for lg in lgs]]
    elif model == "lg_lambda_lg_m":
        names = [*[f"lambda_{lg}" for lg in lgs], *[f"m_{lg}" for lg in lgs]]
    else:
        names = ["lambda", "m"]
    for name, value in zip(names, params):
        row = {"fit_label": label, "model": model, "parameter": name, "estimate": float(value), "log_likelihood": fit["log_likelihood"], "n_observations": fit["n_observations"]}
        if extra:
            row.update(extra)
        out.append(row)
    return out


def make_cv_assignments(pair_rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], dict[str, int]]:
    windows = defaultdict(list)
    for row in pair_rows:
        if row["region_class"] == "inside":
            windows[(str(row["lg"]), str(row["window_id"]))].append(row)
    assignments = {}
    out = []
    for lg in LGS:
        ordered = sorted((key for key in windows if key[0] == lg), key=lambda key: (min(r["start"] for r in windows[key]), key[1]))
        n = len(ordered)
        for index, key in enumerate(ordered):
            fold = min(4, (5 * index) // n)
            assignments[key[1]] = fold
            first = windows[key][0]
            out.append({"lg": lg, "window_id": key[1], "start": first["start"], "end": first["end"], "rank_within_lg": index, "n_inside_windows_lg": n, "fold": fold, "assignment_rule": "sorted_by_start_contiguous_equal_count"})
    return out, assignments


def observed_scores() -> dict[str, dict[str, object]]:
    return {row["window_id"]: {**row, "D": float(row["D_arrangement_minus_baseline"])} for row in read_tsv(SCORE_PATH)}


def split_set(split: str) -> tuple[frozenset[str], frozenset[str]]:
    left, right = split.split("|")
    return frozenset(left.split(",")), frozenset(right.split(","))


def split_equal(a: str, b: str) -> bool:
    x, y = split_set(a), split_set(b)
    return x == y or x == y[::-1]


def mrca_time(root, ages: dict[str, float], a: str, b: str) -> float:
    node = S8.mrca_node(root, a, b)
    return float(ages[node.id])


def quartet_interval(root, ages: dict[str, float], baseline_split: str) -> float:
    left, right = split_set(baseline_split)
    if len(left) != 2 or len(right) != 2:
        raise ValueError("Expected a four-taxon 2:2 split")
    left_pair = tuple(sorted(left))
    right_pair = tuple(sorted(right))
    between_a = next(iter(left))
    between_b = next(iter(right))
    root_age = mrca_time(root, ages, between_a, between_b)
    child_age = max(mrca_time(root, ages, *left_pair), mrca_time(root, ages, *right_pair))
    return max(0.0, root_age - child_age)


def A_m(t: float, lam: float, m: float) -> float:
    if t <= 0:
        return 0.0
    if m <= 1e-12:
        return 1.0 - math.exp(-2.0 * lam * t)
    rate = lam + 2.0 * m
    return lam / rate * (-math.expm1(-2.0 * rate * t))


def predict_quartets_for_window(rows: list[dict[str, str]], root, ages: dict[str, float], lam: float, m: float) -> dict[str, object]:
    counts = Counter()
    pred = []
    exclusions = Counter()
    for row in rows:
        states = [row[f"state_{i}"] for i in range(1, 5)] if "state_1" in row else []
        if not states:
            # The local-quartet table does not carry states; reconstruct the
            # 2:2 pattern from the independent arrangement split.
            arrangement_left, arrangement_right = split_set(row["arrangement_split"])
            taxa = [row[f"taxon{i}"] for i in range(1, 5)]
            states = ["0" if taxon in arrangement_left else "1" for taxon in taxa]
        if sorted(Counter(states).values()) != [2, 2]:
            exclusions["arrangement_not_2x2"] += 1
            continue
        if split_equal(row["arrangement_split"], row["baseline_split"]):
            exclusions["arrangement_equals_baseline"] += 1
            continue
        if row["informative_for_competing_predictions"] != "true":
            exclusions["not_informative_flag_other"] += 1
            continue
        try:
            t = quartet_interval(root, ages, row["baseline_split"])
        except (KeyError, ValueError):
            exclusions["baseline_interval_unavailable"] += 1
            continue
        asymmetry = A_m(t, lam, m)
        arrangement_prob = (1.0 + 2.0 * asymmetry) / 3.0
        minor_prob = (1.0 - asymmetry) / 3.0
        local_class = row["local_class"]
        if local_class not in {"arrangement", "baseline", "third"}:
            exclusions["observed_unresolved"] += 1
            continue
        counts[local_class] += 1
        pred.append({
            "asymmetry": asymmetry, "arrangement_prob": arrangement_prob,
            "baseline_prob": minor_prob, "third_prob": minor_prob,
            "local_class": local_class, "t": t,
        })
    if not pred:
        return {"n_eligible": 0, "n_excluded": sum(exclusions.values()), "exclusions": exclusions}
    observed_d = (counts["arrangement"] - counts["baseline"]) / len(pred)
    predicted_d = float(np.mean([r["asymmetry"] for r in pred]))
    return {
        "n_eligible": len(pred), "n_excluded": sum(exclusions.values()), "exclusions": exclusions,
        "n_arrangement": counts["arrangement"], "n_baseline": counts["baseline"], "n_third": counts["third"],
        "observed_D_eligible": observed_d, "predicted_D_MSRC": predicted_d,
        "predicted_q_arrangement": float(np.mean([r["arrangement_prob"] for r in pred])),
        "predicted_q_baseline": float(np.mean([r["baseline_prob"] for r in pred])),
        "predicted_q_third": float(np.mean([r["third_prob"] for r in pred])),
        "mean_t": float(np.mean([r["t"] for r in pred])),
        "predictions": pred,
    }


def score_metrics(rows: list[dict[str, object]], observed_key: str = "observed_D_eligible", predicted_key: str = "predicted_D_MSRC") -> dict[str, float]:
    rows = [r for r in rows if r.get(observed_key) is not None and math.isfinite(float(r[observed_key]))]
    if not rows:
        return {"n_windows": 0, "mae": math.nan, "rmse": math.nan, "correlation": math.nan, "mean_error": math.nan}
    obs = np.asarray([float(r[observed_key]) for r in rows])
    pred = np.asarray([float(r[predicted_key]) for r in rows])
    error = pred - obs
    corr = float(np.corrcoef(obs, pred)[0, 1]) if np.std(obs) > 0 and np.std(pred) > 0 else math.nan
    return {"n_windows": len(rows), "mae": float(np.mean(np.abs(error))), "rmse": float(np.sqrt(np.mean(error * error))), "correlation": corr, "mean_error": float(np.mean(error))}


def ordinary_msc_prediction(t: float, lam: float) -> tuple[float, float, float, float]:
    no_coal = math.exp(-lam * t)
    base = 1.0 - 2.0 * no_coal / 3.0
    minor = no_coal / 3.0
    return base, minor, minor, minor - base


def prediction_rows_for_fold(pair_rows, local_rows, assignment, fold, fit, divergence, baseline_fit):
    params = np.asarray(fit["params"], dtype=float)
    lam = float(params[0])
    m_by_lg = {lg: float(params[i + 1]) for i, lg in enumerate(LGS)}
    score = observed_scores()
    by_window = defaultdict(list)
    for row in local_rows:
        if row["window_id"] in assignment and assignment[row["window_id"]] == fold:
            by_window[row["window_id"]].append(row)
    out = []
    for window_id in sorted(by_window, key=lambda wid: (by_window[wid][0]["lg"], int(by_window[wid][0]["start"]))):
        first = by_window[window_id][0]
        result = predict_quartets_for_window(by_window[window_id], baseline_fit.rooted_root, baseline_fit.ages, lam, m_by_lg[str(first["lg"])])
        if not result.get("n_eligible"):
            continue
        eligible_pred = result["predictions"]
        msc_d = []
        for q in eligible_pred:
            msc_d.append(ordinary_msc_prediction(float(q["t"]), lam)[3])
        observed_full = score.get(window_id, {}).get("D", math.nan)
        row = {
            "fold": fold, "lg": first["lg"], "window_id": window_id, "start": first["start"], "end": first["end"],
            "lambda_fit": lam, "m_fit": m_by_lg[str(first["lg"])], "n_eligible_quartets": result["n_eligible"],
            "n_excluded_quartets": result["n_excluded"], "exclusion_counts": json.dumps(dict(result["exclusions"]), sort_keys=True),
            "observed_D_full": observed_full, "observed_D_eligible": result["observed_D_eligible"],
            "predicted_D_MSRC": result["predicted_D_MSRC"], "predicted_D_MSC": float(np.mean(msc_d)),
            "predicted_q_arrangement": result["predicted_q_arrangement"], "predicted_q_baseline": result["predicted_q_baseline"],
            "predicted_q_third": result["predicted_q_third"], "mean_baseline_interval_t": result["mean_t"],
        }
        row["MSRC_log_predictive_probability"] = float(np.mean([math.log(max(q["arrangement_prob"] if p["local_class"] == "arrangement" else q["baseline_prob"] if p["local_class"] == "baseline" else q["third_prob"], 1e-300)) for p, q in zip(eligible_pred, eligible_pred)]))
        # Reuse the same symmetric probabilities for the MSRC log score and
        # the ordinary MSC probabilities for the comparator.
        msrc_logs, msc_logs = [], []
        for q in eligible_pred:
            msrc_p = q["arrangement_prob"] if q["local_class"] == "arrangement" else q["baseline_prob"] if q["local_class"] == "baseline" else q["third_prob"]
            msc_base, msc_arr, msc_third, _ = ordinary_msc_prediction(float(q["t"]), lam)
            msc_p = msc_arr if q["local_class"] == "arrangement" else msc_base if q["local_class"] == "baseline" else msc_third
            msrc_logs.append(math.log(max(msrc_p, 1e-300)))
            msc_logs.append(math.log(max(msc_p, 1e-300)))
        row["MSRC_log_predictive_probability"] = float(np.mean(msrc_logs))
        row["MSC_log_predictive_probability"] = float(np.mean(msc_logs))
        out.append(row)
    return out


def block_chunks(rows: list[dict[str, object]], n_chunks: int = 12) -> dict[str, list[list[str]]]:
    by_lg = defaultdict(list)
    for row in rows:
        by_lg[str(row["lg"])].append((str(row["window_id"]), int(row["start"])))
    chunks = {}
    for lg in LGS:
        ordered = [wid for wid, _ in sorted(set(by_lg[lg]), key=lambda x: x[1])]
        if not ordered:
            continue
        n = min(n_chunks, len(ordered))
        chunks[lg] = [ordered[i * len(ordered) // n:(i + 1) * len(ordered) // n] for i in range(n)]
    return chunks


def bootstrap_parameters(obs_rows, full_fit, lgs):
    rng = random.Random(SEED)
    chunks = block_chunks(obs_rows)
    by_window = defaultdict(list)
    for row in obs_rows:
        by_window[(str(row["lg"]), str(row["window_id"]))].append(row)
    output = []
    for replicate in range(24):
        selected = []
        for lg in LGS:
            lg_chunks = chunks.get(lg, [])
            if not lg_chunks:
                continue
            for _ in range(len(lg_chunks)):
                selected_windows = rng.choice(lg_chunks)
                selected.extend(row for wid in selected_windows for row in by_window[(lg, wid)])
        if not selected:
            continue
        fit = fit_model(selected, lgs, "shared_lambda_lg_m", start=np.asarray(full_fit["log_params"]), max_iter=450)
        values = np.asarray(fit["params"], dtype=float)
        output.append({"replicate": replicate, "lambda": values[0], **{f"m_{lg}": values[i + 1] for i, lg in enumerate(lgs)}})
    return output


def bootstrap_prediction_metrics(predictions):
    rng = random.Random(SEED + 1)
    chunks = block_chunks(predictions, n_chunks=10)
    by_window = defaultdict(list)
    for row in predictions:
        by_window[(str(row["lg"]), str(row["window_id"]))].append(row)
    output = []
    for replicate in range(200):
        sample = []
        for lg in LGS:
            lg_chunks = chunks.get(lg, [])
            if not lg_chunks:
                continue
            for _ in range(len(lg_chunks)):
                for wid in rng.choice(lg_chunks):
                    sample.extend(by_window[(lg, wid)])
        metrics = score_metrics(sample)
        if metrics["n_windows"]:
            output.append({"replicate": replicate, **metrics})
    return output


def simulation_draw(rng: random.Random, lam: float, m: float, start_state: str) -> float:
    state = start_state
    elapsed = 0.0
    while True:
        if state == "S":
            total = lam + 2.0 * m
            elapsed += rng.expovariate(total)
            if rng.random() < lam / total:
                return elapsed
            state = "D"
        else:
            if m <= 0:
                return math.inf
            elapsed += rng.expovariate(2.0 * m)
            state = "S"


def run_simulation_validation() -> list[dict[str, object]]:
    rng = random.Random(SEED)
    scenarios = [(0.5, 0.2, "S"), (0.5, 0.2, "D"), (1.0, 0.01, "S"), (1.0, 0.01, "D"), (2.0, 5.0, "S"), (2.0, 5.0, "D")]
    rows = []
    for lam, m, state in scenarios:
        simulated = [simulation_draw(rng, lam, m, state) for _ in range(1800)]
        fit_rows = [{"T": value, "d": 0.0, "state_sign": 1.0 if state == "S" else -1.0, "lg": "SIM"} for value in simulated]
        fit = fit_model(fit_rows, ("SIM",), "one_lambda_one_m", max_iter=550)
        estimates = np.asarray(fit["params"])
        true_ll = log_likelihood(np.log([lam, m]), observations_to_arrays(fit_rows, ("SIM",)), ("SIM",), "one_lambda_one_m")
        rows.append({
            "validation": "synthetic_ctmc_recovery", "lambda_true": lam, "m_true": m, "start_state": state,
            "lambda_estimate": estimates[0], "m_estimate": estimates[1], "lambda_relative_error": estimates[0] / lam - 1.0,
            "m_relative_error": estimates[1] / m - 1.0, "loglik_at_estimate": fit["log_likelihood"], "loglik_at_truth": true_ll,
            "likelihood_improves_near_truth": fit["log_likelihood"] >= true_ll - 1e-5, "positive_estimates": bool(np.all(estimates > 0)),
            "n": len(simulated), "seed": SEED,
        })
    # Both starting states are needed to separate lambda from m in the
    # small-m regime. This is an explicit identifiability check.
    for lam, m in ((0.5, 0.2), (1.0, 0.01), (2.0, 5.0)):
        mixed = []
        for state, sign in (("S", 1.0), ("D", -1.0)):
            mixed.extend({"T": simulation_draw(rng, lam, m, state), "d": 0.0, "state_sign": sign, "lg": "SIM"} for _ in range(1800))
        fit = fit_model(mixed, ("SIM",), "one_lambda_one_m", start=np.log([lam, m]), max_iter=750)
        estimates = np.asarray(fit["params"])
        true_ll = log_likelihood(np.log([lam, m]), observations_to_arrays(mixed, ("SIM",)), ("SIM",), "one_lambda_one_m")
        rows.append({
            "validation": "synthetic_ctmc_mixed_start_recovery", "lambda_true": lam, "m_true": m, "start_state": "S+D",
            "lambda_estimate": estimates[0], "m_estimate": estimates[1], "lambda_relative_error": estimates[0] / lam - 1.0,
            "m_relative_error": estimates[1] / m - 1.0, "loglik_at_estimate": fit["log_likelihood"], "loglik_at_truth": true_ll,
            "likelihood_improves_near_truth": fit["log_likelihood"] >= true_ll - 1e-5, "positive_estimates": bool(np.all(estimates > 0)),
            "n": len(mixed), "seed": SEED,
        })
    for m in (1e-8, 1e6):
        rows.append({"validation": "quartet_limit", "lambda_true": 1.0, "m_true": m, "start_state": "1010_special_case", "lambda_estimate": "", "m_estimate": "", "lambda_relative_error": "", "m_relative_error": "", "loglik_at_estimate": "", "loglik_at_truth": "", "likelihood_improves_near_truth": A_m(0.4, 1.0, m), "positive_estimates": "", "n": "", "seed": SEED})
    return rows


def make_summary(predictions, pred_bootstrap, parameter_bootstrap, full_fit, sensitivity_fit, audit, cv_assignments):
    summary = []
    all_metrics = score_metrics(predictions)
    msc_metrics = score_metrics(predictions, predicted_key="predicted_D_MSC")
    summary.append({"summary_type": "all_LGs", "lg": "ALL", "n_windows": all_metrics["n_windows"], "mae_MSRC": all_metrics["mae"], "rmse_MSRC": all_metrics["rmse"], "correlation_MSRC": all_metrics["correlation"], "mean_error_MSRC": all_metrics["mean_error"], "mae_MSC": msc_metrics["mae"], "rmse_MSC": msc_metrics["rmse"], "correlation_MSC": msc_metrics["correlation"], "mean_error_MSC": msc_metrics["mean_error"], "MSRC_mean_log_predictive_probability": float(np.mean([r["MSRC_log_predictive_probability"] for r in predictions])), "MSC_mean_log_predictive_probability": float(np.mean([r["MSC_log_predictive_probability"] for r in predictions])), "bootstrap_mae_low": np.quantile([r["mae"] for r in pred_bootstrap], .025) if pred_bootstrap else math.nan, "bootstrap_mae_high": np.quantile([r["mae"] for r in pred_bootstrap], .975) if pred_bootstrap else math.nan, "bootstrap_rmse_low": np.quantile([r["rmse"] for r in pred_bootstrap], .025) if pred_bootstrap else math.nan, "bootstrap_rmse_high": np.quantile([r["rmse"] for r in pred_bootstrap], .975) if pred_bootstrap else math.nan})
    for lg in LGS:
        p = [r for r in predictions if r["lg"] == lg]
        a = score_metrics(p)
        m = score_metrics(p, predicted_key="predicted_D_MSC")
        summary.append({"summary_type": "linkage_group", "lg": lg, "n_windows": a["n_windows"], "mae_MSRC": a["mae"], "rmse_MSRC": a["rmse"], "correlation_MSRC": a["correlation"], "mean_error_MSRC": a["mean_error"], "mae_MSC": m["mae"], "rmse_MSC": m["rmse"], "correlation_MSC": m["correlation"], "mean_error_MSC": m["mean_error"], "MSRC_mean_log_predictive_probability": float(np.mean([r["MSRC_log_predictive_probability"] for r in p])), "MSC_mean_log_predictive_probability": float(np.mean([r["MSC_log_predictive_probability"] for r in p])), "bootstrap_mae_low": math.nan, "bootstrap_mae_high": math.nan, "bootstrap_rmse_low": math.nan, "bootstrap_rmse_high": math.nan})
    summary.append({"summary_type": "audit", "lg": "ALL", "n_windows": len(cv_assignments), "mae_MSRC": "", "rmse_MSRC": "", "correlation_MSRC": "", "mean_error_MSRC": "", "mae_MSC": "", "rmse_MSC": "", "correlation_MSC": "", "mean_error_MSC": "", "MSRC_mean_log_predictive_probability": "", "MSC_mean_log_predictive_probability": "", "bootstrap_mae_low": "", "bootstrap_mae_high": "", "bootstrap_rmse_low": "", "bootstrap_rmse_high": "", "n_incompatible": sum(1 for r in audit if r["incompatible"]), "n_audit": len(audit), "n_parameter_bootstrap": len(parameter_bootstrap), "full_primary_log_likelihood": full_fit["log_likelihood"], "sensitivity_log_likelihood": sensitivity_fit["log_likelihood"]})
    return summary


def plot_figure(predictions, parameter_bootstrap, summary):
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.8), constrained_layout=True)
    ax = axes[0, 0]
    param_values = {"lambda": [float(r["lambda"]) for r in parameter_bootstrap]}
    for lg in LGS:
        param_values[f"m_{lg}"] = [float(r[f"m_{lg}"]) for r in parameter_bootstrap]
    estimates = {"lambda": float(summary[0].get("lambda_estimate", math.nan)) if "lambda_estimate" in summary[0] else math.nan}
    # The parameter table is easier to consume than this internal dictionary;
    # the plotted point estimates are supplied by the caller through bootstrap
    # medians when the table is not passed here.
    labels = ["lambda", *[f"m_{lg}" for lg in LGS]]
    medians = [np.median(param_values[k]) if param_values[k] else math.nan for k in labels]
    lows = [np.quantile(param_values[k], .025) if param_values[k] else math.nan for k in labels]
    highs = [np.quantile(param_values[k], .975) if param_values[k] else math.nan for k in labels]
    ax.errorbar(range(len(labels)), medians, yerr=[np.asarray(medians) - lows, np.asarray(highs) - medians], fmt="o", color="#1f4e79", capsize=3)
    ax.set_xticks(range(len(labels)), ["λ", "m01", "m02", "m07", "m12"])
    ax.set_yscale("log"); ax.set_ylabel("effective rate"); ax.set_title("A  block-bootstrap parameter estimates")
    ax.grid(axis="y", alpha=.25)

    ax = axes[0, 1]
    for lg in LGS:
        p = [r for r in predictions if r["lg"] == lg]
        ax.scatter([r["observed_D_eligible"] for r in p], [r["predicted_D_MSRC"] for r in p], s=12, alpha=.55, label=lg)
    ax.axline((0, 0), slope=1, color="black", lw=.8, ls="--")
    ax.set_xlabel("observed eligible quartet contrast D"); ax.set_ylabel("MSRC-predicted D"); ax.set_title("B  held-out quartet contrast")
    ax.legend(frameon=False, fontsize=8); ax.grid(alpha=.2)

    ax = axes[1, 0]
    for lg in LGS:
        p = sorted([r for r in predictions if r["lg"] == lg], key=lambda r: int(r["start"]))
        if not p: continue
        x = np.asarray([int(r["start"]) for r in p])
        ax.plot(x, [r["observed_D_eligible"] for r in p], lw=.8, alpha=.75, label=f"{lg} observed")
        ax.plot(x, [r["predicted_D_MSRC"] for r in p], lw=1.2, label=f"{lg} MSRC")
    ax.set_xlabel("window start (bp)"); ax.set_ylabel("quartet contrast D"); ax.set_title("C  genomic held-out tracks")
    ax.legend(frameon=False, ncol=2, fontsize=7); ax.grid(alpha=.2)

    ax = axes[1, 1]
    all_row = next(r for r in summary if r["summary_type"] == "all_LGs")
    metric_names = ["MAE", "RMSE"]
    msrc = [float(all_row["mae_MSRC"]), float(all_row["rmse_MSRC"])]
    msc = [float(all_row["mae_MSC"]), float(all_row["rmse_MSC"])]
    x = np.arange(2); width = .34
    ax.bar(x - width / 2, msrc, width, label="MSRC", color="#1f4e79")
    ax.bar(x + width / 2, msc, width, label="MSC", color="#a65f2b")
    ax.set_xticks(x, metric_names); ax.set_ylabel("held-out prediction error"); ax.set_title("D  MSRC versus ordinary MSC")
    ax.legend(frameon=False); ax.grid(axis="y", alpha=.2)
    fig.suptitle("Atlantic cod Stage 10: coalescence-time fit to held-out quartet prediction", fontsize=13)
    fig.savefig(FIG_PDF); fig.savefig(FIG_PNG, dpi=300); plt.close(fig)


def write_text_outputs(summary, audit, parameter_rows_all, sim_rows, fit, sensitivity_fit, exclusion_counts):
    incompatible = sum(1 for r in audit if r["incompatible"])
    total = len(audit)
    all_row = next(r for r in summary if r["summary_type"] == "all_LGs")
    METHODS_OUT.write_text(textwrap.dedent(f"""
    # Stage 10 methods

    Stage 10 fits a simplified symmetric two-background MSRC to published SNAPP pairwise MRCA times. The topology is the independently frozen collinear population tree in `data/atlantic_cod/processed/baseline_population_tree_12taxa.nwk`; its original branch lengths are ignored. Outside-supergene windows (`region_class=outside`) provide pairwise median MRCA targets for a constrained monotone node-age projection. The root edge is selected by minimum outside-only projection loss, and the resulting time-calibrated tree is frozen in `stage10_baseline_time_tree.nwk`.

    For a pair entering the shared population at divergence time d, present-day states are propagated backward for d using the symmetric two-state switching process. If the extant states agree, the probability of entering in S is `(1 + exp(-4 m d))/2`; if they differ it is `(1 - exp(-4 m d))/2`. After entry, the transient generator is `Q=[[-(lambda+2m),2m],[2m,-2m]]` and c=`[lambda,0]`. The density is `alpha exp(Q (T-d)) c`, evaluated by the exact analytic 2-by-2 matrix exponential. Rates are optimized on the log scale.

    The fit is a pairwise composite likelihood because the same local genealogy contributes to many population pairs. The primary model has one shared lambda and linkage-group-specific m values for LG01, LG02, LG07, and LG12. Boundary windows are excluded. The {incompatible}/{total} ({incompatible / total:.3%}) observations with T < d are explicitly audited and excluded from the primary fit; a separately labelled clipped sensitivity fit replaces them by d + 1e-9 only to show the effect of this incompatibility rule.

    Cross-validation divides inside-supergene windows into five deterministic contiguous blocks within each linkage group. Each fold fits coalescence times from all outside windows plus four inside blocks. Quartet counts are read only after the fold fit. The primary quartet prediction uses the exact repository-supported symmetric 2:2 special case `A_m(t) = lambda/(lambda+2m) * [1-exp(-2(lambda+2m)t)]`, where t is the fixed baseline quartet internal interval. It predicts `(q_arrangement,q_baseline,q_third)=((1+2A_m)/3,(1-A_m)/3,(1-A_m)/3)`. Other background configurations are excluded and reported by reason.

    Uncertainty uses a deterministic contiguous-window bootstrap with seed {SEED}. The effective switching parameter m is reported in inverse SNAPP time units and is not interpreted as a meiotic recombination rate. Stage 7B remains an interpretive decomposition of observed coalescence-time contrast; Stage 10 estimates parameters from coalescence times and predicts independent held-out quartet behavior.
    """).strip() + "\n")
    RESULTS_TEXT_OUT.write_text(textwrap.dedent(f"""
    # Stage 10 results

    The primary full-data composite fit used {fit['n_observations']} compatible non-boundary pair-window observations. The held-out MSRC prediction used {all_row['n_windows']} windows and achieved MAE {float(all_row['mae_MSRC']):.4g}, RMSE {float(all_row['rmse_MSRC']):.4g}, and correlation {float(all_row['correlation_MSRC']):.4g} between predicted and observed eligible-quartet contrasts. The ordinary MSC comparator achieved MAE {float(all_row['mae_MSC']):.4g}, RMSE {float(all_row['rmse_MSC']):.4g}, and correlation {float(all_row['correlation_MSC']):.4g}; mean window-level log predictive probabilities were {float(all_row['MSRC_mean_log_predictive_probability']):.4g} for MSRC and {float(all_row['MSC_mean_log_predictive_probability']):.4g} for MSC.

    The fit has {len(parameter_rows_all)} reported parameter rows, including the primary and linkage-group-specific sensitivity models. Profile values and block-bootstrap intervals are in the TSV outputs. Simulation validation covered six known CTMC parameter/state settings; the detailed rows are in `stage10_simulation_validation.tsv`.

    Quartet eligibility was restricted to symmetric 2:2 background patterns supported by the exact available formula. Exclusion counts are {json.dumps(dict(exclusion_counts), sort_keys=True)}. The result is therefore a parameter-to-prediction test for that mathematically supported quartet subset, rather than a claim of a full 12-population MSRC quartet distribution.
    """).strip() + "\n")
    REPORT_OUT.write_text(textwrap.dedent(f"""
    # Atlantic cod Stage 10 report

    ## Implemented analysis

    Stage 10 freezes a time-calibrated version of the independent collinear population-tree topology using outside-supergene SNAPP MRCA times, fits a two-background MSRC pairwise composite likelihood, and uses the fitted coalescence-time parameters to predict held-out quartet contrast in contiguous genomic blocks.

    ## Estimates and uncertainty

    The primary estimates and sensitivity estimates are in `stage10_msrc_parameter_estimates.tsv`. The likelihood profiles are in `stage10_msrc_parameter_profiles.tsv`; bootstrap parameter replicates use seed {SEED}. Interpret m as an effective ancestry-background switching rate in inverse SNAPP time units.

    ## Identifiability and prediction

    The reported profile shapes and parameter bootstrap intervals should be used to judge whether lambda and each m are separately estimable. Held-out prediction summaries are in `stage10_prediction_summary.tsv`, with window-level predictions in `stage10_cv_predictions.tsv`. The primary result is the predictive comparison, not a p-value.

    ## Limitations

    {incompatible}/{total} pair-window observations satisfy T < d and were excluded from the primary fit. The clipped fit is an explicit sensitivity analysis. Quartet prediction is limited to the exact symmetric 2:2 `A_m(t)` case available in the repository theory. No unsupported full quartet distribution was imputed. A good fit would support the statement that coalescence-time-estimated MSRC parameters predict held-out arrangement-associated quartet bias; it would not establish MSRC uniquely.
    """).strip() + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true", help="run internal deterministic checks and exit")
    args = parser.parse_args()
    if args.self_test:
        assert math.isclose(A_m(0.4, 1.0, 0.0), 1.0 - math.exp(-0.8))
        assert A_m(0.4, 1.0, 1e6) < 1e-5
        assert matrix_exponential_density(np.array([1.0]), np.array([0.0]), np.array([1.0]), 1.0, 0.0)[0] > 0
        print("stage10 self-tests passed")
        return

    pair_rows = parse_pair_rows()
    baseline_fit, divergence, outside_targets, root_audit = fit_baseline_time_tree(pair_rows)
    compatible_rows, audit_rows = make_observation_rows(pair_rows, divergence, clipped=False)
    clipped_rows, _ = make_observation_rows(pair_rows, divergence, clipped=True)
    write_tsv(INCOMPAT_OUT, audit_rows, ["lg", "window_id", "start", "end", "population1", "population2", "observed_mrca_time", "baseline_divergence_time", "incompatibility_amount", "incompatible", "primary_status"])

    primary_fit = fit_model(compatible_rows, LGS, "shared_lambda_lg_m")
    sensitivity_fit = fit_model(compatible_rows, LGS, "lg_lambda_lg_m")
    clipped_fit = fit_model(clipped_rows, LGS, "shared_lambda_lg_m")
    param_rows = parameter_rows(primary_fit, "full_primary_incompatible_excluded", LGS)
    param_rows += parameter_rows(sensitivity_fit, "sensitivity_lg_specific_incompatible_excluded", LGS)
    param_rows += parameter_rows(clipped_fit, "sensitivity_shared_lambda_clipped_incompatible", LGS)
    profile_rows = []
    mle = np.asarray(primary_fit["log_params"], dtype=float)
    names = ["lambda", *[f"m_{lg}" for lg in LGS]]
    arrays = primary_fit["arrays"]
    for index, name in enumerate(names):
        for delta in np.linspace(-3.0, 3.0, 31):
            p = mle.copy(); p[index] += delta
            profile_rows.append({"model": "full_primary", "parameter": name, "value": math.exp(float(p[index])), "log_likelihood": log_likelihood(p, arrays, LGS, "shared_lambda_lg_m"), "delta_log_likelihood": log_likelihood(p, arrays, LGS, "shared_lambda_lg_m") - primary_fit["log_likelihood"], "other_parameters_fixed_at_mle": True})
    write_tsv(PARAM_OUT, param_rows, ["fit_label", "model", "parameter", "estimate", "log_likelihood", "n_observations"])
    write_tsv(PROFILE_OUT, profile_rows, ["model", "parameter", "value", "log_likelihood", "delta_log_likelihood", "other_parameters_fixed_at_mle"])

    assignments, assignment_map = make_cv_assignments(pair_rows)
    write_tsv(CV_ASSIGN_OUT, assignments, ["lg", "window_id", "start", "end", "rank_within_lg", "n_inside_windows_lg", "fold", "assignment_rule"])
    local_rows = read_tsv(QUARTET_PATH)
    predictions = []
    cv_param_rows = []
    outside_all = [r for r in compatible_rows if r["region_class"] == "outside"]
    for fold in range(5):
        train = outside_all + [r for r in compatible_rows if r["region_class"] == "inside" and assignment_map.get(str(r["window_id"])) != fold]
        fold_fit = fit_model(train, LGS, "shared_lambda_lg_m", start=np.asarray(primary_fit["log_params"]), max_iter=650)
        cv_param_rows += parameter_rows(fold_fit, f"cv_fold_{fold}", LGS, {"fold": fold})
        predictions.extend(prediction_rows_for_fold(pair_rows, local_rows, assignment_map, fold, fold_fit, divergence, baseline_fit))
    param_rows += cv_param_rows
    write_tsv(PARAM_OUT, param_rows, ["fit_label", "model", "parameter", "estimate", "log_likelihood", "n_observations", "fold"])
    write_tsv(CV_PRED_OUT, predictions, ["fold", "lg", "window_id", "start", "end", "lambda_fit", "m_fit", "n_eligible_quartets", "n_excluded_quartets", "exclusion_counts", "observed_D_full", "observed_D_eligible", "predicted_D_MSRC", "predicted_D_MSC", "predicted_q_arrangement", "predicted_q_baseline", "predicted_q_third", "mean_baseline_interval_t", "MSRC_log_predictive_probability", "MSC_log_predictive_probability"])

    parameter_bootstrap = bootstrap_parameters(compatible_rows, primary_fit, LGS)
    prediction_bootstrap = bootstrap_prediction_metrics(predictions)
    summary = make_summary(predictions, prediction_bootstrap, parameter_bootstrap, primary_fit, sensitivity_fit, audit_rows, assignments)
    # Add full-fit estimates and bootstrap intervals to the parameter table.
    boot_by = {"lambda": [float(r["lambda"]) for r in parameter_bootstrap], **{f"m_{lg}": [float(r[f"m_{lg}"]) for r in parameter_bootstrap] for lg in LGS}}
    for row in param_rows:
        values = boot_by.get(str(row["parameter"]), [])
        if values:
            row["bootstrap_low"] = float(np.quantile(values, .025)); row["bootstrap_high"] = float(np.quantile(values, .975))
        else:
            row["bootstrap_low"] = math.nan; row["bootstrap_high"] = math.nan
    write_tsv(PARAM_OUT, param_rows, ["fit_label", "model", "parameter", "estimate", "log_likelihood", "n_observations", "fold", "bootstrap_low", "bootstrap_high"])
    write_tsv(PRED_SUMMARY_OUT, summary, ["summary_type", "lg", "n_windows", "mae_MSRC", "rmse_MSRC", "correlation_MSRC", "mean_error_MSRC", "mae_MSC", "rmse_MSC", "correlation_MSC", "mean_error_MSC", "MSRC_mean_log_predictive_probability", "MSC_mean_log_predictive_probability", "bootstrap_mae_low", "bootstrap_mae_high", "bootstrap_rmse_low", "bootstrap_rmse_high", "n_incompatible", "n_audit", "n_parameter_bootstrap", "full_primary_log_likelihood", "sensitivity_log_likelihood"])
    sim_rows = run_simulation_validation()
    write_tsv(SIM_OUT, sim_rows, ["validation", "lambda_true", "m_true", "start_state", "lambda_estimate", "m_estimate", "lambda_relative_error", "m_relative_error", "loglik_at_estimate", "loglik_at_truth", "likelihood_improves_near_truth", "positive_estimates", "n", "seed"])
    exclusion_counts = Counter()
    for row in predictions:
        exclusion_counts.update(json.loads(row["exclusion_counts"]))
    write_text_outputs(summary, audit_rows, param_rows, sim_rows, primary_fit, sensitivity_fit, exclusion_counts)
    plot_figure(predictions, parameter_bootstrap, summary)

    outputs = [BASELINE_TREE_OUT, PROJECTION_OUT, PARAM_OUT, PROFILE_OUT, CV_ASSIGN_OUT, CV_PRED_OUT, PRED_SUMMARY_OUT, INCOMPAT_OUT, SIM_OUT, METHODS_OUT, RESULTS_TEXT_OUT, REPORT_OUT, FIG_PDF, FIG_PNG]
    manifest = {
        "stage": 10, "analysis": "Atlantic cod simplified two-background MSRC parameter-to-prediction test", "seed": SEED,
        "inputs": {str(path.relative_to(ROOT)): sha256(path) for path in [BASELINE_TOPOLOGY, MRCA_PATH, QUARTET_PATH, SCORE_PATH, ARRANGEMENT_PATH, REGION_PATH]},
        "outputs": {str(path.relative_to(ROOT)): sha256(path) for path in outputs},
        "n_pair_rows_raw": len(pair_rows), "n_pair_rows_primary": len(compatible_rows), "n_incompatible": sum(1 for r in audit_rows if r["incompatible"]),
        "n_cv_predictions": len(predictions), "n_cv_assignment_windows": len(assignments), "quartet_prediction_theory": "symmetric_2x2_A_m_special_case_only",
    }
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"primary_log_likelihood": primary_fit["log_likelihood"], "n_primary": len(compatible_rows), "n_incompatible": sum(1 for r in audit_rows if r["incompatible"]), "n_predictions": len(predictions), "outputs": len(outputs)}, indent=2))


if __name__ == "__main__":
    main()
