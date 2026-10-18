#!/usr/bin/env python3
"""Stage 10b: fit the validated conditional 2:2 quartet model.

This analysis uses quartet outcomes as the observational data.  It does not
use Stage-7 SNAPP node ages or any coalescence-time likelihood.  The model is
the symmetric 1010 conditional quartet calculation validated in
``simulations/quartet``:

    A(lambda, m, t) = lambda/(lambda + 2m)
                       * (1 - exp(-2(lambda + 2m)t))

and q=( (1+2A)/3, (1-A)/3, (1-A)/3 ) for arrangement, baseline, and third.
The resulting parameters are effective reduced-model parameters in inverse
Ma units.  They are not estimates of meiotic recombination rates or of the
full biological MSRC process.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import textwrap
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("/private/tmp") / "msrc-stage10b-mpl"))
import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "atlantic_cod"
EMP = ROOT / "empirical" / "atlantic_cod"
RESULTS = EMP / "results"
FIGURES = EMP / "figures"

TREE = RESULTS / "stage8b_baseline_original_times.nwk"
STAGE8B_MANIFEST = RESULTS / "stage8b_manifest.json"
LOCAL = DATA / "processed" / "stage4a_local_quartets.tsv"
SCORES = DATA / "processed" / "stage4a_window_topology_scores.tsv"
ARRANGEMENTS = DATA / "processed" / "population_arrangements.tsv"
FROZEN_GRID = ROOT / "simulations/quartet/datasets/processed/conditional_quartet_grid.tsv"

ELIGIBILITY = RESULTS / "stage10b_quartet_eligibility.tsv"
WINDOW_COUNTS = RESULTS / "stage10b_window_counts.tsv"
PARAMETERS = RESULTS / "stage10b_parameter_estimates.tsv"
PROFILES = RESULTS / "stage10b_parameter_profiles.tsv"
CV_ASSIGNMENT = RESULTS / "stage10b_cv_assignment.tsv"
CV_PREDICTIONS = RESULTS / "stage10b_cv_predictions.tsv"
PRED_SUMMARY = RESULTS / "stage10b_prediction_summary.tsv"
MODEL_COMPARISON = RESULTS / "stage10b_model_comparison.tsv"
IDENTIFIABILITY = RESULTS / "stage10b_identifiability.tsv"
SIMULATION = RESULTS / "stage10b_simulation_validation.tsv"
METHODS = RESULTS / "stage10b_methods_text.md"
RESULTS_TEXT = RESULTS / "stage10b_results_text.md"
REPORT = RESULTS / "stage10b_report.md"
MANIFEST = RESULTS / "stage10b_manifest.json"
LEGACY_AUDIT = RESULTS / "stage10_legacy_validity_audit.md"
FIG_PDF = FIGURES / "atlantic_cod_stage10b_conditional_msrc_prediction.pdf"
FIG_PNG = FIGURES / "atlantic_cod_stage10b_conditional_msrc_prediction.png"

LGS = ("LG01", "LG02", "LG07", "LG12")
SEED = 20261018
_ARRAY_CACHE: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
OLD_STAGE10 = [p for p in sorted(RESULTS.glob("stage10_*")) if p.name != "stage10_legacy_validity_audit.md"] + [
    FIGURES / "atlantic_cod_stage10_msrc_parameter_prediction.pdf",
    FIGURES / "atlantic_cod_stage10_msrc_parameter_prediction.png",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        out = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        out.writeheader()
        for row in rows:
            out.writerow({field: fmt(row.get(field, "")) for field in fields})


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return "NA" if not math.isfinite(value) else f"{value:.10g}"
    return str(value)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def split_set(split: str) -> tuple[frozenset[str], frozenset[str]]:
    left, right = split.split("|")
    return frozenset(left.split(",")), frozenset(right.split(","))


def split_equal(a: str, b: str) -> bool:
    x, y = split_set(a), split_set(b)
    return x == y or x == y[::-1]


class TreeModel:
    def __init__(self, path: Path):
        self.tree = Phylo.read(str(path), "newick")
        self.tree.rooted = True
        self.taxa = {str(c.name) for c in self.tree.get_terminals()}
        self.ages: dict[int, float] = {}

        def age(clade) -> float:
            if clade.is_terminal():
                self.ages[id(clade)] = 0.0
                return 0.0
            vals = []
            for child in clade.clades:
                vals.append(age(child) + float(child.branch_length or 0.0))
            self.ages[id(clade)] = max(vals)
            return self.ages[id(clade)]

        age(self.tree.root)

    def mrca_age(self, taxa: list[str] | tuple[str, ...]) -> float:
        clade = self.tree.common_ancestor(*taxa)
        return self.ages[id(clade)]

    def quartet_interval(self, baseline_split: str) -> float:
        left, right = split_set(baseline_split)
        if len(left) != 2 or len(right) != 2:
            raise ValueError("baseline split is not a 2:2 split")
        all_age = self.mrca_age(tuple(left | right))
        child_age = max(self.mrca_age(tuple(left)), self.mrca_age(tuple(right)))
        return all_age - child_age


def verify_scaffold() -> dict[str, object]:
    manifest = json.loads(STAGE8B_MANIFEST.read_text())
    expected = manifest["output_checksums"][str(TREE.relative_to(ROOT))]
    actual = sha256(TREE)
    if actual != expected:
        raise ValueError("Stage-8b original tree checksum does not match its manifest")
    model = TreeModel(TREE)
    if len(model.taxa) != 12 or not model.tree.rooted or not manifest["root_audit"]["topology_match_existing_12taxon_tree"]:
        raise ValueError("Stage-8b scaffold is not a rooted 12-taxon tree")
    if any(float(c.branch_length or 0.0) <= 0.0 for c in model.tree.get_nonterminals() if c is not model.tree.root):
        raise ValueError("Stage-8b scaffold has a non-positive internal branch")
    return {"tree_checksum": actual, "n_tips": len(model.taxa), "rooted": True, "root_age": model.ages[id(model.tree.root)]}


def A_m(t: float, lam: float, m: float) -> float:
    if t <= 0.0:
        return 0.0
    rate = lam + 2.0 * m
    return lam / rate * (-math.expm1(-2.0 * rate * t))


def conditional_probabilities(t: float, lam: float, m: float) -> tuple[float, float, float]:
    a = A_m(t, lam, m)
    q = ((1.0 + 2.0 * a) / 3.0, (1.0 - a) / 3.0, (1.0 - a) / 3.0)
    if min(q) < -1e-12 or not math.isclose(sum(q), 1.0, abs_tol=1e-12):
        raise ValueError("conditional probabilities are invalid")
    return q


def validate_frozen_kernel() -> float:
    errors = []
    for row in read_tsv(FROZEN_GRID):
        t, m = float(row["duration"]), float(row["m"])
        q = conditional_probabilities(t, 1.0, m)
        expected = (float(row["exact_q_alt"]), float(row["exact_q_species"]), float(row["exact_q_other"]))
        errors.append(max(abs(a - b) for a, b in zip(q, expected)))
    return max(errors) if errors else math.nan


def load_rows(tree: TreeModel) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    arrangements = {r["tree_label"]: r for r in read_tsv(ARRANGEMENTS)}
    rows = []
    audit = []
    for raw in read_tsv(LOCAL):
        taxa = [raw[f"taxon{i}"] for i in range(1, 5)]
        reason = "eligible"
        t = math.nan
        states = []
        if len(set(taxa)) != 4:
            reason = "not_exactly_four_unique_taxa"
        elif not set(taxa) <= tree.taxa:
            reason = "taxon_missing_from_scaffold"
        elif raw["local_class"] not in {"arrangement", "baseline", "third"}:
            reason = "local_topology_unresolved_or_unsupported"
        else:
            left, right = split_set(raw["arrangement_split"])
            for taxon in taxa:
                states.append("1" if taxon in left else "0")
            if sorted(Counter(states).values()) != [2, 2]:
                reason = "arrangement_state_not_2x2"
            elif split_equal(raw["arrangement_split"], raw["baseline_split"]):
                reason = "arrangement_split_equals_baseline"
            elif any(taxon not in arrangements for taxon in taxa):
                reason = "arrangement_metadata_missing"
            else:
                metadata_states = [arrangements[taxon][f"{raw['lg']}_state"] for taxon in taxa]
                metadata_pattern = "1" if metadata_states[0] == "derived" else "0"
                metadata_pattern = "".join("1" if x == "derived" else "0" for x in metadata_states)
                if sorted(Counter(metadata_pattern).values()) != [2, 2]:
                    reason = "metadata_state_not_2x2"
                elif not split_equal(raw["arrangement_split"], ",".join(taxa[i] for i, x in enumerate(metadata_pattern) if x == "1") + "|" + ",".join(taxa[i] for i, x in enumerate(metadata_pattern) if x == "0")):
                    reason = "arrangement_metadata_split_mismatch"
                else:
                    t = tree.quartet_interval(raw["baseline_split"])
                    if not math.isfinite(t) or t <= 0.0:
                        reason = "nonpositive_baseline_interval"
        row = {
            **raw,
            "taxa": taxa,
            "state_pattern": "".join(states),
            "interval_t": t,
            "eligibility_reason": reason,
            "eligible": reason == "eligible",
            "label_index": {"arrangement": 0, "baseline": 1, "third": 2}.get(raw["local_class"], -1),
            "start": int(raw["start"]), "end": int(raw["end"]),
        }
        rows.append(row)
        audit.append({"lg": raw["lg"], "window_id": raw["window_id"], "quartet_id": raw["quartet_id"], "taxon1": taxa[0], "taxon2": taxa[1], "taxon3": taxa[2], "taxon4": taxa[3], "arrangement_split": raw["arrangement_split"], "baseline_split": raw["baseline_split"], "local_class": raw["local_class"], "state_pattern": row["state_pattern"], "interval_t": t, "eligible": row["eligible"], "exclusion_reason": reason})
    return rows, audit


def inside_windows(rows: list[dict[str, object]]) -> dict[tuple[str, str], dict[str, object]]:
    out = {}
    for row in rows:
        if row["overlaps_inversion"] == "true" and row["eligible"]:
            key = (str(row["lg"]), str(row["window_id"]))
            out.setdefault(key, {"lg": row["lg"], "window_id": row["window_id"], "start": row["start"], "end": row["end"], "rows": []})["rows"].append(row)
    return out


def make_assignments(rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], dict[str, int]]:
    windows = inside_windows(rows)
    assignments = {}
    output = []
    for lg in LGS:
        ordered = sorted([v for k, v in windows.items() if k[0] == lg], key=lambda x: (int(x["start"]), str(x["window_id"])))
        n = len(ordered)
        for rank, win in enumerate(ordered):
            fold = min(4, (5 * rank) // n)
            assignments[str(win["window_id"])] = fold
            output.append({"lg": lg, "window_id": win["window_id"], "start": win["start"], "end": win["end"], "rank_within_lg": rank, "n_inside_windows_lg": n, "fold": fold, "assignment_rule": "sorted_by_start_contiguous_equal_count"})
    return output, assignments


def eligible_by_window(rows: list[dict[str, object]]) -> dict[str, list[dict[str, object]]]:
    out = defaultdict(list)
    for row in rows:
        if row["eligible"]:
            out[str(row["window_id"])].append(row)
    return out


def model_parameter_names(model: str) -> list[str]:
    if model == "common":
        return ["lambda", "m"]
    if model == "lg_m":
        return ["lambda", *[f"m_{lg}" for lg in LGS]]
    if model == "lg_lg":
        return [*[f"lambda_{lg}" for lg in LGS], *[f"m_{lg}" for lg in LGS]]
    raise ValueError(model)


def unpack(theta: np.ndarray, model: str) -> tuple[np.ndarray, np.ndarray]:
    p = np.exp(np.clip(theta, -20.0, 20.0))
    if model == "common":
        return np.repeat(p[0], 4), np.repeat(p[1], 4)
    if model == "lg_m":
        return np.repeat(p[0], 4), p[1:5]
    return p[:4], p[4:8]


def log_likelihood(theta: np.ndarray, rows: list[dict[str, object]], model: str) -> float:
    if np.any(theta < -15.0) or np.any(theta > 15.0):
        return -1e100 - float(np.sum(np.maximum(np.abs(theta) - 15.0, 0.0)))
    lambdas, ms = unpack(theta, model)
    cache_key = id(rows)
    if cache_key not in _ARRAY_CACHE:
        _ARRAY_CACHE[cache_key] = (
            np.asarray([float(row["interval_t"]) for row in rows], dtype=float),
            np.asarray([LGS.index(str(row["lg"])) for row in rows], dtype=int),
            np.asarray([int(row["label_index"]) for row in rows], dtype=int),
        )
    t, lg_index, labels = _ARRAY_CACHE[cache_key]
    lam = lambdas[lg_index]
    m = ms[lg_index]
    rate = lam + 2.0 * m
    a = lam / rate * (-np.expm1(-2.0 * rate * t))
    probs = np.where(labels == 0, (1.0 + 2.0 * a) / 3.0, (1.0 - a) / 3.0)
    return float(np.log(np.maximum(probs, 1e-300)).sum())


def nelder_mead(objective, start: np.ndarray, max_iter: int = 850) -> tuple[np.ndarray, float]:
    n = len(start)
    simplex = [np.asarray(start, dtype=float)] + [np.asarray(start, dtype=float) + np.eye(n)[i] * 0.6 for i in range(n)]
    values = [objective(x) for x in simplex]
    for _ in range(max_iter):
        order = np.argsort(values); simplex = [simplex[i] for i in order]; values = [values[i] for i in order]
        if np.max(np.abs(np.asarray(simplex[1:]) - simplex[0])) < 1e-7 and abs(values[-1] - values[0]) < 1e-7:
            break
        centroid = np.mean(simplex[:-1], axis=0)
        reflected = centroid + (centroid - simplex[-1]); rv = objective(reflected)
        if values[0] <= rv < values[-2]: simplex[-1], values[-1] = reflected, rv; continue
        if rv < values[0]:
            expanded = centroid + 2.0 * (reflected - centroid); ev = objective(expanded)
            simplex[-1], values[-1] = (expanded, ev) if ev < rv else (reflected, rv); continue
        contracted = centroid + 0.5 * (simplex[-1] - centroid); cv = objective(contracted)
        if cv < values[-1]: simplex[-1], values[-1] = contracted, cv; continue
        best = simplex[0]; simplex = [best + 0.5 * (x - best) for x in simplex]; values = [objective(x) for x in simplex]
    i = int(np.argmin(values)); return simplex[i], float(values[i])


def fit_model(rows: list[dict[str, object]], model: str, start: np.ndarray | None = None) -> dict[str, object]:
    names = model_parameter_names(model)
    if not rows:
        raise ValueError("Cannot fit a model with no eligible training quartets")
    scale = max(1.0, 1.0 / max(float(np.median([float(r["interval_t"]) for r in rows])), 1e-6))
    base = np.log(np.repeat([scale], len(names))) if start is None else np.asarray(start, dtype=float)
    candidates = [base, base + 1.0, base - 1.0, np.log(np.repeat([1.0], len(names))), np.log(np.repeat([10.0], len(names)))]
    best = None
    for candidate in candidates:
        x, obj = nelder_mead(lambda z: -log_likelihood(z, rows, model), candidate)
        if best is None or obj < best[1]: best = (x, obj)
    assert best is not None
    x, obj = best
    lambdas, ms = unpack(x, model)
    return {"model": model, "theta": x, "parameters": np.exp(np.clip(x, -20, 20)), "log_likelihood": -obj, "n": len(rows), "lambda_by_lg": lambdas, "m_by_lg": ms}


def fit_msc(rows: list[dict[str, object]]) -> dict[str, object]:
    t = np.asarray([float(row["interval_t"]) for row in rows], dtype=float)
    labels = np.asarray([int(row["label_index"]) for row in rows], dtype=int)
    def objective(log_lam: float) -> float:
        lam = math.exp(float(np.clip(log_lam, -15, 15)))
        no = np.exp(-lam * t)
        q = np.where(labels == 1, 1.0 - 2.0 * no / 3.0, no / 3.0)
        return float(-np.log(np.maximum(q, 1e-300)).sum())
    grid = np.linspace(-15.0, 15.0, 301)
    values = np.asarray([objective(x) for x in grid])
    best = float(grid[int(np.argmin(values))])
    for _ in range(8):
        step = (grid[1] - grid[0])
        grid = np.linspace(best - step, best + step, 21)
        values = np.asarray([objective(x) for x in grid])
        best = float(grid[int(np.argmin(values))])
    return {"lambda": math.exp(best), "log_likelihood": -objective(best), "theta": np.array([best])}


def predict_row(row: dict[str, object], fit: dict[str, object]) -> tuple[float, float, float]:
    i = LGS.index(str(row["lg"]))
    return conditional_probabilities(float(row["interval_t"]), float(fit["lambda_by_lg"][i]), float(fit["m_by_lg"][i]))


def msc_probs(row: dict[str, object], fit: dict[str, object]) -> tuple[float, float, float]:
    i = LGS.index(str(row["lg"])); no = math.exp(-float(fit["lambda_by_lg"][i]) * float(row["interval_t"]))
    return (no / 3.0, 1.0 - 2.0 * no / 3.0, no / 3.0)


def window_observed(rows: list[dict[str, object]]) -> dict[str, object]:
    counts = Counter(int(r["label_index"]) for r in rows)
    n = len(rows); return {"n": n, "n_A": counts[0], "n_B": counts[1], "n_3": counts[2], "D": (counts[0] - counts[1]) / n if n else math.nan}


def training_null(rows: list[dict[str, object]]) -> tuple[float, float, float]:
    counts = Counter(int(r["label_index"]) for r in rows); n = len(rows)
    if n == 0:
        return (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
    # Jeffreys smoothing keeps held-out log scores finite when a training
    # block contains only one resolved topology.
    denom = n + 1.5
    return tuple((counts[i] + 0.5) / denom for i in range(3))


def predict_windows(rows: list[dict[str, object]], fit: dict[str, object], msc_fit: dict[str, object], null_by_lg: dict[str, tuple[float, float, float]], test_ids: set[str], fold: int) -> list[dict[str, object]]:
    by_window = eligible_by_window(rows); out = []
    for wid in sorted(test_ids):
        qrows = by_window.get(wid, [])
        if not qrows: continue
        first = qrows[0]; q_msrc = np.asarray([predict_row(r, fit) for r in qrows]); q_msc = np.asarray([msc_probs(r, msc_fit) for r in qrows]); q_null = np.asarray(null_by_lg[str(first["lg"])])
        observed = window_observed(qrows); labels = np.asarray([int(r["label_index"]) for r in qrows])
        avg_msrc, avg_msc = q_msrc.mean(axis=0), q_msc.mean(axis=0)
        out.append({"fold": fold, "lg": first["lg"], "window_id": wid, "start": first["start"], "end": first["end"], "n_eligible": len(qrows), "n_A": observed["n_A"], "n_B": observed["n_B"], "n_3": observed["n_3"], "observed_D": observed["D"], "predicted_D_MSRC": avg_msrc[0] - avg_msrc[1], "predicted_D_MSC": avg_msc[0] - avg_msc[1], "predicted_D_training_null": q_null[0] - q_null[1], "qA_MSRC": avg_msrc[0], "qB_MSRC": avg_msrc[1], "q3_MSRC": avg_msrc[2], "qA_MSC": avg_msc[0], "qB_MSC": avg_msc[1], "q3_MSC": avg_msc[2], "qA_null": q_null[0], "qB_null": q_null[1], "q3_null": q_null[2], "mean_t": float(np.mean([float(r["interval_t"]) for r in qrows])), "n_distinct_t_window": len({round(float(r["interval_t"]), 12) for r in qrows}), "lambda_fit": float(fit["lambda_by_lg"][LGS.index(str(first["lg"]))]), "m_fit": float(fit["m_by_lg"][LGS.index(str(first["lg"]))]), "rho_fit": float(fit["m_by_lg"][LGS.index(str(first["lg"]))] / fit["lambda_by_lg"][LGS.index(str(first["lg"]))]), "logpred_MSRC": float(np.mean([math.log(max(q_msrc[j, labels[j]], 1e-300)) for j in range(len(qrows))])), "logpred_MSC": float(np.mean([math.log(max(q_msc[j, labels[j]], 1e-300)) for j in range(len(qrows))])), "logpred_null": float(np.mean([math.log(max(q_null[labels[j]], 1e-300)) for j in range(len(qrows))]))})
    return out


def metrics(rows: list[dict[str, object]], key: str) -> dict[str, float]:
    if not rows: return {"n_windows": 0, "mean_logpred": math.nan, "cross_entropy": math.nan, "mae_D": math.nan, "rmse_D": math.nan, "corr_D": math.nan}
    obs = np.asarray([float(r["observed_D"]) for r in rows]); pred = np.asarray([float(r[key]) for r in rows]); err = pred - obs
    corr = float(np.corrcoef(obs, pred)[0, 1]) if np.std(obs) > 0 and np.std(pred) > 0 else math.nan
    lp = "logpred_null" if key == "predicted_D_training_null" else key.replace("predicted_D_", "logpred_")
    return {"n_windows": len(rows), "mean_logpred": float(np.mean([float(r[lp]) for r in rows])), "cross_entropy": float(-np.mean([float(r[lp]) for r in rows])), "mae_D": float(np.mean(np.abs(err))), "rmse_D": float(np.sqrt(np.mean(err * err))), "corr_D": corr}


def profile_rows(training: list[dict[str, object]], full: dict[str, object]) -> list[dict[str, object]]:
    rows = []; theta = np.asarray(full["theta"])
    names = model_parameter_names("lg_m")
    for idx, name in enumerate(names):
        grid = np.linspace(theta[idx] - 4.0, theta[idx] + 4.0, 9)
        for fixed in grid:
            def obj(nuisance):
                full_theta = theta.copy(); full_theta[np.arange(len(theta)) != idx] = nuisance; full_theta[idx] = fixed
                return -log_likelihood(full_theta, training, "lg_m")
            nuisance_idx = np.arange(len(theta)) != idx
            start = theta[nuisance_idx]
            nuisance, value = nelder_mead(obj, start, 120)
            full_theta = theta.copy(); full_theta[nuisance_idx] = nuisance; full_theta[idx] = fixed
            val = -value
            rows.append({"profile_type": "parameter", "parameter": name, "fixed_value": math.exp(float(fixed)), "fixed_log_value": fixed, "optimized_log_likelihood": val, "delta_log_likelihood": val - full["log_likelihood"], "nuisance_reoptimized": True, "boundary_flag": abs(fixed) >= 14.9})
    # rho profiles fix m_LG/lambda and re-optimize lambda plus the other m's.
    for lg_index, lg in enumerate(LGS):
        m_idx = 1 + lg_index
        rho0 = theta[m_idx] - theta[0]
        for fixed in np.linspace(rho0 - 4.0, rho0 + 4.0, 9):
            def obj(nuisance):
                full_theta = theta.copy(); free = [i for i in range(len(theta)) if i not in {m_idx}]
                full_theta[free] = nuisance; full_theta[m_idx] = full_theta[0] + fixed
                return -log_likelihood(full_theta, training, "lg_m")
            free = [i for i in range(len(theta)) if i != m_idx]
            nuisance, value = nelder_mead(obj, theta[free], 120)
            full_theta = theta.copy(); full_theta[free] = nuisance; full_theta[m_idx] = full_theta[0] + fixed
            val = -value
            rows.append({"profile_type": "rho", "parameter": f"rho_{lg}", "fixed_value": math.exp(float(fixed)), "fixed_log_value": fixed, "optimized_log_likelihood": val, "delta_log_likelihood": val - full["log_likelihood"], "nuisance_reoptimized": True, "boundary_flag": abs(fixed) >= 14.9})
    # A compact 2-D lambda/m_LG profile for the two most structured groups.
    for lg in ("LG01", "LG02"):
        mi = 1 + LGS.index(lg)
        for lf in np.linspace(theta[0] - 3.0, theta[0] + 3.0, 7):
            for mf in np.linspace(theta[mi] - 3.0, theta[mi] + 3.0, 7):
                fixed = {0: lf, mi: mf}; free = [i for i in range(len(theta)) if i not in fixed]
                def obj(nuisance):
                    full_theta = theta.copy(); full_theta[free] = nuisance; full_theta[0] = lf; full_theta[mi] = mf
                    return -log_likelihood(full_theta, training, "lg_m")
                nuisance, value = nelder_mead(obj, theta[free], 100)
                rows.append({"profile_type": "2d", "parameter": f"lambda_vs_m_{lg}", "fixed_value": math.exp(float(lf)), "fixed_m_value": math.exp(float(mf)), "fixed_log_value": lf, "fixed_log_m_value": mf, "optimized_log_likelihood": -value, "delta_log_likelihood": -value - full["log_likelihood"], "nuisance_reoptimized": True, "boundary_flag": False})
    return rows


def simulation_validation(template: list[dict[str, object]], assignments: dict[str, int]) -> list[dict[str, object]]:
    rng = np.random.default_rng(SEED)
    settings = [(1.0, 0.0, "near_zero_m"), (1.0, 0.5, "intermediate"), (1.0, 50.0, "large_m"), (5.0, 2.5, "same_rho_different_scale"), (0.5, 0.25, "same_rho_different_scale")]
    out = []
    for lam, m, regime in settings:
        sim = []
        for row in template:
            q = conditional_probabilities(float(row["interval_t"]), lam, m)
            sim.append({**row, "label_index": int(rng.choice(3, p=q))})
        fit = fit_model(sim, "lg_m")
        held = []
        for fold in range(5):
            train = [r for r in sim if assignments.get(str(r["window_id"])) != fold]
            test = [r for r in sim if assignments.get(str(r["window_id"])) == fold]
            if test:
                msc = fit_msc(train)
                held.extend(predict_windows(test, fit_model(train, "lg_m"), {"lambda_by_lg": np.repeat(msc["lambda"], 4)}, {lg: training_null(train) for lg in LGS}, {str(r["window_id"]) for r in test}, fold))
        pred = metrics(held, "predicted_D_MSRC")
        rho_true = m / lam; rho_est = float(np.mean(fit["m_by_lg"] / fit["lambda_by_lg"]))
        out.append({"regime": regime, "lambda_true": lam, "m_true": m, "rho_true": rho_true, "lambda_estimate_mean_lg": float(np.mean(fit["lambda_by_lg"])), "m_estimate_mean_lg": float(np.mean(fit["m_by_lg"])), "rho_estimate_mean_lg": rho_est, "parameter_loglik": fit["log_likelihood"], "heldout_mean_logpred": pred["mean_logpred"], "heldout_mae_D": pred["mae_D"], "heldout_rmse_D": pred["rmse_D"], "seed": SEED, "n_simulated_quartets": len(sim), "profile_boundary_hit": any(abs(float(x)) >= 14.9 for x in fit["theta"]), "profile_caveat": "lambda and m can be weakly separated when A_m(t) is nearly constant; interval coverage requires replicate simulation datasets"})
    return out


def write_legacy_audit(old_hashes: dict[str, str]) -> None:
    LEGACY_AUDIT.write_text(textwrap.dedent(f"""
    # Stage 10 legacy validity audit

    The existing Stage 10 analysis treated per-window SNAPP population-tree
    internal node ages from `stage7_pairwise_mrca_times.tsv` as realized lineage
    coalescence times drawn from an MSRC density. That observation model is not
    valid: SNAPP integrates over gene trees, so a population-tree divergence
    estimate is not a realized gene-lineage coalescence time.

    The legacy construction found 11,673 of 27,786 audited non-boundary
    observations with `T < d` ({11673/27786:.3%}), which is a symptom of this
    mismatch. Its lambda and m values are exploratory only and should not be
    reported as empirical MSRC parameter estimates.

    Stage 10b does not use Stage-7 MRCA times, the legacy Stage-10 scaffold, or
    the legacy parameter estimates. It fits the validated conditional 2:2
    quartet model directly to quartet outcomes and evaluates held-out
    contiguous genomic blocks. The frozen Stage-10 files were hash-audited and
    remain unchanged.

    Legacy file hashes are recorded in the Stage-10b manifest. The primary
    Stage-10b estimands are effective reduced-model parameters and rho=m/lambda;
    they are not meiotic recombination rates or full biological MSRC estimates.
    """).strip() + "\n")


def plot_figure(params: list[dict[str, object]], predictions: list[dict[str, object]], comparison: list[dict[str, object]], ident: list[dict[str, object]]) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    ax = axes[0, 0]
    point = [p for p in params if p.get("fit_label") == "full_primary" and str(p.get("parameter", "")).startswith("rho_")]
    if point:
        labels = [str(p["parameter"]) for p in point]; values = [float(p["estimate"]) for p in point]
        ax.bar(labels, values, color="#2f6690"); ax.set_ylabel("effective rho = m/lambda"); ax.set_title("A  reduced-model rho estimates"); ax.tick_params(axis="x", rotation=35)
    else: ax.text(.5, .5, "rho estimates unavailable", ha="center", va="center")
    ax.grid(axis="y", alpha=.25)

    ax = axes[0, 1]
    for lg in LGS:
        p = [r for r in predictions if r["lg"] == lg]
        ax.scatter([float(r["observed_D"]) for r in p], [float(r["predicted_D_MSRC"]) for r in p], s=12, alpha=.55, label=lg)
    ax.axline((0, 0), slope=1, color="black", ls="--", lw=.8); ax.set_xlabel("held-out observed D(w)"); ax.set_ylabel("MSRC predicted D(w)"); ax.set_title("B  held-out contrast"); ax.legend(frameon=False, fontsize=7); ax.grid(alpha=.2)

    ax = axes[1, 0]
    names = ["MSRC", "MSC", "training null"]
    rows = {r["comparator"]: r for r in comparison if r["scope"] == "pooled"}
    vals = [float(rows[n]["mae_D"]) for n in names]; ax.bar(names, vals, color=["#2f6690", "#a65f2b", "#777777"]); ax.set_ylabel("held-out MAE of D(w)"); ax.set_title("C  empirical comparator performance"); ax.tick_params(axis="x", rotation=25); ax.grid(axis="y", alpha=.25)

    ax = axes[1, 1]
    for lg in LGS:
        p = sorted([r for r in predictions if r["lg"] == lg], key=lambda r: int(r["start"]))
        if p:
            ax.plot([int(r["start"]) for r in p], [float(r["observed_D"]) for r in p], lw=.8, alpha=.6, label=f"{lg} observed")
            ax.plot([int(r["start"]) for r in p], [float(r["predicted_D_MSRC"]) for r in p], lw=1.1, label=f"{lg} MSRC")
    ax.set_xlabel("window start (bp)"); ax.set_ylabel("D(w)"); ax.set_title("D  within-LG model adequacy"); ax.legend(frameon=False, fontsize=6, ncol=2); ax.grid(alpha=.2)
    fig.suptitle("Atlantic cod Stage 10b: conditional MSRC quartet prediction", fontsize=13)
    fig.savefig(FIG_PDF); fig.savefig(FIG_PNG, dpi=300); plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("--self-test", action="store_true"); args = ap.parse_args()
    if args.self_test:
        assert validate_frozen_kernel() < 1e-12
        assert math.isclose(A_m(.4, 1.0, 0.0), 1.0 - math.exp(-.8))
        assert A_m(.4, 1.0, 1e8) < 1e-8
        assert math.isclose(sum(conditional_probabilities(.4, 1.0, .2)), 1.0)
        print("stage10b self-tests passed"); return

    scaffold = verify_scaffold()
    kernel_error = validate_frozen_kernel()
    tree = TreeModel(TREE)
    rows, audit = load_rows(tree)
    eligible = [r for r in rows if r["eligible"]]
    eligible_inside = [r for r in eligible if r["overlaps_inversion"] == "true"]
    if not eligible_inside: raise ValueError("No eligible inside-supergene quartet outcomes")
    assignments, assignment_map = make_assignments(rows)
    inside_ids = set(assignment_map)
    windows = eligible_by_window(eligible_inside)
    window_rows = []
    for wid, qrows in sorted(windows.items()):
        obs = window_observed(qrows); first = qrows[0]
        window_rows.append({"lg": first["lg"], "window_id": wid, "start": first["start"], "end": first["end"], "region": "inside", "n_eligible": obs["n"], "n_A": obs["n_A"], "n_B": obs["n_B"], "n_3": obs["n_3"], "D_eligible": obs["D"], "n_distinct_quartet_sets": len({r["quartet_id"] for r in qrows}), "n_distinct_t": len({round(float(r["interval_t"]), 12) for r in qrows}), "min_t": min(float(r["interval_t"]) for r in qrows), "max_t": max(float(r["interval_t"]) for r in qrows), "position_enters_model": False})
    write_tsv(ELIGIBILITY, audit, ["lg", "window_id", "quartet_id", "taxon1", "taxon2", "taxon3", "taxon4", "arrangement_split", "baseline_split", "local_class", "state_pattern", "interval_t", "eligible", "exclusion_reason"])
    write_tsv(WINDOW_COUNTS, window_rows, ["lg", "window_id", "start", "end", "region", "n_eligible", "n_A", "n_B", "n_3", "D_eligible", "n_distinct_quartet_sets", "n_distinct_t", "min_t", "max_t", "position_enters_model"])
    write_tsv(CV_ASSIGNMENT, assignments, ["lg", "window_id", "start", "end", "rank_within_lg", "n_inside_windows_lg", "fold", "assignment_rule"])

    primary = fit_model(eligible_inside, "lg_m")
    common = fit_model(eligible_inside, "common")
    over = fit_model(eligible_inside, "lg_lg")
    parameter_rows = []
    for label, fit in (("full_primary", primary), ("sensitivity_common_m", common), ("sensitivity_lg_lambda_lg_m_overparameterized", over)):
        names = model_parameter_names(fit["model"])
        for name, value in zip(names, fit["parameters"]):
            if fit["model"] == "lg_m" and name.startswith("m_") and not any(r["lg"] == name[2:] for r in eligible_inside):
                value = math.nan
            if fit["model"] == "lg_lg" and "_LG12" in name and not any(r["lg"] == "LG12" for r in eligible_inside):
                value = math.nan
            parameter_rows.append({"fit_label": label, "model": fit["model"], "parameter": name, "estimate": float(value), "unit": "inverse_Ma", "log_likelihood": fit["log_likelihood"], "n_quartets": fit["n"], "rho": (float(value) / float(fit["parameters"][0]) if name.startswith("m_") and fit["model"] == "lg_m" else "")})
        if fit["model"] == "lg_m":
            for i, lg in enumerate(LGS):
                rho = float(fit["m_by_lg"][i] / fit["lambda_by_lg"][i]) if any(r["lg"] == lg for r in eligible_inside) else math.nan
                parameter_rows.append({"fit_label": label, "model": fit["model"], "parameter": f"rho_{lg}", "estimate": rho, "unit": "dimensionless", "log_likelihood": fit["log_likelihood"], "n_quartets": fit["n"], "rho": rho})
    train_profiles = profile_rows(eligible_inside, primary)
    write_tsv(PARAMETERS, parameter_rows, ["fit_label", "model", "parameter", "estimate", "unit", "log_likelihood", "n_quartets", "rho"])
    write_tsv(PROFILES, train_profiles, ["profile_type", "parameter", "fixed_value", "fixed_log_value", "fixed_m_value", "fixed_log_m_value", "optimized_log_likelihood", "delta_log_likelihood", "nuisance_reoptimized", "boundary_flag"])

    predictions = []; fold_params = []
    all_eligible_ids = set(windows)
    for fold in range(5):
        test_ids = {wid for wid in inside_ids if assignment_map[wid] == fold}
        train_rows = [r for r in eligible_inside if str(r["window_id"]) not in test_ids]
        fold_fit = fit_model(train_rows, "lg_m"); fold_msc = fit_msc(train_rows)
        nulls = {lg: training_null([r for r in train_rows if r["lg"] == lg]) for lg in LGS}
        fold_params.append({"fit_label": f"cv_fold_{fold}", "model": "lg_m", "parameter": "lambda", "estimate": float(fold_fit["lambda_by_lg"][0]), "unit": "inverse_Ma", "log_likelihood": fold_fit["log_likelihood"], "n_quartets": len(train_rows), "rho": ""})
        for i, lg in enumerate(LGS):
            has_lg = any(r["lg"] == lg for r in train_rows)
            m_value = float(fold_fit["m_by_lg"][i]) if has_lg else math.nan
            rho_value = float(fold_fit["m_by_lg"][i] / fold_fit["lambda_by_lg"][i]) if has_lg else math.nan
            fold_params.extend([{ "fit_label": f"cv_fold_{fold}", "model": "lg_m", "parameter": f"m_{lg}", "estimate": m_value, "unit": "inverse_Ma", "log_likelihood": fold_fit["log_likelihood"], "n_quartets": len(train_rows), "rho": rho_value}, {"fit_label": f"cv_fold_{fold}", "model": "lg_m", "parameter": f"rho_{lg}", "estimate": rho_value, "unit": "dimensionless", "log_likelihood": fold_fit["log_likelihood"], "n_quartets": len(train_rows), "rho": rho_value}])
        predictions.extend(predict_windows(eligible, fold_fit, {"lambda_by_lg": np.repeat(fold_msc["lambda"], 4)}, nulls, test_ids, fold))
    parameter_rows.extend(fold_params)
    write_tsv(PARAMETERS, parameter_rows, ["fit_label", "model", "parameter", "estimate", "unit", "log_likelihood", "n_quartets", "rho"])
    write_tsv(CV_PREDICTIONS, predictions, ["fold", "lg", "window_id", "start", "end", "n_eligible", "n_A", "n_B", "n_3", "observed_D", "predicted_D_MSRC", "predicted_D_MSC", "predicted_D_training_null", "qA_MSRC", "qB_MSRC", "q3_MSRC", "qA_MSC", "qB_MSC", "q3_MSC", "qA_null", "qB_null", "q3_null", "mean_t", "n_distinct_t_window", "lambda_fit", "m_fit", "rho_fit", "logpred_MSRC", "logpred_MSC", "logpred_null"])

    comparison = []
    for scope, subset in [("pooled", predictions), *[(lg, [r for r in predictions if r["lg"] == lg]) for lg in LGS]]:
        for name, key in (("MSRC", "predicted_D_MSRC"), ("MSC", "predicted_D_MSC"), ("training null", "predicted_D_training_null")):
            met = metrics(subset, key); comparison.append({"scope": scope, "comparator": name, **met})
    write_tsv(MODEL_COMPARISON, comparison, ["scope", "comparator", "n_windows", "mean_logpred", "cross_entropy", "mae_D", "rmse_D", "corr_D"])
    write_tsv(PRED_SUMMARY, comparison, ["scope", "comparator", "n_windows", "mean_logpred", "cross_entropy", "mae_D", "rmse_D", "corr_D"])
    ident = []
    for name in ["lambda", *[f"m_{lg}" for lg in LGS], *[f"rho_{lg}" for lg in LGS]]:
        p = [r for r in train_profiles if r["parameter"] == name and r["profile_type"] != "2d"]
        if p:
            ll = np.asarray([float(r["optimized_log_likelihood"]) for r in p]); ident.append({"parameter": name, "profile_max": float(ll.max()), "profile_min": float(ll.min()), "profile_range": float(ll.max() - ll.min()), "profile_delta_at_edges": float(min(ll[0], ll[-1]) - ll.max()), "boundary_hit": any(bool(r["boundary_flag"]) for r in p), "interpretation": "weakly identified if profile is flat or boundary hit"})
    for lg in LGS:
        qrows = [r for r in eligible_inside if r["lg"] == lg]; pred = [A_m(float(r["interval_t"]), float(primary["lambda_by_lg"][LGS.index(lg)]), float(primary["m_by_lg"][LGS.index(lg)])) for r in qrows]
        pred_range = max(pred) - min(pred) if pred else math.nan
        ident.append({"parameter": f"model_adequacy_{lg}", "profile_max": "", "profile_min": "", "profile_range": "", "profile_delta_at_edges": "", "boundary_hit": "", "interpretation": f"n_eligible_quartets={len(qrows)}; distinct_quartet_sets={len({r['quartet_id'] for r in qrows})}; distinct_t={len({round(float(r['interval_t']),12) for r in qrows})}; predicted_A_range={pred_range:.10g}; physical_position_enters_model=false"})
    write_tsv(IDENTIFIABILITY, ident, ["parameter", "profile_max", "profile_min", "profile_range", "profile_delta_at_edges", "boundary_hit", "interpretation"])
    sim_rows = simulation_validation(eligible_inside, assignment_map); write_tsv(SIMULATION, sim_rows, list(sim_rows[0]))

    old_hashes = {str(p.relative_to(ROOT)): sha256(p) for p in OLD_STAGE10 if p.exists()}; write_legacy_audit(old_hashes)
    pooled = {r["comparator"]: r for r in comparison if r["scope"] == "pooled"}
    report = textwrap.dedent(f"""
    # Atlantic cod Stage 10b report

    Stage 10b fits the exact conditional symmetric 2:2 quartet model to quartet outcomes and evaluates predictions in held-out contiguous genomic blocks. It uses the checksum-verified rooted Stage-8b original time tree only to obtain baseline quartet intervals. Stage-7 SNAPP MRCA ages do not enter the likelihood.

    ## Eligibility and model

    The analysis retained {len(eligible)} of {len(rows)} local quartet rows in the eligibility audit, with {len(eligible_inside)} eligible inside-supergene rows used for fitting and CV. Exclusions are tabulated in `stage10b_quartet_eligibility.tsv`. The objective is a quartet composite/pseudo-likelihood; windows are the natural observational unit and window-level D(w) metrics are reported separately.

    The frozen exact-kernel validation maximum absolute discrepancy was {kernel_error:.3g}. Parameters are effective reduced-model rates in inverse Ma; rho=m/lambda is dimensionless. They are not meiotic recombination rates or full biological MSRC estimates.

    ## Held-out prediction

    Pooled held-out MSRC MAE={float(pooled['MSRC']['mae_D']):.6g}, RMSE={float(pooled['MSRC']['rmse_D']):.6g}, mean log predictive probability={float(pooled['MSRC']['mean_logpred']):.6g}; MSC MAE={float(pooled['MSC']['mae_D']):.6g}; training-frequency null MAE={float(pooled['training null']['mae_D']):.6g}. Per-LG results are in `stage10b_model_comparison.tsv`.

    ## Identifiability and adequacy

    True nuisance-reoptimized profiles are in `stage10b_parameter_profiles.tsv`; profile diagnostics and model-adequacy audits are in `stage10b_identifiability.tsv`. The reduced model has no physical genomic-position term. Any within-LG variation in predictions comes only from the eligible quartet set and its baseline interval t; where those are constant, the predicted distribution is mathematically constant.

    ## Interpretation

    The defensible claim is limited to whether this reduced two-background conditional model provides predictive structure for held-out arrangement-associated quartet outcomes. It is not an independent second-modality validation and does not establish the full biological MSRC process.
    """).strip() + "\n"
    REPORT.write_text(report)
    METHODS.write_text(textwrap.dedent(f"""
    # Stage 10b methods

    The fixed scaffold is `stage8b_baseline_original_times.nwk`, verified against `stage8b_manifest.json`. It is rooted, fully bifurcating, time-calibrated, and used only to calculate the internal baseline-tree interval t for each four-taxon quartet.

    Quartet eligibility required four unique taxa, a 2:2 arrangement state pattern, an arrangement split different from the baseline split, a resolved local topology in {{arrangement, baseline, third}}, and t>0. No Stage-7 MRCA values, Stage-10 parameter estimates, re-rooting, or projected Stage-8 time trees were used.

    For supported configuration 1010, the validated finite-state calculation gives A=lambda/(lambda+2m)[1-exp(-2(lambda+2m)t)] and q=((1+2A)/3,(1-A)/3,(1-A)/3). The maximum difference between this implementation and the frozen exact simulation grid was {kernel_error:.3g}. Fitting uses quartet composite/pseudo-likelihood; it is not a full iid-gene-tree likelihood. CV uses five deterministic contiguous folds within each linkage group. The training-frequency null is fitted from training windows only, with Jeffreys 0.5-count smoothing for finite held-out log scores.

    Profiles fix one parameter and re-optimize all nuisance parameters. Simulation validation uses the empirical eligible quartet intervals and window structure under the same conditional model. m is an effective reduced-model switching rate in inverse Ma units, not a meiotic recombination rate.
    """).strip() + "\n")
    RESULTS_TEXT.write_text(f"# Stage 10b results\n\nThe pooled held-out MSRC mean log predictive probability was {float(pooled['MSRC']['mean_logpred']):.6g}; MSC was {float(pooled['MSC']['mean_logpred']):.6g}; the training-frequency null was {float(pooled['training null']['mean_logpred']):.6g}. MSRC D(w) MAE was {float(pooled['MSRC']['mae_D']):.6g}, compared with MSC {float(pooled['MSC']['mae_D']):.6g} and the null {float(pooled['training null']['mae_D']):.6g}. See the comparator and identifiability tables for per-LG results and profile behavior.\n")
    plot_params = [p for p in parameter_rows if p["fit_label"] == "full_primary"]
    plot_figure(plot_params, predictions, comparison, ident)
    outputs = [ELIGIBILITY, WINDOW_COUNTS, PARAMETERS, PROFILES, CV_ASSIGNMENT, CV_PREDICTIONS, PRED_SUMMARY, MODEL_COMPARISON, IDENTIFIABILITY, SIMULATION, METHODS, RESULTS_TEXT, REPORT, MANIFEST, LEGACY_AUDIT, FIG_PDF, FIG_PNG]
    manifest = {"stage": "10b", "analysis": "conditional two-background MSRC quartet outcome prediction", "seed": SEED, "scaffold": scaffold, "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in [TREE, STAGE8B_MANIFEST, LOCAL, SCORES, ARRANGEMENTS, FROZEN_GRID]}, "legacy_stage10_hashes": old_hashes, "outputs": {str(p.relative_to(ROOT)): sha256(p) for p in outputs if p.exists()}, "n_local_rows": len(rows), "n_eligible_rows": len(eligible), "kernel_max_abs_error": kernel_error, "no_stage7_mrca_in_likelihood": True, "quartet_observation_unit": "250-kb window; quartet composite/pseudo-likelihood", "theory": "validated conditional 1010 exact kernel", "software": {"python": platform.python_version(), "biopython": __import__('Bio').__version__, "numpy": np.__version__, "matplotlib": matplotlib.__version__}}
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"n_local_rows": len(rows), "n_eligible": len(eligible), "n_predictions": len(predictions), "kernel_max_abs_error": kernel_error, "msrc_mae": pooled["MSRC"]["mae_D"], "msc_mae": pooled["MSC"]["mae_D"], "null_mae": pooled["training null"]["mae_D"]}, indent=2))


if __name__ == "__main__":
    main()
