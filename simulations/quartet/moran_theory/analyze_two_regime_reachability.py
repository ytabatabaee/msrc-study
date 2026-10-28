"""Bounded exact search of unconditional two-regime quartet reachability."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from .full_species_tree import full_species_tree_quartet_probabilities
from .two_regime_species_tree import (
    RegimeBranch, TwoRegimeModel, balanced_two_regime_quartet_probabilities,
    unbalanced_two_regime_quartet_probabilities,
)


def _metrics(q):
    q = np.asarray(q)
    return float(max(q[1], q[2]) - q[0]), float(abs(q[1] - q[2]))


def _evaluate(rng, shape, index):
    p_low = float(rng.uniform(.01, .49)); p_high = float(rng.uniform(.51, .99))
    rho = float(10 ** rng.uniform(-2, .4)); a = float(10 ** rng.uniform(-2.5, .5)); b = float(10 ** rng.uniform(-2.5, .5))
    M = int(rng.choice([4, 6, 10, 20])); pi = np.array([rng.uniform(.05, .95), 0.0]); pi[1] = 1 - pi[0]
    model = TwoRegimeModel(p_low, p_high, rho)
    def branch(t): return RegimeBranch(float(t), a, b, M)
    root = float(rng.uniform(0, 2)); left = float(rng.uniform(0, 1.5)); right = float(rng.uniform(0, 1.5))
    terminal = tuple(branch(rng.uniform(0, .8)) for _ in range(4))
    # Random stochastic inheritance is part of the model class.  Most draws
    # are near identity to keep the search biologically interpretable.
    def kernel():
        if rng.random() < .5:
            return np.eye(2)
        x = rng.uniform(.05, .95, size=2)
        return np.array([[x[0], 1 - x[0]], [x[1], 1 - x[1]]])
    if shape == "balanced":
        kernels = [kernel() for _ in range(6)]
        result = balanced_two_regime_quartet_probabilities(
            model=model, root_branch=branch(root), left_pair_branch=branch(left), right_pair_branch=branch(right),
            terminal_branches=terminal, initial_regime=pi,
            H_root_left=kernels[0], H_root_right=kernels[1], H_left_1=kernels[2], H_left_2=kernels[3], H_right_3=kernels[4], H_right_4=kernels[5],
        )
    else:
        kernels = [kernel() for _ in range(6)]
        result = unbalanced_two_regime_quartet_probabilities(
            model=model, root_branch=branch(root), intermediate_branch=branch(left), sister_branch=branch(right),
            terminal_branches=terminal, initial_regime=pi,
            H_intermediate_pair=kernels[0], H_intermediate_3=kernels[1], H_root_intermediate=kernels[2], H_root_4=kernels[3], H_pair_1=kernels[4], H_pair_2=kernels[5],
        )
    anomaly, asymmetry = _metrics(result.quartet_probabilities)
    return {"shape": shape, "index": index, "M": M, "p_low": p_low, "p_high": p_high,
            "rho": rho, "a": a, "b": b, "pi_low": pi[0], "root": root, "branch1": left,
            "branch2": right, "terminal_lengths": ",".join(str(x.duration) for x in terminal),
            "q1": result.quartet_probabilities[0], "q2": result.quartet_probabilities[1],
            "q3": result.quartet_probabilities[2], "anomaly": anomaly, "asymmetry": asymmetry,
            "normalization_error": float(result.quartet_probabilities.sum() - 1),
            "inheritance": json.dumps([k.tolist() for k in kernels], separators=(",", ":"))}


def _from_row(row, with_joint=False):
    model = TwoRegimeModel(float(row["p_low"]), float(row["p_high"]), float(row["rho"]))
    a, b, M = float(row["a"]), float(row["b"]), int(row["M"])
    branch = lambda t: RegimeBranch(float(t), a, b, M)
    terminal = tuple(branch(x) for x in map(float, row["terminal_lengths"].split(",")))
    kernels = [np.asarray(x, dtype=float) for x in json.loads(row["inheritance"])]
    pi = np.array([float(row["pi_low"]), 1.0 - float(row["pi_low"])])
    if row["shape"] == "balanced":
        return balanced_two_regime_quartet_probabilities(
            model=model, root_branch=branch(row["root"]), left_pair_branch=branch(row["branch1"]),
            right_pair_branch=branch(row["branch2"]), terminal_branches=terminal, initial_regime=pi,
            H_root_left=kernels[0], H_root_right=kernels[1], H_left_1=kernels[2], H_left_2=kernels[3],
            H_right_3=kernels[4], H_right_4=kernels[5], with_joint=with_joint)
    return unbalanced_two_regime_quartet_probabilities(
        model=model, root_branch=branch(row["root"]), intermediate_branch=branch(row["branch1"]),
        sister_branch=branch(row["branch2"]), terminal_branches=terminal, initial_regime=pi,
        H_intermediate_pair=kernels[0], H_intermediate_3=kernels[1], H_root_intermediate=kernels[2],
        H_root_4=kernels[3], H_pair_1=kernels[4], H_pair_2=kernels[5], with_joint=with_joint)


def _dense_recheck(row):
    import scipy.linalg
    from . import two_regime_species_tree as impl
    old = impl.expm_multiply
    impl.expm_multiply = lambda A, v: scipy.linalg.expm(A.toarray() if hasattr(A, "toarray") else np.asarray(A)) @ v
    try:
        stored = np.array([float(row[f"q{i}"]) for i in (1, 2, 3)])
        return float(np.max(np.abs(_from_row(row).quartet_probabilities - stored)))
    finally:
        impl.expm_multiply = old


def jacobian_rank():
    model = TwoRegimeModel(.2, .8, .17)
    base = dict(root_branch=RegimeBranch(.5, .2, .4, 10), left_pair_branch=RegimeBranch(.3, .2, .4, 10), right_pair_branch=RegimeBranch(.1, .2, .4, 10), terminal_branches=tuple(RegimeBranch(x, .2, .4, 10) for x in (.2, .4, .1, .3)), initial_regime=np.array([.6, .4]))
    def q(theta):
        m = TwoRegimeModel(theta[0], theta[1], theta[2])
        return balanced_two_regime_quartet_probabilities(model=m, **base).quartet_probabilities[:2]
    theta = np.array([.2, .8, .17]); steps = np.array([1e-4, 1e-4, 1e-4]); J = np.column_stack([(q(theta + np.eye(3)[j] * steps[j]) - q(theta - np.eye(3)[j] * steps[j])) / (2 * steps[j]) for j in range(3)])
    return int(np.linalg.matrix_rank(J, tol=1e-7)), np.linalg.svd(J, compute_uv=False)


def run(outdir: Path, n: int, seed: int):
    outdir.mkdir(parents=True, exist_ok=True); rng = np.random.default_rng(seed)
    rows = [_evaluate(rng, shape, i) for shape in ("balanced", "unbalanced") for i in range(n)]
    fields = list(rows[0])
    with (outdir / "search.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    best_anomaly = max(rows, key=lambda r: r["anomaly"]); best_asym = max(rows, key=lambda r: r["asymmetry"])
    rank, singular = jacobian_rank()
    dense_error = _dense_recheck(best_anomaly)
    with (outdir / "summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t"); writer.writerow(("metric", "value")); writer.writerow(("n_per_shape", n)); writer.writerow(("best_anomaly", best_anomaly)); writer.writerow(("best_asymmetry", best_asym)); writer.writerow(("jacobian_rank", rank)); writer.writerow(("jacobian_singular_values", ",".join(map(str, singular)))); writer.writerow(("best_anomaly_dense_recheck_error", dense_error))
    _write_model_comparison(outdir)
    _plot(outdir, rows, best_anomaly)
    print("best anomaly", best_anomaly); print("best asymmetry", best_asym); print("jacobian rank", rank, singular)


def _xy(q):
    # top q1, lower-left q2, lower-right q3
    return .5 * (2 * q[0] + q[2]), (np.sqrt(3) / 2) * q[0]


def _write_model_comparison(outdir):
    model = TwoRegimeModel(.2, .8, .17)
    b = lambda t: RegimeBranch(t, .2, .4, 10)
    two = balanced_two_regime_quartet_probabilities(
        model=model, root_branch=b(.5), left_pair_branch=b(.3), right_pair_branch=b(.1),
        terminal_branches=tuple(b(x) for x in (.2, .4, .1, .3)), initial_regime=np.array([.5, .5]))
    moran = full_species_tree_quartet_probabilities(
        M_root=10, M_left=10, M_right=10, M_terminal=(10,) * 4, rho=.17,
        root_length=.5, sister_lengths=(.3, .1), terminal_lengths=(.2, .4, .1, .3), pi_old=np.eye(11)[5])
    with (outdir / "model_comparison.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t"); writer.writerow(("model", "q1", "q2", "q3")); writer.writerow(("two_regime", *two.quartet_probabilities)); writer.writerow(("finite_M_moran", *moran.quartet_probabilities))


def _plot(outdir, rows, best_row):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(6.6, 6.0), constrained_layout=True)
    colors = {"balanced": "#0072B2", "unbalanced": "#D55E00"}
    for ax, shape in zip(axes[0], ("balanced", "unbalanced")):
        tri = np.array([_xy((1, 0, 0)), _xy((0, 1, 0)), _xy((0, 0, 1)), _xy((1, 0, 0))])
        ax.plot(tri[:, 0], tri[:, 1], color="#555555", lw=.8)
        pts = [(_xy((r["q1"], r["q2"], r["q3"]))) for r in rows if r["shape"] == shape]
        if pts: ax.scatter(*np.array(pts).T, s=12, alpha=.45, color=colors[shape], edgecolors="none")
        ax.plot([1, .25], [np.sqrt(3)/2, 0], ls="--", color="#999999", lw=.7)
        ax.set_aspect("equal"); ax.set_title(shape); ax.axis("off")
        ax.text(.5, np.sqrt(3)/2 + .04, "q1", ha="center", fontsize=8); ax.text(-.04, -.04, "q2", ha="right", fontsize=8); ax.text(1.04, -.04, "q3", ha="left", fontsize=8)
    candidate = _from_row(best_row, with_joint=True)
    ax = axes[1, 0]
    tri = np.array([_xy((1, 0, 0)), _xy((0, 1, 0)), _xy((0, 0, 1)), _xy((1, 0, 0))])
    ax.plot(tri[:, 0], tri[:, 1], color="#555555", lw=.8)
    for i, q in enumerate(candidate.conditional_quartet_probabilities):
        ax.scatter(*_xy(q), s=25 + 250 * candidate.pattern_probabilities[i], alpha=.65, color=plt.cm.viridis(i / 15), edgecolors="none")
    ax.plot([1, .25], [np.sqrt(3)/2, 0], ls="--", color="#999999", lw=.7)
    ax.set_aspect("equal"); ax.set_title("conditional: anomaly witness"); ax.axis("off")
    ax.text(.5, np.sqrt(3)/2 + .04, "q1", ha="center", fontsize=8); ax.text(-.04, -.04, "q2", ha="right", fontsize=8); ax.text(1.04, -.04, "q3", ha="left", fontsize=8)
    ax = axes[1, 1]
    ax.plot(tri[:, 0], tri[:, 1], color="#555555", lw=.8)
    for shape in ("balanced", "unbalanced"):
        pts = [_xy((r["q1"], r["q2"], r["q3"])) for r in rows if r["shape"] == shape]
        ax.scatter(*np.array(pts).T, s=12, alpha=.35, color=colors[shape], edgecolors="none", label=shape)
    ax.scatter(*_xy(tuple(float(best_row[f"q{i}"]) for i in (1, 2, 3))), marker="*", s=70, color="black", label="anomaly witness")
    ax.plot([1, .25], [np.sqrt(3)/2, 0], ls="--", color="#999999", lw=.7)
    ax.set_aspect("equal"); ax.set_title("unconditional search"); ax.axis("off"); ax.legend(frameon=False, fontsize=7, loc="lower left")
    ax.text(.5, np.sqrt(3)/2 + .04, "q1", ha="center", fontsize=8); ax.text(-.04, -.04, "q2", ha="right", fontsize=8); ax.text(1.04, -.04, "q3", ha="left", fontsize=8)
    fig.savefig(outdir / "two_regime_reachability_simplex.pdf")
    fig.savefig(outdir / "two_regime_reachability_simplex.svg")
    fig.savefig(outdir / "two_regime_reachability_simplex.png", dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results" / "two_regime_reachability"); parser.add_argument("--n-per-shape", type=int, default=30); parser.add_argument("--seed", type=int, default=20261028); args = parser.parse_args(); run(args.output, args.n_per_shape, args.seed)


if __name__ == "__main__": main()
