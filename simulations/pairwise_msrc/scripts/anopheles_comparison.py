#!/usr/bin/env python3
"""Matched-site Anopheles summaries and illustrative MSRC parameter surfaces."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pairwise_msrc import Segment, expected_pair_times, predicted_statistics  # noqa: E402

TRACK = ROOT / "empirical/anopheles_2la/analysis/joint_genomic_tracks/results/anopheles_2la_dxy_fst_windows.tsv"
OUT = ROOT / "simulations/pairwise_msrc/results"
FIG = ROOT / "simulations/pairwise_msrc/figures"
INVERSION_START, INVERSION_END = 20_524_058, 42_165_532


def write_frame(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(OUT / name, sep="\t", index=False)


def block_bootstrap(frame: pd.DataFrame, metric: str, reps: int = 10_000,
                    seed: int = 84291) -> tuple[float, float]:
    """Ratio-of-sums estimate and resampling CI over 500-kb physical blocks."""
    blocks = frame.groupby("block", sort=False)[["dxy_sum", "delta_sum", "pi_sum"]].sum().to_numpy()
    # Each row carries sufficient totals and callable length; resample whole blocks.
    vals = np.empty(reps)
    rng = np.random.default_rng(seed)
    length = frame.groupby("block", sort=False)["callable_sites"].sum().to_numpy()
    n = len(blocks)
    if metric == "Delta":
        numerator = frame.groupby("block", sort=False)["delta_sum"].sum().to_numpy()
        for i in range(reps):
            ids = rng.integers(n, size=n)
            vals[i] = numerator[ids].sum() / length[ids].sum()
    else:
        pi_sum = frame.groupby("block", sort=False)["pi_sum"].sum().to_numpy()
        dxy = frame.groupby("block", sort=False)["dxy_sum"].sum().to_numpy()
        for i in range(reps):
            ids = rng.integers(n, size=n)
            vals[i] = 2 * dxy[ids].sum() / pi_sum[ids].sum()
    return tuple(np.quantile(vals, [0.025, 0.975]))


def empirical() -> pd.DataFrame:
    x = pd.read_csv(TRACK, sep="\t")
    # Exclude boundary-overlap windows exactly as the frozen regional comparison does.
    x = x.loc[(x.region_class != "boundary") & x.comparison.isin(
        ["primary_within_gambiae", "primary_within_coluzzii"]
    )].copy()
    x["callable_sites"] = x.callable_sites.astype(float)
    x["hudson_fst_sites"] = x.hudson_fst_sites.astype(float)
    if not np.array_equal(x.callable_sites.to_numpy(), x.hudson_fst_sites.to_numpy()):
        raise ValueError("dXY and Hudson FST denominators differ")
    x["delta_per_site"] = x.dxy_per_callable_site - (x.pi_X + x.pi_Y) / 2
    x["R_window"] = 2 * x.dxy_per_callable_site / (x.pi_X + x.pi_Y)
    x["delta_sum"] = x.dxy_sum - ((x.pi_X + x.pi_Y) / 2) * x.callable_sites
    x["pi_sum"] = (x.pi_X + x.pi_Y) * x.callable_sites
    x["block"] = (x.start.astype(int) - 1) // 500_000
    x["midpoint_mb"] = (x.start + x.end) / 2e6
    write_frame(x, "anopheles_matched_delta_R_windows.tsv")

    rows = []
    for comp, dat in x.groupby("comparison", sort=True):
        for reg, part in dat.groupby("region_class", sort=True):
            n = part.callable_sites.sum()
            dxy = part.dxy_sum.sum() / n
            pi_x = (part.pi_X * part.callable_sites).sum() / n
            pi_y = (part.pi_Y * part.callable_sites).sum() / n
            delta = dxy - (pi_x + pi_y) / 2
            ratio = 2 * dxy / (pi_x + pi_y)
            d_ci = block_bootstrap(part, "Delta", seed=13 + len(rows))
            r_ci = block_bootstrap(part, "R", seed=69 + len(rows))
            rows.append({"comparison": comp, "region": reg, "n_windows": len(part),
                         "n_500kb_blocks": part.block.nunique(), "callable_bases": int(n),
                         "pi_X": pi_x, "pi_Y": pi_y, "dXY": dxy,
                         "Delta_AB": delta, "Delta_block_bootstrap_low": d_ci[0],
                         "Delta_block_bootstrap_high": d_ci[1], "R_AB": ratio,
                         "R_block_bootstrap_low": r_ci[0], "R_block_bootstrap_high": r_ci[1],
                         "pooled_Hudson_FST": (part.hudson_fst * part.dxy_sum).sum()
                         / part.dxy_sum.sum()})
    out = pd.DataFrame(rows)
    write_frame(out, "anopheles_matched_delta_R_inside_outside.tsv")

    contrasts = []
    rng = np.random.default_rng(70129)
    for comp, dat in x.groupby("comparison", sort=True):
        grouped = {reg: part.groupby("block", sort=False)[
            ["delta_sum", "dxy_sum", "pi_sum", "callable_sites"]].sum().to_numpy()
                   for reg, part in ((r, dat.loc[dat.region_class == r])
                                     for r in ("inside", "outside"))}
        region_replicates = {}
        for reg, a in grouped.items():
            sampled = a[rng.integers(len(a), size=(10_000, len(a)))].sum(axis=1)
            region_replicates[reg] = {
                "Delta": sampled[:, 0] / sampled[:, 3],
                "R": 2 * sampled[:, 1] / sampled[:, 2],
            }
        item = {"comparison": comp,
                "uncertainty": "independent 500-kb physical block bootstrap within region"}
        for metric in ("Delta", "R"):
            contrast = region_replicates["inside"][metric] - region_replicates["outside"][metric]
            item[f"{metric}_inside_minus_outside"] = contrast.mean()
            item[f"{metric}_contrast_block_bootstrap_low"] = np.quantile(contrast, .025)
            item[f"{metric}_contrast_block_bootstrap_high"] = np.quantile(contrast, .975)
        contrasts.append(item)
    write_frame(pd.DataFrame(contrasts), "anopheles_matched_delta_R_contrasts.tsv")

    # Join existing A(w), B-like, and quartet tracks by frozen genomic window ID.
    gene = pd.read_csv(ROOT / "empirical/anopheles_2la/results/stage12_distance/stage12_window_distance.tsv",
                       sep="\t")
    quartet = pd.read_csv(ROOT / "empirical/anopheles_2la/results/stage2r_quartet_topology_signal.tsv",
                          sep="\t")
    joined = x.merge(gene[["window_id", "A_distance", "mean_pair_distance", "tree_diameter"]],
                     on="window_id", how="left").merge(
                         quartet[["window_id", "q_species", "q_arrangement", "D"]],
                         on="window_id", how="left")
    write_frame(joined, "anopheles_pairwise_and_genealogy_join.tsv")
    corr_rows = []
    for (comp, region), dat in joined.groupby(["comparison", "region_class"], sort=True):
        for predictor in ("A_distance", "mean_pair_distance", "tree_diameter",
                          "q_arrangement", "D"):
            z = dat[["delta_per_site", "R_window", predictor, "block"]].dropna()
            if len(z) < 5:
                continue
            block_ids = z.block.unique()
            rng = np.random.default_rng(1903 + len(corr_rows))
            boot = {"Delta": [], "R": []}
            raw_predictor = z[predictor].to_numpy()
            raw_responses = {"Delta": z.delta_per_site.to_numpy(),
                             "R": z.R_window.to_numpy()}
            for _ in range(1000):
                chosen = rng.choice(block_ids, size=len(block_ids), replace=True)
                indices = np.concatenate([np.flatnonzero(z.block.to_numpy() == b) for b in chosen])
                for response in ("Delta", "R"):
                    a = rankdata(raw_responses[response][indices])
                    b = rankdata(raw_predictor[indices])
                    if np.std(a) and np.std(b):
                        boot[response].append(np.corrcoef(a, b)[0, 1])
            row = {"comparison": comp, "region": region,
                   "predictor": predictor, "n_windows": len(z),
                   "uncertainty": "500-kb block bootstrap; exploratory associations"}
            for response in ("Delta", "R"):
                value = spearmanr(z["delta_per_site" if response == "Delta" else "R_window"],
                                  z[predictor]).statistic
                lo, hi = np.nanquantile(boot[response], [.025, .975])
                row[f"spearman_{response}"] = value
                row[f"{response}_block_bootstrap_low"] = lo
                row[f"{response}_block_bootstrap_high"] = hi
            corr_rows.append(row)
    write_frame(pd.DataFrame(corr_rows), "anopheles_genealogy_track_correlations.tsv")

    # Figure: four matched-site tracks; vertical lines are published inversion bounds.
    fig, axes = plt.subplots(4, 1, figsize=(12, 9), sharex=True, constrained_layout=True)
    for comp, color in [("primary_within_gambiae", "#2563a6"),
                        ("primary_within_coluzzii", "#c45a2a")]:
        dat = x.loc[x.comparison == comp].sort_values("midpoint_mb")
        label = "Gambiae" if "gambiae" in comp else "Coluzzii"
        for ax, col in zip(axes, ["delta_per_site", "R_window", "pi_X", "pi_Y"]):
            ax.plot(dat.midpoint_mb, dat[col], color=color, lw=.65, alpha=.8,
                    label=label)
    ylabels = [r"$\Delta_{AB}$ per callable bp", r"$R_{AB}$", r"$\pi_X$ per callable bp",
               r"$\pi_Y$ per callable bp"]
    for ax, label in zip(axes, ylabels):
        for pos in (INVERSION_START / 1e6, INVERSION_END / 1e6):
            ax.axvline(pos, color="black", ls="--", lw=1)
        ax.set_ylabel(label)
        ax.grid(alpha=.2)
    axes[0].legend(ncol=2, frameon=False)
    axes[-1].set_xlabel("AgamP4 chromosome 2L position (Mb; 1-based window midpoints)")
    fig.suptitle("Anopheles 2La: matched callable-site population statistics")
    fig.savefig(FIG / "anopheles_matched_statistics.pdf")
    fig.savefig(FIG / "anopheles_matched_statistics.png", dpi=300)
    plt.close(fig)
    return out


def theoretical() -> pd.DataFrame:
    rows = []
    p_grid = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95]
    gamma_grid = [0, .03, .1, .3, 1, 3, 10, 30, 100]
    tau_grid = [.1, .25, .5, 1, 2, 4, 8]
    ne = 10_000.0
    for ne_ratio in (.5, 1, 2, 5, 10):
        ne_in = ne_ratio * ne
        for p in p_grid:
            for gamma in gamma_grid:
                for tau_c in tau_grid:
                    rho = gamma / (2 * ne_in)
                    t = expected_pair_times([Segment(tau_c * 2 * ne_in, p, ne_in, rho)],
                                            ancestral_ne=ne)
                    s = predicted_statistics(t, 1.0)
                    rows.append({"Ne_in_over_Ne_out": ne_ratio, "p_B": p,
                                 "gamma_2Ne_in_rho": gamma,
                                 "tau_over_2Ne_in": tau_c,
                                 **{f"T_{k}": v for k, v in t.items()},
                                 "Delta_over_2muNe_out": s["E_Delta_AB"] / (2 * ne),
                                 "R_ratio_of_expectations": s["R_ratio_of_expectations"],
                                 "FST_ratio_of_expectations": s["FST_ratio_of_expectations"]})
    grid = pd.DataFrame(rows)
    write_frame(grid, "theoretical_parameter_grid.tsv")
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for p in (.1, .5, .9):
        sub = grid.loc[(grid.Ne_in_over_Ne_out == 1) & (grid.p_B == p)
                       & (grid.tau_over_2Ne_in == 2)].sort_values("gamma_2Ne_in_rho")
        lab = f"p={p:g}"
        axes[0, 0].plot(sub.gamma_2Ne_in_rho, sub.T_AA / (2 * ne), marker="o", label=lab)
        axes[0, 0].plot(sub.gamma_2Ne_in_rho, sub.T_AB / (2 * ne), marker="s", ls="--")
        axes[0, 0].plot(sub.gamma_2Ne_in_rho, sub.T_BB / (2 * ne), marker="^", ls=":")
        axes[0, 1].plot(sub.gamma_2Ne_in_rho, sub.R_ratio_of_expectations, marker="o", label=lab)
        axes[1, 0].plot(sub.gamma_2Ne_in_rho, sub.Delta_over_2muNe_out, marker="o", label=lab)
    axes[0, 0].set_ylabel("Pairwise E[T] / 2Ne")
    axes[0, 0].set_title("Solid AA, dashed AB, dotted BB; tau=4Ne")
    axes[0, 1].set_ylabel("R (ratio of expected dXY and pi)")
    axes[0, 1].set_title("R declines across the tested switching range")
    axes[1, 0].set_ylabel(r"$\Delta/(2\mu N_e)$")
    axes[1, 1].axis("off")
    axes[0, 1].legend(frameon=False, fontsize=8)
    axes[1, 0].legend(frameon=False, fontsize=8)
    for ax in axes.flat[:3]:
        ax.set_xscale("symlog", linthresh=.03)
        ax.set_xlabel(r"$\gamma=2N_{e,in}\rho$")
        ax.grid(alpha=.25)
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "pairwise_msrc_parameter_sensitivity.pdf")
    fig.savefig(FIG / "pairwise_msrc_parameter_sensitivity.png", dpi=300)
    plt.close(fig)
    return grid


def compare_observed_to_grid(observed: pd.DataFrame, grid: pd.DataFrame) -> pd.DataFrame:
    """Moment-distance sensitivity fit; retain all near-best ridges."""
    fits = []
    for comp, dat in observed.groupby("comparison", sort=True):
        d = dat.set_index("region")
        outside = d.loc["outside"]
        baseline = (outside.pi_X + outside.pi_Y + outside.dXY) / 3
        inside = d.loc["inside"]
        target = np.asarray([inside.pi_X, inside.pi_Y, inside.dXY]) / baseline
        pred = np.column_stack([grid.T_AA / (2 * 10_000), grid.T_BB / (2 * 10_000),
                                grid.T_AB / (2 * 10_000)])
        score = np.square(pred - target).sum(axis=1)
        best = grid.assign(
            comparison=comp, outside_reference_mean=baseline,
            observed_pi_X_over_outside=target[0],
            observed_pi_Y_over_outside=target[1],
            observed_dXY_over_outside=target[2],
            predicted_pi_A_over_outside=pred[:, 0],
            predicted_pi_B_over_outside=pred[:, 1],
            predicted_dXY_over_outside=pred[:, 2], equal_weight_squared_error=score,
            fitting_note="descriptive equal-weight moment distance; no likelihood or demographic prior",
        ).sort_values("equal_weight_squared_error").head(25)
        fits.append(best)
    result = pd.concat(fits, ignore_index=True)
    write_frame(result, "anopheles_grid_moment_fit_top25.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True, constrained_layout=True)
    for ax, (comp, dat) in zip(axes, observed.groupby("comparison", sort=True)):
        best = result.loc[result.comparison == comp].iloc[0]
        d = dat.set_index("region")
        outside = d.loc["outside"]
        baseline = (outside.pi_X + outside.pi_Y + outside.dXY) / 3
        target = [d.loc["inside", "pi_X"] / baseline,
                  d.loc["inside", "pi_Y"] / baseline,
                  d.loc["inside", "dXY"] / baseline]
        prediction = [best.predicted_pi_A_over_outside,
                      best.predicted_pi_B_over_outside,
                      best.predicted_dXY_over_outside]
        pos = np.arange(3)
        ax.scatter(pos - .08, target, label="observed 2La / outside reference", s=42)
        ax.scatter(pos + .08, prediction, marker="D", label="best grid point", s=36)
        ax.set_xticks(pos, [r"$\pi_X$", r"$\pi_Y$", r"$d_{XY}$"])
        ax.set_title("Gambiae" if "gambiae" in comp else "Coluzzii")
        ax.grid(axis="y", alpha=.25)
        ax.set_ylabel("inside statistic / mean outside statistic")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Anopheles moment-distance sensitivity fit (not parameter inference)")
    fig.savefig(FIG / "anopheles_observed_vs_grid.pdf")
    fig.savefig(FIG / "anopheles_observed_vs_grid.png", dpi=300)
    plt.close(fig)
    return result


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    summaries = empirical()
    grid = theoretical()
    fit = compare_observed_to_grid(summaries, grid)
    print("best grid points:")
    print(fit.groupby("comparison", sort=True).head(1)[[
        "comparison", "Ne_in_over_Ne_out", "p_B", "gamma_2Ne_in_rho",
        "tau_over_2Ne_in", "equal_weight_squared_error"]].to_string(index=False))
    print(summaries.to_string(index=False))
    # A compact nonmonotonicity audit across p and duration.
    changing = 0
    for _, sub in grid.groupby(["Ne_in_over_Ne_out", "p_B", "tau_over_2Ne_in"]):
        r = sub.sort_values("gamma_2Ne_in_rho").R_ratio_of_expectations.to_numpy()
        d = np.diff(r)
        changing += bool(np.any(d > 1e-9) and np.any(d < -1e-9))
    print(f"p,tau combinations with nonmonotone R across gamma: {changing}")
