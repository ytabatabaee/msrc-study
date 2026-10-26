"""Validate arrangement-conditional pectinate quartet probabilities."""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import chi2

from .full_species_tree import full_species_tree_quartet_probabilities
from .moran_quartet_operator import PATTERNS
from .unbalanced_species_tree import unbalanced_species_tree_conditional_quartet_probabilities


def condition_grid() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    def add(label: str, *, p0=.5, rho=.17, S=.3, I=.2, R=.5,
            M_root=10, M_intermediate=10, M_sister=10,
            M_terminal=(10, 10, 10, 10), selection=0.0):
        row = dict(condition=label, p0=p0, rho=rho, sister_length=S,
                   intermediate_length=I, root_length=R, M_root=M_root,
                   M_intermediate=M_intermediate, M_sister=M_sister,
                   M_terminal=tuple(M_terminal), selection=selection)
        row["k0"] = int(round(M_root * p0))
        rows.append(row)

    add("base")
    add("rho_0", rho=0.0)
    add("rho_high", rho=1.0)
    add("p0_low", p0=.1)
    add("p0_high", p0=.9)
    add("selection_positive", selection=.25)
    add("selection_negative", selection=-.25)
    add("sister_long", S=.8, I=.1, R=.3)
    add("intermediate_long", S=.1, I=.8, R=.3)
    add("zero_branches", S=0.0, I=0.0, R=0.0)
    add("asymmetric_sizes", p0=.25, rho=.5, M_root=20, M_intermediate=12,
        M_sister=8, M_terminal=(8, 10, 12, 16))
    add("asymmetric_sizes_selection", p0=.75, rho=.2, S=.8, I=.1,
        M_root=12, M_intermediate=20, M_sister=16,
        M_terminal=(10, 14, 18, 12), selection=-.25)
    return rows


def _pi(condition: dict[str, object]) -> np.ndarray:
    pi = np.zeros(int(condition["M_root"]) + 1)
    pi[int(condition["k0"])] = 1.0
    return pi


def unbalanced_result(condition: dict[str, object]):
    S, I = float(condition["sister_length"]), float(condition["intermediate_length"])
    return unbalanced_species_tree_conditional_quartet_probabilities(
        M_root=int(condition["M_root"]), M_intermediate=int(condition["M_intermediate"]),
        M_sister=int(condition["M_sister"]), M_terminal=tuple(condition["M_terminal"]),
        rho=float(condition["rho"]), root_length=float(condition["root_length"]),
        intermediate_length=I, sister_length=S, terminal_lengths=(0.0, 0.0, S, S + I),
        pi_old=_pi(condition), selection=float(condition["selection"]),
    )


def balanced_result(condition: dict[str, object]):
    """Matched time-depth comparison: left split at S, right split at S+I."""
    S, I = float(condition["sister_length"]), float(condition["intermediate_length"])
    return full_species_tree_quartet_probabilities(
        M_root=int(condition["M_root"]), M_left=int(condition["M_sister"]),
        M_right=int(condition["M_intermediate"]), M_terminal=tuple(condition["M_terminal"]),
        rho=float(condition["rho"]), root_length=float(condition["root_length"]),
        sister_lengths=(S, S + I), terminal_lengths=(0.0, 0.0, 0.0, 0.0),
        pi_old=_pi(condition), selection=float(condition["selection"]),
    )


def analytic_rows() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    effects: list[dict[str, object]] = []
    for condition in condition_grid():
        unbalanced = unbalanced_result(condition)
        balanced = balanced_result(condition)
        for index, pattern in enumerate(PATTERNS):
            q = unbalanced.conditional_quartet_probabilities[index]
            qb = balanced.conditional_quartet_probabilities[index]
            pz = float(unbalanced.pattern_probabilities[index])
            rows.append({
                "condition": condition["condition"], "pattern": pattern,
                "pattern_probability": pz, "q1": q[0], "q2": q[1], "q3": q[2],
                "delta_q2_minus_q3": q[1] - q[2],
                "gamma_anomaly_margin": max(q[1], q[2]) - q[0],
                "balanced_q1": qb[0], "balanced_q2": qb[1], "balanced_q3": qb[2],
                "balanced_delta": qb[1] - qb[2],
                "balanced_gamma": max(qb[1], qb[2]) - qb[0],
                "delta_tree_difference": (q[1] - q[2]) - (qb[1] - qb[2]),
                "gamma_tree_difference": (max(q[1], q[2]) - q[0]) - (max(qb[1], qb[2]) - qb[0]),
            })
        effects.append({
            "condition": condition["condition"],
            "max_abs_delta": float(np.max(np.abs(unbalanced.conditional_quartet_probabilities[:, 1] - unbalanced.conditional_quartet_probabilities[:, 2]))),
            "max_gamma": float(np.max(np.max(unbalanced.conditional_quartet_probabilities[:, 1:], axis=1) - unbalanced.conditional_quartet_probabilities[:, 0])),
            "pattern_at_max_abs_delta": PATTERNS[int(np.argmax(np.abs(unbalanced.conditional_quartet_probabilities[:, 1] - unbalanced.conditional_quartet_probabilities[:, 2])))],
            "pattern_at_max_gamma": PATTERNS[int(np.argmax(np.max(unbalanced.conditional_quartet_probabilities[:, 1:], axis=1) - unbalanced.conditional_quartet_probabilities[:, 0]))],
            "pattern_probability_at_max_delta": float(unbalanced.pattern_probabilities[int(np.argmax(np.abs(unbalanced.conditional_quartet_probabilities[:, 1] - unbalanced.conditional_quartet_probabilities[:, 2])))]),
            "pattern_probability_at_max_gamma": float(unbalanced.pattern_probabilities[int(np.argmax(np.max(unbalanced.conditional_quartet_probabilities[:, 1:], axis=1) - unbalanced.conditional_quartet_probabilities[:, 0]))]),
        })
    return rows, effects


def parameter_rows() -> list[dict[str, object]]:
    """Small one-pattern analytical sweeps used for the parameter figure."""
    rows: list[dict[str, object]] = []
    for rho in (0.0, .05, .2, .5, 1.0):
        condition = {"condition": f"rho_{rho:g}", "p0": .5, "rho": rho,
                     "sister_length": .3, "intermediate_length": .2,
                     "root_length": .5, "M_root": 10, "M_intermediate": 10,
                     "M_sister": 10, "M_terminal": (10, 10, 10, 10),
                     "selection": 0.0, "k0": 5}
        result = unbalanced_result(condition); i = PATTERNS.index("0101")
        for qi, topology in enumerate(("q1", "q2", "q3")):
            rows.append({"parameter": "rho", "value": rho, "pattern": "0101",
                         "topology": topology, "probability": result.conditional_quartet_probabilities[i, qi]})
    for intermediate in (0.0, .1, .2, .5, .8):
        condition = {"condition": f"I_{intermediate:g}", "p0": .5, "rho": .17,
                     "sister_length": .3, "intermediate_length": intermediate,
                     "root_length": .5, "M_root": 10, "M_intermediate": 10,
                     "M_sister": 10, "M_terminal": (10, 10, 10, 10),
                     "selection": 0.0, "k0": 5}
        result = unbalanced_result(condition); i = PATTERNS.index("0101")
        for qi, topology in enumerate(("q1", "q2", "q3")):
            rows.append({"parameter": "intermediate_length", "value": intermediate, "pattern": "0101",
                         "topology": topology, "probability": result.conditional_quartet_probabilities[i, qi]})
    return rows


def _tree(msrc_sim: Path, condition: dict[str, object]):
    sys.path.insert(0, str(msrc_sim))
    from msrcsim.species_tree import SpeciesTree

    S, I = float(condition["sister_length"]), float(condition["intermediate_length"])
    M_root = int(condition["M_root"])
    newick = f"(((1:0,2:0)N12:{S},3:{S})N123:{I},4:{S + I})ROOT;"
    sizes = {"ROOT": M_root, "N123": int(condition["M_intermediate"]),
             "N12": int(condition["M_sister"]), "1": condition["M_terminal"][0],
             "2": condition["M_terminal"][1], "3": condition["M_terminal"][2],
             "4": condition["M_terminal"][3]}
    selection = float(condition["selection"])
    params = {name: {"effective_population_size": int(M) // 2,
                     "selection_coefficient": selection} for name, M in sizes.items()}
    return SpeciesTree(newick, M_root // 2, float(condition["root_length"]), branch_parameters=params)


def simulate_joint(msrc_sim: Path, condition: dict[str, object], replicates: int, seed: int) -> np.ndarray:
    sys.path.insert(0, str(msrc_sim))
    from msrcsim.moran import simulate_frequency_history
    from msrcsim.rearrangement import Rearrangement
    from msrcsim.structured_coalescent import simulate_genealogy

    tree = _tree(msrc_sim, condition)
    rng = np.random.default_rng(seed)
    rearrangement = Rearrangement("conditional_validation", "inversion", "ROOT", 0,
                                  int(condition["k0"]), float(condition["selection"]))
    counts = np.zeros((16, 3), dtype=int)
    for locus in range(replicates):
        history = simulate_frequency_history(tree, rearrangement, rng)
        sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
        pattern = "".join(str(sampled[taxon]) for taxon in tree.taxa)
        result = simulate_genealogy(locus, tree, history, sampled, float(condition["rho"]), 1.0, rng, record_events=False)
        counts[PATTERNS.index(pattern), int(result.topology_index)] += 1
    return counts


def validation_rows(msrc_sim: Path, replicates: int, seed: int) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for ci, condition in enumerate(condition_grid()[:3]):
        analytic = unbalanced_result(condition)
        observed_counts = simulate_joint(msrc_sim, condition, replicates, seed + ci * 1009)
        expected = analytic.joint * replicates
        chi = float(np.sum((observed_counts - expected) ** 2 / expected, where=expected > 0))
        for zi, pattern in enumerate(PATTERNS):
            for qi, topology in enumerate(("q1", "q2", "q3")):
                p = float(analytic.joint[zi, qi])
                obs = float(observed_counts[zi, qi] / replicates)
                se = float(np.sqrt(max(p * (1.0 - p), 0.0) / replicates))
                output.append({
                    "condition": condition["condition"], "pattern": pattern, "topology": topology,
                    "analytic_joint": p, "observed_joint": obs,
                    "count": int(observed_counts[zi, qi]), "mc_se": se,
                    "absolute_error": abs(obs - p),
                    "standardized_error": (obs - p) / se if se else np.nan,
                    "pattern_count": int(observed_counts[zi].sum()), "replicates": replicates,
                    "joint_pearson_chi2": chi, "joint_pearson_pvalue": float(chi2.sf(chi, 47)),
                })
    return output


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--msrc-sim", type=Path, default=Path("/Users/ytabatabaee/Desktop/msrc-sim"))
    parser.add_argument("--replicates", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=20261025)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "results")
    args = parser.parse_args()
    start = time.perf_counter()
    rows, effects = analytic_rows()
    parameters = parameter_rows()
    validation = validation_rows(args.msrc_sim, args.replicates, args.seed)
    write_tsv(args.output_dir / "unbalanced_conditional_analytic.tsv", rows)
    write_tsv(args.output_dir / "unbalanced_conditional_effects.tsv", effects)
    write_tsv(args.output_dir / "unbalanced_conditional_parameters.tsv", parameters)
    write_tsv(args.output_dir / "unbalanced_conditional_simulation_validation.tsv", validation)
    max_error = max(float(row["absolute_error"]) for row in validation)
    standardized = [abs(float(row["standardized_error"])) for row in validation if np.isfinite(float(row["standardized_error"]))]
    (args.output_dir / "unbalanced_conditional_validation_summary.tsv").write_text(
        "metric\tvalue\n" + "\n".join([
            f"analytic_conditions\t{len(condition_grid())}",
            f"simulation_conditions\t3",
            f"replicates_per_simulation_condition\t{args.replicates}",
            f"max_joint_absolute_error\t{max_error}",
            f"joint_cells_within_3_se\t{float(np.mean(np.asarray(standardized) <= 3.0))}",
            f"runtime_seconds\t{time.perf_counter() - start}",
        ]) + "\n"
    )
    print(f"wrote {len(rows)} analytic and {len(validation)} simulation-validation rows")


if __name__ == "__main__":
    main()
