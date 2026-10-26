"""Targeted extended validation of the exact pectinate quartet recursion.

This driver is deliberately manifest based: the compact condition set is frozen
in the script, every locus gets a fresh Moran history, and the simulator is
kept independent of the operator calculation.  It does not submit jobs.
"""

from __future__ import annotations

import argparse
import csv
import resource
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import chi2

try:
    from .unbalanced_species_tree import unbalanced_species_tree_quartet_probabilities
except ImportError:  # pragma: no cover
    from unbalanced_species_tree import unbalanced_species_tree_quartet_probabilities


TOPOLOGIES = ("q1", "q2", "q3")


def _tree(msrc_sim: Path, condition: dict[str, object]):
    sys.path.insert(0, str(msrc_sim))
    from msrcsim.species_tree import SpeciesTree

    M_root = int(condition["M_root"])
    M_intermediate = int(condition["M_intermediate"])
    M_sister = int(condition["M_sister"])
    M_terminal = tuple(int(x) for x in condition["M_terminal"])
    S = float(condition["sister_length"])
    I = float(condition["intermediate_length"])
    R = float(condition["root_length"])
    selection = float(condition["selection"])
    if any(m <= 0 or m % 2 for m in (M_root, M_intermediate, M_sister, *M_terminal)):
        raise ValueError("all chromosome-copy counts must be positive and even")
    # The Newick tree is ultrametric: tip 3 spans the sister branch and tip 4
    # spans the sister plus intermediate branches before entering the root.
    newick = f"(((1:0,2:0)N12:{S},3:{S})N123:{I},4:{S + I})ROOT;"
    branch_ne = {
        "ROOT": M_root // 2, "N123": M_intermediate // 2,
        "N12": M_sister // 2,
        "1": M_terminal[0] // 2, "2": M_terminal[1] // 2,
        "3": M_terminal[2] // 2, "4": M_terminal[3] // 2,
    }
    parameters = {
        name: {"effective_population_size": ne, "selection_coefficient": selection}
        for name, ne in branch_ne.items()
    }
    return SpeciesTree(newick, M_root // 2, R, branch_parameters=parameters)


def simulate_condition(msrc_sim: Path, condition: dict[str, object], replicates: int, seed: int) -> np.ndarray:
    sys.path.insert(0, str(msrc_sim))
    from msrcsim.moran import simulate_frequency_history
    from msrcsim.rearrangement import Rearrangement
    from msrcsim.structured_coalescent import simulate_genealogy

    tree = _tree(msrc_sim, condition)
    M_root = int(condition["M_root"])
    k0 = int(condition["k0"])
    rho = float(condition["rho"])
    rng = np.random.default_rng(seed)
    rearrangement = Rearrangement("pectinate_extended_validation", "inversion", "ROOT", 0, k0, float(condition["selection"]))
    counts = np.zeros(3, dtype=int)
    for locus in range(replicates):
        history = simulate_frequency_history(tree, rearrangement, rng)
        sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
        result = simulate_genealogy(locus, tree, history, sampled, rho, 1.0, rng, record_events=False)
        counts[int(result.topology_index)] += 1
    return counts


def conditions() -> list[dict[str, object]]:
    """Return the compact, non-Cartesian validation grid."""
    out: list[dict[str, object]] = []

    def add(label: str, *, M_root=10, M_intermediate=10, M_sister=10,
            M_terminal=(10, 10, 10, 10), p0=.5, rho=.17,
            sister=.3, intermediate=.2, root=.5, selection=0.0):
        out.append({
            "condition": label, "M_root": M_root, "M_intermediate": M_intermediate,
            "M_sister": M_sister, "M_terminal": tuple(M_terminal),
            "p0": p0, "rho": rho, "sister_length": sister,
            "intermediate_length": intermediate, "root_length": root,
            "selection": selection,
        })

    # Starting-frequency boundaries and interior values.
    for p in (0.0, .1, .25, .5, .75, .9, 1.0):
        add(f"p0_{p:g}", p0=p)
    # Independent branch-duration asymmetries, including zero-duration limits.
    for label, S, I, R in (
        ("dur_zero", 0.0, 0.0, 0.0),
        ("dur_I_long", .1, .8, .3),
        ("dur_S_long", .8, .1, .3),
        ("dur_all_long", .5, .5, 1.0),
    ):
        add(label, sister=S, intermediate=I, root=R)
    # Switching and selection, keeping the branch geometry fixed.
    for rho in (0.0, .05, .2, .5, 1.0):
        add(f"rho_{rho:g}", rho=rho)
    add("selection_positive", selection=.25)
    add("selection_negative", selection=-.25)
    # Genuinely asymmetric chromosome-copy counts.
    add("M_asym_1", M_root=20, M_intermediate=12, M_sister=8, M_terminal=(8, 10, 12, 16), p0=.25, rho=.2)
    add("M_asym_2", M_root=12, M_intermediate=20, M_sister=16, M_terminal=(10, 14, 18, 12), p0=.75, rho=.5)
    add("M_asym_3", M_root=16, M_intermediate=8, M_sister=20, M_terminal=(12, 8, 10, 18), p0=.5, rho=.05, sister=.8, intermediate=.1)
    add("M_asym_4", M_root=8, M_intermediate=16, M_sister=12, M_terminal=(8, 16, 12, 20), p0=.1, rho=1.0, selection=.25)
    add("M_asym_5", M_root=20, M_intermediate=16, M_sister=8, M_terminal=(20, 10, 8, 16), p0=.9, rho=.2, selection=-.25)
    for row in out:
        row["k0"] = int(round(int(row["M_root"]) * float(row["p0"])))
    return out


def _theory(condition: dict[str, object]):
    M_root = int(condition["M_root"])
    pi = np.zeros(M_root + 1)
    pi[int(condition["k0"])] = 1.0
    # These are the terminal lengths in the simulator's ultrametric tree.
    terminal = (0.0, 0.0, float(condition["sister_length"]),
                float(condition["sister_length"]) + float(condition["intermediate_length"]))
    return unbalanced_species_tree_quartet_probabilities(
        M_root=M_root, M_intermediate=int(condition["M_intermediate"]),
        M_sister=int(condition["M_sister"]), M_terminal=tuple(condition["M_terminal"]),
        rho=float(condition["rho"]), root_length=float(condition["root_length"]),
        intermediate_length=float(condition["intermediate_length"]),
        sister_length=float(condition["sister_length"]), terminal_lengths=terminal,
        pi_old=pi, selection=float(condition["selection"]),
    )


def run(msrc_sim: Path, replicates: int, seed: int) -> tuple[list[dict[str, object]], dict[str, object], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    condition_rows: list[dict[str, object]] = []
    start = time.perf_counter()
    for index, condition in enumerate(conditions()):
        analytic_result = _theory(condition)
        analytic = np.asarray(analytic_result.quartet_probabilities, dtype=float)
        counts = simulate_condition(msrc_sim, condition, replicates, seed + 1009 * index)
        observed = counts / replicates
        errors = observed - analytic
        se = np.sqrt(np.maximum(analytic * (1.0 - analytic), 0.0) / replicates)
        standardized = np.divide(errors, se, out=np.full(3, np.nan), where=se > 0)
        expected = replicates * analytic
        pearson = float(np.sum((counts - expected) ** 2 / expected, where=expected > 0))
        p_value = float(chi2.sf(pearson, 2))
        max_error = float(np.max(np.abs(errors)))
        asymmetry = max(condition["M_root"], condition["M_intermediate"], condition["M_sister"], *condition["M_terminal"]) / min(condition["M_root"], condition["M_intermediate"], condition["M_sister"], *condition["M_terminal"])
        condition_rows.append({
            "condition_index": index + 1, "condition": condition["condition"],
            "M_root": condition["M_root"], "M_intermediate": condition["M_intermediate"],
            "M_sister": condition["M_sister"], "M_terminal": ",".join(map(str, condition["M_terminal"])),
            "p0": condition["p0"], "k0": condition["k0"], "rho": condition["rho"],
            "selection": condition["selection"], "sister_length": condition["sister_length"],
            "intermediate_length": condition["intermediate_length"], "root_length": condition["root_length"],
            "population_size_ratio": asymmetry, "replicates": replicates,
            "pearson_chi2": pearson, "pearson_pvalue": p_value,
            "max_absolute_error": max_error, "probability_sum_error": float(analytic.sum() - 1.0),
        })
        for r, topology in enumerate(TOPOLOGIES):
            rows.append({
                "condition_index": index + 1, "condition": condition["condition"],
                "M_root": condition["M_root"], "M_intermediate": condition["M_intermediate"],
                "M_sister": condition["M_sister"], "M_terminal": ",".join(map(str, condition["M_terminal"])),
                "p0": condition["p0"], "k0": condition["k0"], "rho": condition["rho"],
                "selection": condition["selection"], "sister_length": condition["sister_length"],
                "intermediate_length": condition["intermediate_length"], "root_length": condition["root_length"],
                "population_size_ratio": asymmetry, "topology": topology,
                "analytic": analytic[r], "observed": observed[r], "count": int(counts[r]),
                "absolute_error": abs(float(errors[r])), "mc_se": float(se[r]),
                "standardized_error": float(standardized[r]), "replicates": replicates,
            })
    elapsed = time.perf_counter() - start
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes; Linux reports KiB.
    rss_mib = rss / (1024 * 1024) if sys.platform == "darwin" else rss / 1024
    all_errors = np.asarray([float(r["absolute_error"]) for r in rows])
    all_std = np.asarray([abs(float(r["standardized_error"])) for r in rows if np.isfinite(float(r["standardized_error"]))])
    summary = {
        "conditions": len(condition_rows), "replicates_per_condition": replicates,
        "total_genealogies": len(condition_rows) * replicates,
        "topology_comparisons": len(rows), "rmse": float(np.sqrt(np.mean(all_errors ** 2))),
        "maximum_absolute_error": float(np.max(all_errors)),
        "fraction_topology_points_within_2_se": float(np.mean(all_std <= 2.0)),
        "fraction_topology_points_within_3_se": float(np.mean(all_std <= 3.0)),
        "fraction_cells_pearson_p_ge_0_05": float(np.mean([float(r["pearson_pvalue"]) >= .05 for r in condition_rows])),
        "runtime_seconds": elapsed, "peak_rss_mib": rss_mib,
    }
    return rows, summary, condition_rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--msrc-sim", type=Path, default=Path("/Users/ytabatabaee/Desktop/msrc-sim"))
    parser.add_argument("--replicates", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261025)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "results")
    args = parser.parse_args()
    rows, summary, condition_rows = run(args.msrc_sim, args.replicates, args.seed)
    write_tsv(args.output_dir / "unbalanced_species_tree_extended_validation.tsv", rows)
    write_tsv(args.output_dir / "unbalanced_species_tree_extended_conditions.tsv", condition_rows)
    (args.output_dir / "unbalanced_species_tree_extended_summary.tsv").write_text(
        "metric\tvalue\n" + "\n".join(f"{key}\t{value}" for key, value in summary.items()) + "\n"
    )
    print(summary)


if __name__ == "__main__":
    main()
