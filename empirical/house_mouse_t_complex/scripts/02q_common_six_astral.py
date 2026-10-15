#!/usr/bin/env python3
"""Run a common-six-population ASTRAL comparison for all three treatments."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter

SCRIPT_DIR = __import__("pathlib").Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import (  # noqa: E402
    RESULTS,
    STAGE2_DATA,
    find_astral,
    load_stage1,
    prune_newick,
    rel,
    rf_distance,
    run_astral,
    tree_bipartitions,
    write_map,
    write_tree_file,
    write_tsv,
)


OUT_DIR = RESULTS / "stage2_common6_population_astral"
SUMMARY = RESULTS / "stage2_common6_population_astral.tsv"
COVERAGE = RESULTS / "stage2_common6_population_coverage.tsv"
PROVENANCE = RESULTS / "stage2_common6_population_provenance.json"
EXCLUDED_POPULATION = "SPRE"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=2)
    return parser.parse_args()


def run(args: argparse.Namespace) -> int:
    meta, trees, mapping = load_stage1()
    usable = [row for row in mapping if row.get("mapping_confidence") in {"strong", "exact"} and row.get("population") != EXCLUDED_POPULATION]
    populations = sorted({row["population"] for row in usable})
    if len(populations) != 6:
        raise RuntimeError(f"Expected six common populations, observed {populations}")
    by_treatment = {
        "C6_STANDARD": {row["tree_tip"] for row in usable if row["t_status"] == "standard_noncarrier"},
        "C6_T_ONLY": {row["tree_tip"] for row in usable if row["t_status"] == "pseudo-t_haplotype"},
        "C6_MIXED": {row["tree_tip"] for row in usable},
    }
    for name, tips in by_treatment.items():
        observed = {row["population"] for row in usable if row["tree_tip"] in tips}
        if observed != set(populations):
            raise RuntimeError(f"{name} does not cover all common populations: {observed}")

    astral = find_astral()
    rows = []
    newicks = {}
    for name, tips in by_treatment.items():
        tree_file = STAGE2_DATA / f"{name}.tre"
        map_file = STAGE2_DATA / f"{name}.map"
        write_tree_file(tree_file, [prune_newick(newick, tips) for newick in trees])
        write_map(map_file, usable, tips, "population")
        result = run_astral(astral, tree_file, map_file, OUT_DIR / f"{name}.nwk", threads=args.threads)
        newicks[name] = str(result["newick"])
        rows.append({
            "treatment": name,
            "n_gene_trees": result.get("n_gene_trees", len(trees)),
            "n_populations": result.get("n_species", len(populations)),
            "n_tips": len(tips),
            "CU_length": result.get("CU_length", "NA"),
            "localPP": result.get("localPP", "NA"),
            "q1": result.get("q1", "NA"),
            "q2": result.get("q2", "NA"),
            "q3": result.get("q3", "NA"),
            "newick": result.get("newick", ""),
        })

    baseline = newicks["C6_STANDARD"]
    for row in rows:
        row["rf_distance_vs_C6_STANDARD"] = 0 if row["treatment"] == "C6_STANDARD" else rf_distance(baseline, newicks[row["treatment"]])
        row["changed_splits_vs_C6_STANDARD"] = 0 if row["treatment"] == "C6_STANDARD" else len(tree_bipartitions(baseline) ^ tree_bipartitions(newicks[row["treatment"]]))
    write_tsv(SUMMARY, rows, ["treatment", "n_gene_trees", "n_populations", "n_tips", "rf_distance_vs_C6_STANDARD", "changed_splits_vs_C6_STANDARD", "CU_length", "localPP", "q1", "q2", "q3", "newick"])

    coverage_rows = []
    for name, tips in by_treatment.items():
        for population in populations:
            subset = [row for row in usable if row["tree_tip"] in tips and row["population"] == population]
            coverage_rows.append({
                "treatment": name,
                "population": population,
                "n_tips": len(subset),
                "t_status_counts": ";".join(f"{key}:{value}" for key, value in sorted(Counter(row["t_status"] for row in subset).items())),
                "tips": ",".join(sorted(row["tree_tip"] for row in subset)),
            })
    write_tsv(COVERAGE, coverage_rows, ["treatment", "population", "n_tips", "t_status_counts", "tips"])
    PROVENANCE.write_text(json.dumps({
        "analysis": "common-six-population three-treatment ASTRAL",
        "excluded_population": EXCLUDED_POPULATION,
        "populations": populations,
        "n_gene_trees_input": len(trees),
        "n_gene_trees_retained": len(trees),
        "treatments": {name: {"n_tips": len(tips), "tips": sorted(tips)} for name, tips in by_treatment.items()},
        "outputs": {"summary": rel(SUMMARY), "coverage": rel(COVERAGE), "provenance": rel(PROVENANCE)},
        "criteria": "All three analyses use the same six populations and all 4,046 original local trees; SPRE is excluded from every treatment.",
        "astral_executable": str(astral),
    }, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {SUMMARY}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run(parse_args()))
