#!/usr/bin/env python3
"""Run the missing full-data t-haplotype-only four-group ASTRAL analysis."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import (  # noqa: E402
    PROCESSED,
    RESULTS,
    STAGE2_DATA,
    all_required_groups_present,
    find_astral,
    four_taxon_topology_from_newick,
    load_stage1,
    mapping_by_tip,
    parse_astral_result,
    prune_newick,
    rel,
    run_astral,
    tips_by_species,
    write_map,
    write_tree_file,
    write_tsv,
)


TREE_INPUT = STAGE2_DATA / "T_T_ONLY.tre"
MAP_INPUT = STAGE2_DATA / "T_T_ONLY_subspecies.map"
TREE_OUTPUT = RESULTS / "stage1_aster" / "T_T_ONLY_subspecies.nwk"
COVERAGE_OUTPUT = RESULTS / "stage2_t_only_sample_coverage.tsv"
PROVENANCE_OUTPUT = RESULTS / "stage2_t_only_astral_provenance.json"
SUMMARY_OUTPUT = RESULTS / "stage2_astral_summary.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=2)
    return parser.parse_args()


def biological_q(topology: str, q1: float, q2: float, q3: float) -> dict[str, float]:
    if topology == "Q_SPECIES":
        return {"q_species": q1, "q_t_alt": q2, "q_other": q3}
    if topology == "Q_T_ALT":
        return {"q_species": q2, "q_t_alt": q1, "q_other": q3}
    if topology == "Q_OTHER":
        return {"q_species": q2, "q_t_alt": q3, "q_other": q1}
    raise ValueError(f"Cannot map quartet annotations for topology {topology}")


def coverage_rows(mapping: list[dict[str, str]], tips: set[str]) -> list[dict[str, object]]:
    rows = []
    for species in sorted({row["subspecies"] for row in mapping if row["tree_tip"] in tips}):
        subset = [row for row in mapping if row["tree_tip"] in tips and row["subspecies"] == species]
        rows.append({
            "treatment": "T_T_ONLY",
            "group": species,
            "n_tips": len(subset),
            "t_status_counts": ";".join(f"{key}:{value}" for key, value in sorted(Counter(row["t_status"] for row in subset).items())),
            "tips": ",".join(sorted(row["tree_tip"] for row in subset)),
        })
    return rows


def run(args: argparse.Namespace) -> int:
    meta, trees, mapping = load_stage1()
    groups = tips_by_species(mapping, "T_ONLY")
    if not all_required_groups_present(groups):
        raise RuntimeError(f"T-only four-group coverage is incomplete: {groups}")
    tips = {tip for group in groups.values() for tip in group}
    if not tips:
        raise RuntimeError("No t-only tips were recovered")

    # The tree set is frozen; only the taxon set changes for this treatment.
    write_tree_file(TREE_INPUT, [prune_newick(newick, tips) for newick in trees])
    write_map(MAP_INPUT, mapping, tips, "subspecies")
    astral = find_astral()
    result = run_astral(astral, TREE_INPUT, MAP_INPUT, TREE_OUTPUT, threads=args.threads)
    topology = str(result["topology"])
    q = biological_q(topology, float(result["q1"]), float(result["q2"]), float(result["q3"]))

    write_tsv(COVERAGE_OUTPUT, coverage_rows(mapping, tips), ["treatment", "group", "n_tips", "t_status_counts", "tips"])
    provenance = {
        "treatment": "T_T_ONLY",
        "definition": "pseudo-t_haplotype tips plus the source-derived Mus spretus outgroup",
        "n_gene_trees_input": len(trees),
        "n_gene_trees_retained": result.get("n_gene_trees", len(trees)),
        "n_input_tips": len(tips),
        "groups": {key: sorted(value) for key, value in groups.items()},
        "tree_input": rel(TREE_INPUT),
        "map_input": rel(MAP_INPUT),
        "tree_output": rel(TREE_OUTPUT),
        "coverage_output": rel(COVERAGE_OUTPUT),
        "astral_executable": str(astral),
        "astral_result": result,
        "topology": topology,
        "biological_quartets": q,
        "criteria": "All 4,046 original local trees were pruned to the validated T_ONLY tip set; no tree was removed.",
    }
    PROVENANCE_OUTPUT.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")

    existing = [row for row in __import__("csv").DictReader((SUMMARY_OUTPUT).open(), delimiter="\t")]
    existing = [row for row in existing if row["treatment"] != "T_T_ONLY"]
    existing.append({
        "treatment": "T_T_ONLY",
        "n_gene_trees": str(result.get("n_gene_trees", len(trees))),
        "topology": topology,
        "rf_distance_vs_P0": "NA",
        "changed_splits_vs_P0": "NA",
        "CU_length": str(result["CU_length"]),
        "localPP": str(result["localPP"]),
        "q1": str(result["q1"]),
        "q2": str(result["q2"]),
        "q3": str(result["q3"]),
        "q_species": str(q["q_species"]),
        "q_t_alt": str(q["q_t_alt"]),
        "q_other": str(q["q_other"]),
        "status": "full_data_t_only_focal_split_annotation",
    })
    fields = list(existing[0])
    write_tsv(SUMMARY_OUTPUT, existing, fields)
    print(f"Wrote {TREE_OUTPUT}")
    print(f"T-only topology: {topology}; q_species={q['q_species']:.6f}; q_t_alt={q['q_t_alt']:.6f}; q_other={q['q_other']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run(parse_args()))
