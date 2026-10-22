#!/usr/bin/env python3
"""Write the reproducible origin-placement design without simulating genealogies."""

from __future__ import annotations

import json
from pathlib import Path

from run_origin_discovery import OUT, RHO_GRID, TREES, executable_rhos, load_base, scenario_rows, write_tsv, exposure
from msrcsim import SpeciesTree


def main():
    scenarios = []
    exposures = []
    for condition, replicate in TREES:
        base = load_base(condition, replicate)
        tree = SpeciesTree(base["species_tree"]["newick"], base["species_tree"]["default_effective_population_size"], base["species_tree"]["root_extension"], ultrametric_tolerance=1e-6)
        rows = scenario_rows(condition, replicate, base, tree)
        scenarios.extend(rows)
        exposures.extend(exposure(tree, row) for row in rows)
    write_tsv(OUT / "origin_scenarios.tsv", scenarios)
    write_tsv(OUT / "quartet_exposure.tsv", exposures)
    grid = []
    for row in scenarios:
        for rho in RHO_GRID:
            grid.append({"condition": row["condition"], "replicate": row["replicate"], "scenario": row["scenario"], "origin_branch": row["origin_branch"], "origin_age": row["origin_age"], "descendant_taxa": row["descendant_taxa"], "rho": rho, "n_histories": 1, "loci_per_history": 1, "status": "executed_forward_only" if rho in executable_rhos(row["scenario"]) else "planned_not_executed"})
    write_tsv(OUT / "parameter_grid.tsv", grid)
    (OUT / "design_manifest.json").write_text(json.dumps({"rho_grid": RHO_GRID, "executed_genealogy_rhos": "rho<=0.02 in design; local genealogy run did not complete", "scenario_rules": "see origin_scenarios.tsv", "total_quartets": 27405, "structural_exposure_calculated_for_all_quartets": True}, indent=2) + "\n")


if __name__ == "__main__":
    main()
