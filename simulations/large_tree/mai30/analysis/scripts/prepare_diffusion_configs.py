#!/usr/bin/env python3
"""Prepare, but do not execute, original-scale wf_diffusion configs."""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
MSRC_ROOT = Path("/Users/ytabatabaee/Desktop/msrc-sim")
SELECTED = [("outgroup.0.species.1.5.genes.1.5", "007"), ("outgroup.0.species.1.5.genes.1.5", "042"), ("outgroup.0.species.1.5.genes.1.5", "074"), ("outgroup.0.species.5.genes.5", "005")]


def main() -> None:
    params = {(r["condition"], r["replicate"]): r for r in csv.DictReader((ROOT / "datasets/processed/mai30_replicate_parameters.tsv").open(), delimiter="\t")}
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=MSRC_ROOT, text=True).strip()
    out_root = ROOT / "datasets/configs/mai30_diffusion_experimental"; out_root.mkdir(parents=True, exist_ok=True)
    rows = []
    for condition, rep in SELECTED:
        row = params[(condition, rep)]
        nh = int(float(row["haploid_effective_population_size"]))
        ne = int(round(nh / 2))
        config = {"mode": "arbitrary_tree", "seed": 20261011 + int(rep), "num_loci": 20, "species_tree": {"newick": Path(row["species_tree"]).read_text().strip(), "root_extension": 40 * ne, "default_effective_population_size": ne, "ultrametric_tolerance": 1e-6}, "rearrangement": {"id": "mai_diffusion_inv", "type": "inversion", "origin_branch": "ROOT", "origin_time_from_branch_start": 1, "initial_copy_count": ne, "selection": {"model": "genic", "coefficient": 0.0}}, "recombination": {"baseline_rate": 0.02, "effective_cross_arrangement_fraction": 0.5}, "genome": {"length_bp": 20, "rearrangement_interval": {"start": 0, "end": 20, "id": "mai_diffusion_inv"}, "breakpoint_rate_per_bp": 0.0}, "population_process": {"model": "wf_diffusion", "diffusion": {"tolerance": 0.001, "max_step_cu": 0.01}}, "ordinary_msc_outside": True, "output": {"directory": str(ROOT / "datasets/processed/mai30_diffusion_experimental" / f"{condition}__replicate_{rep}"), "plots": False}, "provenance": {"source": "Mai et al. 2017 D1 fixed backbone", "simphy_haploid_ne": nh, "msrc_diploid_ne": ne, "simulator_required_commit": commit, "status": "experimental; do not use in primary figure until genealogy-level diffusion validation passes"}}
        path = out_root / f"{condition}__replicate_{rep}.yaml"; path.write_text(yaml.safe_dump(config, sort_keys=False)); rows.append({"config": str(path), "required_msrc_sim_commit": commit, "status": "prepared_not_executed"})
    with (out_root / "MANIFEST.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    (out_root / "README.md").write_text("Experimental original-scale wf_diffusion configurations. They are prepared only; do not execute or include their results in the primary scaled-WF figure until validation/wf_diffusion/DIFFUSION_VALIDATION_REPORT.md changes its experimental status.\n")
    print(f"wrote {len(rows)} diffusion configs; required commit {commit}")


if __name__ == "__main__": main()
