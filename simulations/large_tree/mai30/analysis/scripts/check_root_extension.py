#!/usr/bin/env python3
"""Check whether ordinary MSC genealogy generation uses root_extension."""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

import numpy as np

MSRC_ROOT = Path("/Users/ytabatabaee/Desktop/msrc-sim")
STUDY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(MSRC_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from msrcsim import SpeciesTree, simulate_msc_genealogy  # noqa: E402
from population_mapping import map_simphy_haploid_to_msrc  # noqa: E402


def main() -> None:
    condition = "outgroup.0.species.1.5.genes.1.5"
    replicate = "042"
    params = {}
    with (STUDY_ROOT / "datasets/processed/mai30_replicate_parameters.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            params[(row["condition"], row["replicate"])] = row
    row = params[(condition, replicate)]
    mapping = map_simphy_haploid_to_msrc(int(float(row["haploid_effective_population_size"])))
    newick = (STUDY_ROOT / row["species_tree"]).read_text().strip()
    seed = 424242
    records = []
    for root_cu in (5, 10, 20, 40):
        tree = SpeciesTree(newick, mapping.msrc_diploid_ne_integer, root_cu * 2 * mapping.msrc_diploid_ne_integer, ultrametric_tolerance=1e-6)
        rng = np.random.default_rng(seed)
        trees = [simulate_msc_genealogy(i, tree, rng, record_events=False).newick for i in range(100)]
        digest = hashlib.sha256("\n".join(trees).encode()).hexdigest()
        records.append({
            "condition": condition,
            "replicate": replicate,
            "root_extension_coalescent_units": root_cu,
            "root_extension_generations": root_cu * 2 * mapping.msrc_diploid_ne_integer,
            "seed": seed,
            "n_trees": len(trees),
            "newick_sha256": digest,
            "multiple_lineages_at_extension_endpoint": "not applicable: ordinary MSC has no finite endpoint",
        })
    out = STUDY_ROOT / "analysis/results/root_extension_same_seed.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    print(out)
    print("identical_hashes", len({r["newick_sha256"] for r in records}) == 1)


if __name__ == "__main__":
    main()
