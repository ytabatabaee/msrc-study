from __future__ import annotations

import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "empirical/neoaves_chr4/scripts/09_full_taxon_background.py"
OUT = ROOT / "empirical/neoaves_chr4/results/stage9_full_taxon_background"


def load_module():
    spec = importlib.util.spec_from_file_location("stage9_full_taxon_background", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_full_taxon_universe_and_outputs_exist():
    cov = rows(OUT / "per_window_taxon_coverage.tsv")
    presence = rows(OUT / "species_presence_summary.tsv")
    assert len(cov) == 6140
    assert len(presence) == 363
    assert max(int(r["n_expected_taxa"]) for r in cov) == 363


def test_stage6_fixed_distance_regression():
    module = load_module()
    new = {r["Gene"]: r for r in rows(OUT / "background_distance_estimators.tsv")}
    old = {r["Gene"]: r for r in rows(module.STAGE6)}
    differences = []
    for gene, row in new.items():
        if row["B47"] not in ("nan", "NaN", "") and old[gene]["B_fixed_mapped"] not in ("nan", "NaN", ""):
            differences.append(abs(float(row["B47"]) - float(old[gene]["B_fixed_mapped"])))
    assert len(differences) == 2707
    assert max(differences) < 1e-12


def test_branch_contribution_mean_matches_pairwise_mean():
    module = load_module()
    tree = module.Tree("(a:1,(b:2,c:3):4);")
    # Pair distances are a-b=7, a-c=8, b-c=5.
    assert abs(tree.subset_pair_mean(["a", "b", "c"]) - (7 + 8 + 5) / 3) < 1e-12


def test_fixed_broad_selection_is_presence_based():
    selection = rows(OUT / "fixed_taxon_selection.tsv")
    broad = next(r for r in selection if r["selection"] == "fixed_broad_95pct")
    assert broad["threshold"] == "presence in >=95% of chr4 local trees"
    assert int(broad["n_taxa"]) == 108
