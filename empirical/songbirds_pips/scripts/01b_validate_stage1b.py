#!/usr/bin/env python3
"""Automated validation for Stage 1B preflight/failure-safe outputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "empirical" / "songbirds_pips" / "results"
PROCESSED = ROOT / "data" / "songbirds_pips" / "processed"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def assert_true(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    manifest = json.loads((RESULTS / "stage1b_manifest.json").read_text())
    design = json.loads((RESULTS / "stage1b_design_validation.json").read_text())
    dxy = read_tsv(PROCESSED / "stage1b_shared_pip_dxy.tsv")
    matrix = read_tsv(PROCESSED / "stage1b_shared_pip_dxy_matrix.tsv")
    sfs = read_tsv(RESULTS / "stage1b_pairwise_sfs_status.tsv")
    samples = read_tsv(ROOT / "data" / "songbirds_pips" / "processed" / "stage1a_shared_pip_samples.tsv")

    sample_ids = [row["sample_id"] for row in samples]
    genotypes = [row["pip_genotype"] for row in samples]
    species = [row["species"] for row in samples]

    assert_true(design["exactly_22_samples"], "frozen design must contain exactly 22 samples")
    assert_true(design["aa_11_bb_11"], "frozen design must contain AA=11 and BB=11")
    assert_true(design["eleven_species"], "frozen design must contain 11 species")
    assert_true(design["one_AA_and_one_BB_per_species"], "frozen design must contain one AA and one BB per species")
    assert_true(design["focal_interval_matches_stage1a"], "focal interval must match frozen Stage 1A interval")
    assert_true(design["reference_assembly_GCA_001746935_2"], "reference assembly must be GCA_001746935.2")
    assert_true(len(dxy) == 231, "dxy pair table must enumerate 231 expected pairs")
    assert_true(len(sfs) == 231, "pairwise SFS status table must enumerate 231 expected pairs")
    assert_true(manifest["design_validation"]["no_background_control_analysis_run"], "background/control analysis must not run in Stage 1B")
    assert_true(manifest["design_validation"]["no_candidate_set_changes"], "candidate set must not be changed in Stage 1B")
    assert_true(not (RESULTS / "stage1b_shared_pip_nj.nwk").exists(), "NJ tree must not exist when dXY was not calculated")
    assert_true(
        not (ROOT / "empirical" / "songbirds_pips" / "figures" / "songbirds_stage1b_shared_pip_nj.pdf").exists(),
        "NJ PDF figure must not exist when dXY was not calculated",
    )
    assert_true(set(genotypes) == {"AA", "BB"}, "genotypes must be AA/BB only")
    assert_true(len(set(sample_ids)) == 22, "sample IDs must be unique")
    assert_true(len(set(species)) == 11, "species names must include 11 unique taxa")

    matrix_names = [field for field in matrix[0].keys() if field != "sample"]
    assert_true(matrix_names == sample_ids, "dXY matrix columns must match frozen sample order")
    assert_true([row["sample"] for row in matrix] == sample_ids, "dXY matrix rows must match frozen sample order")
    for i, row in enumerate(matrix):
        for j, col in enumerate(sample_ids):
            value = row[col]
            if i == j:
                assert_true(value == "0", "dXY matrix diagonal must be zero")
            else:
                assert_true(value == "NA", "blocked dXY matrix off-diagonal values must be NA")
                assert_true(matrix[j][sample_ids[i]] == value, "dXY matrix must be symmetric")

    print("Stage 1B validation passed: frozen design valid; blocked outputs are internally consistent.")


if __name__ == "__main__":
    main()
