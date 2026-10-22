#!/usr/bin/env python3
"""Static and synthetic-fixture tests; this file never imports msrc-sim."""
from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

import mai30_bridges2_driver as driver  # noqa: E402
from mai30_workflow import (  # noqa: E402
    EXPECTED_SIMULATOR_COMMIT,
    atomic_write,
    boundary_observations,
    canonical_quartet_topologies,
    induced_quartet_split,
    output_manifest,
    quartet_topology,
    seed_for,
    sampled_pattern,
    task_manifest,
    validate_history_payload,
    validate_record_states,
    validate_newick,
)


class StaticWorkflowTests(unittest.TestCase):
    def test_manifest_arrays_and_seed_uniqueness(self) -> None:
        driver.check_manifest()
        rows = driver.rows()
        self.assertEqual(sum(r["stage"] == "forward" for r in rows), 8)
        self.assertEqual(sum(r["stage"] == "backward" for r in rows), 48)
        seeds = {seed_for(tree, scenario, history) for tree in range(4) for scenario in range(2) for history in range(100)}
        self.assertEqual(len(seeds), 800)

    def test_newick_taxa_validation(self) -> None:
        validate_newick("((a,b),(c,d));", {"a", "b", "c", "d"})
        with self.assertRaises(ValueError):
            validate_newick("((a,b),(c,c));", {"a", "b", "c", "d"})

    def test_all_unrooted_quartet_topologies_and_taxon_permutations(self) -> None:
        quartet = ("a", "b", "c", "d")
        trees = {
            "((a,b),(c,d));": canonical_quartet_topologies(quartet)[0],
            "((a,c),(b,d));": canonical_quartet_topologies(quartet)[1],
            "((a,d),(b,c));": canonical_quartet_topologies(quartet)[2],
        }
        for newick, expected in trees.items():
            self.assertEqual(quartet_topology(newick, quartet), expected)
            self.assertEqual(quartet_topology(newick, ("d", "b", "a", "c")), expected)
            self.assertIsNotNone(induced_quartet_split(newick, quartet))
        self.assertIsNone(induced_quartet_split("(a,b,c,d);", quartet))

    def test_boundary_persistence_uses_counts_not_terminal_frequency(self) -> None:
        def payload(count: int, total: int, sampled: dict[str, int]) -> dict:
            record = {"branch_id": "N1", "parent_branch_id": "N2", "absolute_age": 4.0, "copy_count_A1": count, "copy_count_A0": total - count, "population_chromosomes": total, "frequency_A1": count / total, "frequency_A0": 1 - count / total, "is_branch_start": True, "is_branch_end": False}
            return {"records": [record], "taxa": sorted(sampled), "sampled": sampled}
        first = boundary_observations(payload(1, 10, {"a": 1, "b": 0}))
        second = boundary_observations(payload(0, 10, {"a": 1, "b": 0}))
        self.assertEqual(first[0]["polymorphic"], 1)
        self.assertEqual(second[0]["polymorphic"], 0)
        self.assertAlmostEqual((first[0]["polymorphic"] + second[0]["polymorphic"]) / 2, 0.5)
        self.assertEqual(sampled_pattern(payload(1, 10, {"b": 0, "a": 1})), "10")

    def test_output_manifest_checksum_validation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            old = os.environ.get("MAI_CLUSTER_OUT")
            os.environ["MAI_CLUSTER_OUT"] = str(root)
            try:
                row = {"stage": "panel", "task_id": "0", "condition": "all", "replicate": "all", "scenario": "all", "rho": "", "history_count": "0"}
                file_path = root / "tasks/panel/0/panel.tsv"
                atomic_write(file_path, b"id\n0\n")
                output_manifest(task_manifest("panel", 0, root), "panel", 0, row, [file_path], 1, 1, ["0"])
                checked = driver.validate_output_manifest("panel", 0, root)
                self.assertEqual(checked["simulator_commit"], EXPECTED_SIMULATOR_COMMIT)
                file_path.write_text("tampered\n")
                with self.assertRaises(SystemExit):
                    driver.validate_output_manifest("panel", 0, root)
            finally:
                if old is None:
                    os.environ.pop("MAI_CLUSTER_OUT", None)
                else:
                    os.environ["MAI_CLUSTER_OUT"] = old

    def test_checkpoint_record_fixture(self) -> None:
        record = {"branch_id": "N8", "parent_branch_id": "N29", "absolute_age": 4.0, "copy_count_A1": 1, "copy_count_A0": 9, "population_chromosomes": 10, "frequency_A1": 0.1, "frequency_A0": 0.9, "is_branch_start": True, "is_branch_end": False}
        validate_record_states([record])
        record["copy_count_A1"] = 2
        with self.assertRaises((ValueError, KeyError)):
            validate_record_states([record])

    def test_refuses_without_slurm(self) -> None:
        saved = {name: os.environ.pop(name, None) for name in ("SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_JOB_PARTITION", "SLURM_CPUS_PER_TASK")}
        try:
            with self.assertRaises(SystemExit):
                driver.live_compute_allocation()
        finally:
            for name, value in saved.items():
                if value is not None:
                    os.environ[name] = value


if __name__ == "__main__":
    unittest.main()
