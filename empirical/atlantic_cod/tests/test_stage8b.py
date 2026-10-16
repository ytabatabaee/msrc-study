import importlib.util
import json
from pathlib import Path
import unittest

from Bio import Phylo


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "empirical/atlantic_cod/scripts/08b_baseline_tree_time_projection.py"
SPEC = importlib.util.spec_from_file_location("stage8b_under_test", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def branch_map(tree):
    out = {}
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            continue
        key = tuple(sorted(t.name for t in clade.get_terminals()))
        out[key] = clade.branch_length
    return out


class Stage8bTests(unittest.TestCase):
    def test_pruning_is_12_tip_and_deterministically_rooted(self):
        first, audit1 = MODULE.prune_to_cod_tree()
        second, audit2 = MODULE.prune_to_cod_tree()
        self.assertEqual(len(first.get_terminals()), 12)
        self.assertEqual(MODULE.rooted_signature_bio(first), MODULE.rooted_signature_bio(second))
        self.assertEqual(audit1["root_operation"], audit2["root_operation"])
        self.assertEqual(audit1["collapsed_nodes"], audit2["collapsed_nodes"])

    def test_no_root_optimization_is_used(self):
        text = SCRIPT.read_text()
        self.assertNotIn("choose_root(", text)
        self.assertNotIn("root_edge", text)

    def test_original_tree_is_ultrametric(self):
        tree, _audit = MODULE.prune_to_cod_tree()
        distances = [tree.distance(tree.root, tip) for tip in tree.get_terminals()]
        self.assertLess(max(distances) - min(distances), 1e-8)

    def test_original_branch_lengths_are_preserved(self):
        source, _audit = MODULE.prune_to_cod_tree()
        output = Phylo.read(MODULE.ORIGINAL_TREE_OUT, "newick")
        source_map = branch_map(source)
        output_map = branch_map(output)
        self.assertEqual(set(source_map), set(output_map))
        for clade in source_map:
            if source_map[clade] is None:
                self.assertIsNone(output_map[clade])
            else:
                self.assertAlmostEqual(source_map[clade], output_map[clade], places=10)

    def test_aggregate_logic_is_the_existing_stage8_logic(self):
        rows, matrices, counts = MODULE.S8.aggregate_matrices()
        frozen = MODULE.read_tsv(MODULE.DATA / "processed/stage8_pairwise_mrca_matrices.tsv")
        frozen_map = {(r["population1"], r["population2"]): r for r in frozen}
        self.assertEqual(len(rows), 66)
        for pair in matrices["all"]:
            self.assertAlmostEqual(matrices["all"][pair], float(frozen_map[pair]["median_all"]), places=10)
            self.assertAlmostEqual(matrices["outside"][pair], float(frozen_map[pair]["median_outside"]), places=10)
            self.assertEqual(counts["all"][pair], int(frozen_map[pair]["n_all"]))
            self.assertEqual(counts["outside"][pair], int(frozen_map[pair]["n_outside"]))

    def test_projection_constraints_and_topology(self):
        tree, _audit = MODULE.prune_to_cod_tree()
        original = MODULE.source_tree_to_rnode(tree)
        _rows, matrices, _counts = MODULE.S8.aggregate_matrices()
        outside = MODULE.S8.fit_node_ages(MODULE.source_tree_to_rnode(tree), matrices["outside"], "test", "outside", "original_root")
        all_fit = MODULE.S8.fit_node_ages(MODULE.source_tree_to_rnode(tree), matrices["all"], "test", "all", "original_root")
        MODULE.validate_rnode(outside.rooted_root)
        MODULE.validate_rnode(all_fit.rooted_root)
        self.assertEqual(MODULE.S8.topology_signature(original), MODULE.S8.topology_signature(outside.rooted_root))
        self.assertEqual(MODULE.S8.topology_signature(outside.rooted_root), MODULE.S8.topology_signature(all_fit.rooted_root))

    def test_stage8_outputs_match_the_frozen_manifest(self):
        manifest = json.loads(MODULE.STAGE8_MANIFEST.read_text())
        for rel, expected in manifest["output_checksums"].items():
            self.assertEqual(MODULE.sha256(ROOT / rel), expected, rel)

    def test_projected_trees_do_not_use_source_branch_lengths_as_constraints(self):
        source, _audit = MODULE.prune_to_cod_tree()
        projected = Phylo.read(MODULE.OUTSIDE_TREE_OUT, "newick")
        source_values = [v for v in branch_map(source).values() if v is not None]
        projected_values = [v for v in branch_map(projected).values() if v is not None]
        self.assertNotEqual(source_values, projected_values)


if __name__ == "__main__":
    unittest.main()
