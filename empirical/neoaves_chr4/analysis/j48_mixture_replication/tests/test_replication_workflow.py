#!/usr/bin/env python3
import csv, importlib.util, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / "empirical" / "neoaves_chr4"
sys.path.insert(0, str(BASE / "scripts"))

def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path); module = importlib.util.module_from_spec(spec); sys.modules[name] = module; spec.loader.exec_module(module); return module

prep = load(BASE / "analysis/j48_mixture_replication/prepare_manifest.py", "rep_prepare")
original = load(BASE / "scripts/05_stage4e_analyze_j48_mixtures.py", "rep_original")

class ReplicationWorkflowTest(unittest.TestCase):
    def test_frozen_manifest_shape_and_seeds(self):
        rows = prep.read_manifest(BASE / "results/stage4e/j48_mixture_replication_manifest.tsv")
        self.assertEqual(len(rows), 25); self.assertEqual(len({r["seed"] for r in rows}), 25)
        self.assertEqual({r["epsilon_requested"] for r in rows}, {"0.020", "0.050", "0.075", "0.100", "0.150"})
        self.assertEqual(sum(int(r["n_chr4_records"]) for r in rows if r["epsilon_requested"] == "0.100") // 5, 6294)

    def test_original_rf_results(self):
        rows = original.read_tsv(BASE / "results/stage4e/j48_mixture_topology_summary.tsv")
        expected = {"0.00": (0, 6), "0.02": (4, 6), "0.05": (2, 4), "0.10": (6, 0), "0.15": (6, 0), "0.20": (6, 0), "0.25": (6, 0), "0.30": (6, 0), "0.40": (6, 0), "0.50": (6, 0), "0.60": (6, 0), "0.75": (6, 0), "1.00": (6, 0)}
        found = {f'{float(r["epsilon"]):.2f}': (int(r["rf_vs_NONCHR4"]), int(r["rf_vs_T0"])) for r in rows}
        self.assertEqual(found, expected)

    def test_low_localpp_audit(self):
        with (BASE / "results/stage4e/j48_mixture_localpp_audit.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        row = next(r for r in rows if r["treatment"] == "J48_MIX_E0p10" and r["relationship"] == "N62")
        self.assertEqual(row["inferred_state"], "q1"); self.assertGreater(float(row["q3"]), float(row["q1"]))
        self.assertAlmostEqual(float(row["localPP"]), float(row["pp1"]), places=8)
        self.assertEqual(row["cu_su_q_equal"], "true"); self.assertEqual(row["cu_su_localPP_equal"], "true"); self.assertEqual(row["tree_equals_T0"], "true")

if __name__ == "__main__": unittest.main()
