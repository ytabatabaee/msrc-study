#!/usr/bin/env python3
import csv, importlib.util, sys, tempfile, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / "empirical" / "neoaves_chr4"
sys.path.insert(0, str(BASE / "scripts"))

def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path); module = importlib.util.module_from_spec(spec); sys.modules[name] = module; spec.loader.exec_module(module); return module

prep = load(BASE / "analysis/j48_mixture_replication/prepare_manifest.py", "rep_prepare")
original = load(BASE / "scripts/05_stage4e_analyze_j48_mixtures.py", "rep_original")
scorer = load(BASE / "analysis/j48_mixture_replication/analyze_fixed_topology_scores.py", "rep_scorer")
rep_analysis = load(BASE / "analysis/j48_mixture_replication/analyze_replicates.py", "rep_analysis")

class ReplicationWorkflowTest(unittest.TestCase):
    def test_frozen_manifest_shape_and_seeds(self):
        rows = prep.read_manifest(BASE / "results/stage4e/j48_mixture_replication_manifest.tsv")
        self.assertEqual(len(rows), 25); self.assertEqual(len({r["seed"] for r in rows}), 25)
        self.assertEqual({r["epsilon_requested"] for r in rows}, {"0.020", "0.050", "0.075", "0.100", "0.150"})
        self.assertEqual(sum(int(r["n_chr4_records"]) for r in rows if r["epsilon_requested"] == "0.100") // 5, 6294)
        self.assertEqual(rows[0]["mixture_id"], "J48_MIX_REP_E0p020_R01")
        self.assertEqual(rows[24]["mixture_id"], "J48_MIX_REP_E0p150_R05")

    def test_array_script_is_manifest_indexed_and_restart_safe(self):
        script = (BASE / "analysis/j48_mixture_replication/cluster/run_replicates.sbatch").read_text()
        self.assertIn("--array=0-24%2", script)
        self.assertIn("SLURM_ARRAY_TASK_ID + 2", script)
        self.assertIn("--mode topology_su", script)
        self.assertIn("--mem=14G", script)
        self.assertIn("refusing overwrite", script)
        self.assertIn("--export=ALL,PYTHON_BIN", (BASE / "analysis/j48_mixture_replication/cluster/runbook.md").read_text())

    def test_manifest_checksum_and_partial_input_guards_are_present(self):
        source = (BASE / "analysis/j48_mixture_replication/prepare_manifest.py").read_text()
        self.assertIn("validate_sources(args.manifest", source)
        self.assertIn("source checksum differs from frozen manifest", source)
        self.assertIn("refusing to resume over incomplete mixture input", source)
        validator = (BASE / "analysis/j48_mixture_replication/validate_results.py").read_text()
        self.assertIn("input checksum mismatch", validator)

    def test_original_rf_results(self):
        rows = original.read_tsv(BASE / "results/stage4e/j48_mixture_topology_summary.tsv")
        expected = {"0.00": (0, 6), "0.02": (4, 6), "0.05": (2, 4), "0.10": (6, 0), "0.15": (6, 0), "0.20": (6, 0), "0.25": (6, 0), "0.30": (6, 0), "0.40": (6, 0), "0.50": (6, 0), "0.60": (6, 0), "0.75": (6, 0), "1.00": (6, 0)}
        found = {f'{float(r["epsilon"]):.2f}': (int(r["rf_vs_NONCHR4"]), int(r["rf_vs_T0"])) for r in rows}
        self.assertEqual(found, expected)

    def test_low_localpp_audit(self):
        with (BASE / "results/stage4e/j48_mixture_localpp_audit.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        row = next(r for r in rows if r["treatment"] == "J48_MIX_E0p10" and r["relationship"] == "N62")
        self.assertEqual(row["biological_focal_state"], "q1"); self.assertEqual(row["branch_local_q_winner"], "q3"); self.assertEqual(row["branch_local_q1_status"], "not_favored"); self.assertGreater(float(row["q3"]), float(row["q1"]))
        self.assertAlmostEqual(float(row["localPP"]), float(row["pp1"]), places=8)
        self.assertEqual(row["cu_su_q_equal"], "true"); self.assertEqual(row["cu_su_localPP_equal"], "true"); self.assertEqual(row["tree_equals_T0"], "true")

    def test_biological_and_branch_local_labels_are_separate(self):
        with (BASE / "results/stage4e/j48_mixture_localpp_audit.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        row = next(r for r in rows if r["treatment"] == "J48_MIX_E0p00" and r["relationship"] == "Columbea")
        self.assertEqual(row["biological_focal_state"], "q2"); self.assertEqual(row["branch_local_q_winner"], "q1"); self.assertEqual(row["branch_local_q1_status"], "favored"); self.assertEqual(row["low_localPP"], "false")

    def test_fixed_score_reference_and_optimized_agreement(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory); t0 = directory / "t0.nwk"; non = directory / "non.nwk"
            t0.write_text("((A,B),(C,D));\n"); non.write_text("((A,C),(B,D));\n")
            lines = ["((A,B),(C,D));\n", "((A,C),(B,D));\n", "((A,B,C),D);\n", "((A,B),(C,E));\n"]
            result = scorer.score_gene_records(lines, t0, non)
            reference = scorer.score_gene_records_bruteforce(lines, t0, non)
            self.assertEqual(result["n_diagnostic_quartets"], 1); self.assertEqual(result["Q_T0"], 1); self.assertEqual(result["Q_NONCHR4"], 1); self.assertEqual(result["eligible_quartet_records"], 2); self.assertEqual(result["unresolved_quartet_records"], 2)
            self.assertEqual((result["Q_T0"], result["Q_NONCHR4"], result["delta_Q"]), (reference["Q_T0"], reference["Q_NONCHR4"], reference["delta_Q"]))
            swapped = scorer.score_gene_records([lines[0]], non, t0)
            self.assertEqual(swapped["delta_Q"], -1)
            equal = scorer.score_gene_records(lines, t0, t0)
            self.assertEqual(equal["n_diagnostic_quartets"], 0); self.assertEqual(equal["delta_Q"], 0)

    def test_completed_replicate_aggregation_and_focal_states(self):
        root = BASE / "results/stage4e/J48_MIX_REPLICATES"
        manifest = BASE / "results/stage4e/j48_mixture_replication_manifest.tsv"
        rows, summary, focal, branches, meta = rep_analysis.analyze(root, manifest)
        self.assertEqual(len(rows), 25)
        counts = {(float(r["epsilon"]), r["state"]): r["n"] for r in summary if r["category"] == "complete_topology"}
        self.assertEqual([counts[(e, "J48_T0")] for e in (0.02, 0.05, 0.075, 0.10, 0.15)], [0, 0, 5, 4, 4])
        self.assertEqual(len(meta["novel"]), 3)
        self.assertEqual({r["relationship"] for r in focal}, {"N61", "N62", "Columbea"})
        self.assertTrue(all(r["homologous_cu_su"] in {"true", "false"} for r in branches))

if __name__ == "__main__": unittest.main()
