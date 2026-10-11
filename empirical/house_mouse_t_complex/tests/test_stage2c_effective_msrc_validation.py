import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/build_stage2c_effective_msrc_validation.py"
SPEC = importlib.util.spec_from_file_location("stage2c", SCRIPT)
stage2c = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(stage2c)


class Stage2CEffectiveMSRCTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.aggregate = stage2c.load_aggregate()
        cls.rows = stage2c.prediction_rows(cls.aggregate)

    def test_symmetric_formula(self):
        self.assertAlmostEqual(stage2c.delta_from_q(0.6856), 0.5284, places=4)
        self.assertAlmostEqual(sum(stage2c.predicted_q(0.5284)), 1.0)

    def test_frozen_aggregate_values(self):
        self.assertAlmostEqual(self.aggregate["STT"]["q_arrangement"], 0.6856, places=4)
        self.assertAlmostEqual(self.aggregate["TST"]["q_arrangement"], 0.6937, places=4)
        self.assertAlmostEqual(self.aggregate["TTS"]["q_arrangement"], 0.7068, places=4)

    def test_leave_one_out_and_stt_only(self):
        loo = [r for r in self.rows if r["analysis"] == "leave_one_configuration_out"]
        stt = [r for r in self.rows if r["analysis"] == "STT_only_prediction"]
        self.assertEqual(len(loo), 3)
        self.assertEqual(len(stt), 2)
        self.assertTrue(all(float(r["max_abs_error"]) < 0.03 for r in loo))
        expected = {"STT": 0.014655, "TST": 0.002526, "TTS": -0.017181}
        for row in loo:
            self.assertAlmostEqual(float(row["arrangement_support_error"]), expected[row["configuration"]], places=5)

    def test_block_assignment_is_complete(self):
        blocks = stage2c.load_block_rows(stage2c.PRIMARY_BLOCK_SIZE)
        self.assertGreater(len(blocks), 5)
        self.assertTrue(all(set(v) == set(stage2c.PATTERNS) for v in blocks.values()))

    def test_bootstrap_uses_whole_blocks_and_is_deterministic(self):
        blocks = stage2c.load_block_rows(stage2c.PRIMARY_BLOCK_SIZE)
        one = stage2c.bootstrap_rows(blocks, stage2c.PRIMARY_BLOCK_SIZE, 5, 123)
        two = stage2c.bootstrap_rows(blocks, stage2c.PRIMARY_BLOCK_SIZE, 5, 123)
        self.assertEqual(one, two)
        self.assertEqual(len(one), 5 * 3)
        self.assertTrue(all(r["n_sampled_blocks"] == len(blocks) for r in one))
        for r in one:
            self.assertAlmostEqual(sum(float(r[f"q_{topo}_predicted"]) for topo in ("species", "t_alt", "other")), 1.0)
            self.assertTrue(all(float(r[f"delta_{p}"]) == float(r[f"delta_{p}"]) for p in stage2c.PATTERNS))

    def test_all_block_sizes_complete(self):
        for size in stage2c.BLOCK_SIZES:
            blocks = stage2c.load_block_rows(size)
            self.assertGreater(len(blocks), 5)
            self.assertTrue(all(set(v) == set(stage2c.PATTERNS) for v in blocks.values()))

    def test_bootstrap_summary_is_frozen(self):
        summary_path = stage2c.RESULTS / "stage2c_effective_msrc_bootstrap_summary.tsv"
        text = summary_path.read_text()
        self.assertIn("arrangement_support_error_median", text.splitlines()[0])
        self.assertIn("1000000\theldout_configuration\tSTT", text)
        self.assertIn("0.0163492324901867", text)

    def test_panel_c_uses_signed_error_and_zero_reference(self):
        source = SCRIPT.read_text()
        panel = source[source.index("    ax = axes[2]"):source.index("    for suffix in (\"png\", \"pdf\")")]
        self.assertIn("arrangement_support_error_median", panel)
        self.assertIn("arrangement_support_error\"]", panel)
        self.assertNotIn("max_abs_error_median", panel)
        self.assertIn("ax.axvline(0", panel)
        self.assertIn("Observed - predicted arrangement support", panel)

    def test_all_spearman_correlations_are_reported(self):
        rows = stage2c.local_heterogeneity_rows(stage2c.load_block_rows(stage2c.PRIMARY_BLOCK_SIZE), stage2c.PRIMARY_BLOCK_SIZE)
        stats = stage2c.local_heterogeneity_summary(rows)
        self.assertAlmostEqual(stats["spearman_STT_TST"], 0.2617, places=3)
        self.assertAlmostEqual(stats["spearman_STT_TTS"], 0.4883, places=3)
        self.assertAlmostEqual(stats["spearman_TST_TTS"], 0.7859, places=3)
        report = (stage2c.RESULTS / "stage2c_effective_msrc_report.md").read_text()
        self.assertIn("TST/TTS", report)

    def test_local_delta_is_not_truncated_and_upstream_is_not_called(self):
        rows = stage2c.local_heterogeneity_rows(stage2c.load_block_rows(stage2c.PRIMARY_BLOCK_SIZE), stage2c.PRIMARY_BLOCK_SIZE)
        values = [float(r["delta_block"]) for r in rows]
        self.assertTrue(any(v < 0 for v in values) or any(v > 1 for v in values) or all(0 <= v <= 1 for v in values))
        self.assertTrue(any(v < 0 for v in values))
        source = SCRIPT.read_text()
        for forbidden in ("subprocess", "Popen", "os.system", "astral4", "raxml"):
            self.assertNotIn(forbidden, source.lower())


if __name__ == "__main__":
    unittest.main()
