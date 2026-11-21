from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from pairwise_msrc import (  # noqa: E402
    Segment, expected_pair_times, ordinary_pair_mean_time, predicted_statistics,
    transient_generator,
)


class TestPairwiseMsrc(unittest.TestCase):
    def test_zero_duration_is_matched_ordinary_coalescent(self):
        out = expected_pair_times([], ancestral_ne=37)
        self.assertEqual(out, {"AA": 74.0, "AB": 74.0, "BB": 74.0})
        stats = predicted_statistics(out, 2e-8)
        self.assertLess(abs(stats["E_Delta_AB"]), 1e-20)
        self.assertEqual(stats["R_ratio_of_expectations"], 1.0)
        self.assertLess(abs(stats["FST_ratio_of_expectations"]), 1e-15)

    def test_zero_switching_matches_closed_form_and_positive_between_excess(self):
        ne, p, rho, tau, na = 100.0, 0.2, 0.0, 140.0, 200.0
        got = expected_pair_times([Segment(tau, p, ne, rho)], ancestral_ne=na)
        floor = 1 / (2 * ne)
        ca = 1 / (2 * ne * max(1 - p, floor))
        cb = 1 / (2 * ne * max(p, floor))

        def closed(rate):
            surv = np.exp(-rate * tau)
            return -np.expm1(-rate * tau) / rate + surv * 2 * na

        self.assertTrue(np.isclose(got["AA"], closed(ca), rtol=1e-12))
        self.assertTrue(np.isclose(got["BB"], closed(cb), rtol=1e-12))
        self.assertTrue(np.isclose(got["AB"], tau + 2 * na, rtol=1e-12))
        self.assertGreater(got["AB"], (got["AA"] + got["BB"]) / 2)

    def test_balanced_symmetry_and_background_label_swap(self):
        balanced = expected_pair_times([Segment(140, 0.5, 35, 0.08)], ancestral_ne=50)
        self.assertTrue(np.isclose(balanced["AA"], balanced["BB"], rtol=1e-12))
        p2 = expected_pair_times([Segment(140, 0.23, 35, 0.08)], ancestral_ne=50)
        p3 = expected_pair_times([Segment(140, 0.77, 35, 0.08)], ancestral_ne=50)
        self.assertTrue(np.isclose(p2["AA"], p3["BB"], rtol=1e-12))
        self.assertTrue(np.isclose(p2["BB"], p3["AA"], rtol=1e-12))
        self.assertTrue(np.isclose(p2["AB"], p3["AB"], rtol=1e-12))

    def test_strong_switching_erases_initial_pair_state(self):
        got = expected_pair_times([Segment(200, 0.37, 80, 10_000)], ancestral_ne=80)
        values = np.array(list(got.values()))
        self.assertLess(values.max() - values.min(), 0.02)

    def test_piecewise_history_and_ne_scale_in_coalescent_units(self):
        seg1 = [Segment(120, 0.2, 60, 0.004), Segment(280, 0.7, 120, 0.004)]
        got1 = expected_pair_times(seg1, [(50, 80)], ancestral_ne=90)
        scale = 5.0
        seg2 = [Segment(s.duration * scale, s.frequency_b, s.effective_population_size * scale,
                        s.rho / scale) for s in seg1]
        got2 = expected_pair_times(seg2, [(50 * scale, 80 * scale)], ancestral_ne=90 * scale)
        for state in ("AA", "AB", "BB"):
            self.assertTrue(np.isclose(got2[state], scale * got1[state], rtol=1e-10))
        self.assertGreater(ordinary_pair_mean_time(seg1, [(50, 80)], ancestral_ne=90), 0)

    def test_finite_copy_boundary_matches_simulator_floor(self):
        ne = 20.0
        q = transient_generator(0.0, ne, 0.0)
        self.assertTrue(np.isclose(-q[0, 0], 1 / (2 * ne)))
        self.assertTrue(np.isclose(-q[2, 2], 1.0))
        self.assertEqual(q[0, 1], 0.0)
        self.assertTrue(np.isclose(q[1, 0], 0.0))
        with self.assertRaises(ZeroDivisionError):
            transient_generator(0.0, ne, 0.0, finite_copy_floor=False)

    def test_hudson_fst_and_relative_divergence_are_algebraically_redundant(self):
        s = predicted_statistics({"AA": 90, "AB": 150, "BB": 110}, 1e-8)
        self.assertTrue(np.isclose(s["FST_ratio_of_expectations"],
                                   1 - 1 / s["R_ratio_of_expectations"]))
        self.assertTrue(np.isclose(s["E_Delta_AB"],
                                   s["E_d_AB"] * s["FST_ratio_of_expectations"]))

    def test_matches_direct_pair_event_times_from_authoritative_msrc_generator(self):
        msrc = ROOT.parent / "msrc-sim"
        self.assertTrue(msrc.is_dir(), f"authoritative simulator expected at {msrc}")
        sys.path.insert(0, str(msrc))
        from msrcsim.species_tree import SpeciesTree
        from msrcsim.structured_coalescent import simulate_genealogy

        class FixedHistory:
            structured_origin_age = 12.0

            def frequency_at(self, branch_id, age):
                return 0.3

            def next_frequency_boundary(self, branch_id, age):
                return float("inf")

        tau, ne, ne_ancestor, p, rho = 12.0, 20, 30, 0.3, 0.01
        tree = SpeciesTree(
            f"((1:0,2:0)NODE:{tau})ROOT;", ne, 1.0,
            {"NODE": {"effective_population_size": ne},
             "ROOT": {"effective_population_size": ne_ancestor}},
        )
        history = FixedHistory()
        exact = expected_pair_times([Segment(tau, p, ne, rho)], ancestral_ne=ne_ancestor)
        rng = np.random.default_rng(19121)
        for state, arr in (("AA", (0, 0)), ("AB", (0, 1)), ("BB", (1, 1))):
            with self.subTest(pair_state=state):
                values = []
                for locus in range(2500):
                    result = simulate_genealogy(
                        locus, tree, history, dict(zip(tree.taxa, arr)), rho, 1.0, rng,
                        record_events=False,
                    )
                    values.append(result.coalescence_times[0])
                values = np.asarray(values)
                se = values.std(ddof=1) / np.sqrt(len(values))
                self.assertLess(abs(values.mean() - exact[state]), 4.5 * se + 0.1)


if __name__ == "__main__":
    unittest.main()
