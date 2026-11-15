import csv
import io
import math
import unittest
from pathlib import Path

from Bio import Phylo

from empirical.neurospora.scripts.plot_branch_lengths import aggregate_bins
from empirical.neurospora.scripts.run_stage1 import quartet_distances


class BranchLengthFigureTests(unittest.TestCase):
    def test_lineage_contrast_and_quartet_mean_use_patristic_distances(self):
        tree = Phylo.read(
            io.StringIO(
                "((CJ03:0.1,CJ04:0.2):0.1,(Ncra:0.1,Hisp:0.2):0.1);"
            ), "newick"
        )
        dm = quartet_distances(tree, ("CJ03", "CJ04", "Ncra", "Hisp"))
        lineage = dm["CJ03", "Hisp"] - dm["CJ04", "Hisp"]
        six_pair_mean = sum(dm[a, b] for i, a in enumerate(("CJ03", "CJ04", "Ncra", "Hisp"))
                            for b in ("CJ03", "CJ04", "Ncra", "Hisp")[i + 1:]) / 6
        self.assertAlmostEqual(lineage, -0.1)
        self.assertAlmostEqual(six_pair_mean, 0.4333333333333333)

    def test_bin_means_ignore_missing_estimates_and_retain_empty_bins(self):
        rows = [
            {"linkage_group": "LGI", "start": "10", "end": "19", "pair1_available": "1",
             "pair1_A_hisp_minus_a_hisp": "0.2", "pair1_B_mean_six_pairwise": "0.6"},
            {"linkage_group": "LGI", "start": "20", "end": "29", "pair1_available": "1",
             "pair1_A_hisp_minus_a_hisp": "NA", "pair1_B_mean_six_pairwise": "0.8"},
            {"linkage_group": "LGI", "start": "41", "end": "49", "pair1_available": "1",
             "pair1_A_hisp_minus_a_hisp": "-0.2", "pair1_B_mean_six_pairwise": "0.4"},
        ]
        bins = aggregate_bins(rows, 10)
        self.assertEqual([r["start"] for r in bins], [10, 20, 30, 40])
        self.assertAlmostEqual(bins[0]["lineage_contrast"], 0.2)
        self.assertAlmostEqual(bins[1]["B"], 0.8)
        self.assertTrue(math.isnan(bins[2]["lineage_contrast"]))
        self.assertTrue(math.isnan(bins[2]["B"]))

    def test_generated_100kb_bins_match_per_window_inventory_means(self):
        root = Path(__file__).resolve().parents[3]
        with (root / "empirical/neurospora/results/branch_length_inventory.tsv").open() as stream:
            source = list(csv.DictReader(stream, delimiter="\t"))
        with (root / "empirical/neurospora/results/lgi_100kb_branch_bins.tsv").open() as stream:
            generated = list(csv.DictReader(stream, delimiter="\t"))
        grouped = {}
        for row in source:
            if row["linkage_group"] != "LGI":
                continue
            start = int(row["start"]) // 100_000 * 100_000
            grouped.setdefault(start, {"pair1_contrast": [], "pair2_contrast": [], "B": []})
            for pair, key, field in (("pair1", "pair1_contrast", "pair1_A_hisp_minus_a_hisp"),
                                     ("pair2", "pair2_contrast", "pair2_A_hisp_minus_a_hisp")):
                if row[f"{pair}_available"] == "1":
                    value = float(row[field])
                    if math.isfinite(value):
                        grouped[start][key].append(value)
            if row["pair1_available"] == "1":
                value = float(row["pair1_B_mean_six_pairwise"])
                if math.isfinite(value):
                    grouped[start]["B"].append(value)
        for row in generated:
            start = int(row["start"])
            bin_values = grouped.get(start, {})
            pair_means = [sum(bin_values[key]) / len(bin_values[key])
                          for key in ("pair1_contrast", "pair2_contrast") if bin_values.get(key)]
            expected = {"lineage_contrast": sum(pair_means) / len(pair_means) if pair_means else math.nan,
                        "B": sum(bin_values.get("B", [])) / len(bin_values["B"]) if bin_values.get("B") else math.nan}
            for key in ("lineage_contrast", "B"):
                observed = math.nan if row[key] == "NA" else float(row[key])
                if math.isfinite(expected[key]):
                    self.assertAlmostEqual(observed, expected[key], places=11)
                else:
                    self.assertTrue(math.isnan(observed))

    def test_other_linkage_groups_are_excluded(self):
        rows = [{"linkage_group": "LGII", "start": "1", "end": "10", "pair1_available": "1",
                 "pair1_A_hisp_minus_a_hisp": "0.1", "pair1_B_mean_six_pairwise": "0.5"}]
        self.assertEqual(aggregate_bins(rows, 10), [])


if __name__ == "__main__":
    unittest.main()
