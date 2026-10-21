import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "empirical/neoaves_chr4/scripts/11c_phlag_representative_robustness.py"
spec = importlib.util.spec_from_file_location("stage11c", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules["stage11c"] = module
spec.loader.exec_module(module)


def test_representative_sets_are_distinct_and_replicate_zero_is_frozen():
    selections = [tuple(row[2:]) for row in module.SETS]
    assert len(selections) == 10
    assert len(set(selections)) == 10
    assert selections[0] == (
        "Columba_livia",
        "Phoenicopterus_ruber",
        "Tauraco_erythrolophus",
        "Gallus_gallus",
    )


def test_species_trees_have_four_taxa_and_preserve_quartet_relationship():
    for rep in module.SETS:
        text = module.species_tree_for(rep, rep[0] == 0)
        tree = module.parse_newick(text)
        assert len(module.leaf_names(tree)) == 4
        assert len(set(module.leaf_names(tree))) == 4
        assert module.topology(tree, tuple(rep[2:])) == "q2"


def test_prepared_inputs_match_coordinate_manifests_and_original_intervals():
    out = ROOT / "empirical/neoaves_chr4/results/stage11c_phlag_robustness"
    for rep in module.SETS:
        manifest = module.read_tsv(out / "raw" / f"replicate_{rep[0]:03d}/inputs/coordinate_manifest.tsv")
        trees = (out / "raw" / f"replicate_{rep[0]:03d}/inputs/gene_trees.nwk").read_text().splitlines()
        assert len(manifest) == len(trees)
        assert all(sorted(module.leaf_names(module.parse_newick(tree))) == sorted(rep[2:]) for tree in trees)
        assert all(int(manifest[i]["start"]) <= int(manifest[i + 1]["start"]) for i in range(len(manifest) - 1))

    intervals = [(int(row["start"]), int(row["end"])) for row in module.read_tsv(module.BASE / "phlag_predicted_intervals.tsv")]
    assert intervals == [(25020000, 32300000), (32670000, 32680000), (33460000, 34480000), (44130000, 46420000), (46770000, 47540000), (47800000, 48580000), (49050000, 56820000)]


def test_phlag_output_parser_rejects_wrong_length_and_nonfinite_values(tmp_path):
    wrong = tmp_path / "wrong.txt"
    wrong.write_text("0,1\n0.5,0.5\n")
    try:
        module.parse_phlag_output(wrong, 3)
    except ValueError:
        pass
    else:
        raise AssertionError("wrong-length Phlag output was accepted")

    nonfinite = tmp_path / "nonfinite.txt"
    nonfinite.write_text("0,1\nnan,0.5\n")
    try:
        module.parse_phlag_output(nonfinite, 2)
    except ValueError:
        pass
    else:
        raise AssertionError("nonfinite Phlag posterior was accepted")


def test_gap_merging_uses_the_frozen_physical_rule():
    manifest = [{"start": "0", "end": "10"}, {"start": "200010", "end": "200020"}, {"start": "450021", "end": "450030"}]
    assert module.merge_flagged(manifest, [1, 1, 1]) == [(0, 200020), (450021, 450030)]


def test_executed_runs_and_window_predictions_have_matching_lengths():
    out = ROOT / "empirical/neoaves_chr4/results/stage11c_phlag_robustness"
    runs = module.read_tsv(out / "run_manifest.tsv")
    windows = module.read_tsv(out / "window_predictions.tsv")
    assert len(runs) == 10
    assert all(row["status"].startswith("valid") for row in runs)
    assert all(row["posterior_finite"] == "True" for row in runs)
    for run in runs:
        rid = run["replicate"]
        assert run["n_states"] == run["n_posteriors"] == run["n_loci"]
        phlag_rows = [row for row in windows if row["replicate"] == rid and row["method"] == "Phlag"]
        detector_rows = [row for row in windows if row["replicate"] == rid and row["method"] == "matched_four_taxon_detector"]
        assert len(phlag_rows) == len(detector_rows) == int(run["n_loci"])
        assert all(0.0 <= float(row["posterior"]) <= 1.0 for row in phlag_rows)
