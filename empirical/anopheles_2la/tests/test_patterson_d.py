#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import importlib.util

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from patterson_d_core import calculate_window, calculate_window_scalar, calculate_window_vectorized, jackknife, jackknife_difference  # noqa: E402

spec = importlib.util.spec_from_file_location("run_patterson_d", Path(__file__).resolve().parents[1] / "scripts" / "run_patterson_d.py")
run_patterson_d = importlib.util.module_from_spec(spec); spec.loader.exec_module(run_patterson_d)

def diploid(*rows):
    return np.asarray(rows, dtype=np.int8)

def test_limits_and_symmetry():
    groups = {"P1": [0], "P2": [1], "P3": [2]}
    d_plus = calculate_window(diploid([[0,0],[1,1],[1,1],[0,0]]), groups, [3])
    assert d_plus.d == 1.0
    d_minus = calculate_window(diploid([[1,1],[0,0],[1,1],[0,0]]), groups, [3])
    assert d_minus.d == -1.0
    equal = calculate_window(diploid([[0,0],[1,1],[1,1],[0,0]], [[1,1],[0,0],[1,1],[0,0]]), groups, [3])
    assert abs(equal.d) < 1e-12
    swapped = calculate_window(diploid([[1,1],[0,0],[1,1],[0,0]]), {"P1":[1],"P2":[0],"P3":[2]}, [3])
    assert swapped.d == 1.0

def test_ref_alt_recode_and_missing_heterozygote():
    gt = diploid([[0,1],[1,1],[1,1],[0,0]], [[-1,-1],[0,1],[1,1],[0,0]])
    recoded = 1 - gt
    recoded[gt < 0] = -1
    groups = {"P1":[0],"P2":[1],"P3":[2]}
    a = calculate_window(gt, groups, [3]); b = calculate_window(recoded, groups, [3])
    assert abs(a.d - b.d) < 1e-12 and a.n_pass == 1

def test_multiple_diploid_samples_are_counted_by_complete_genotype():
    gt = diploid([[0, 0], [1, 1], [1, 1], [0, 0], [0, 1], [1, 1]])
    result = calculate_window(gt, {"P1": [0, 4], "P2": [1], "P3": [2]}, [3])
    assert result.n_pass == 1 and result.n_informative == 1

def test_vectorized_estimator_matches_scalar_reference():
    rng = np.random.default_rng(20261024)
    gt = rng.integers(0, 2, size=(500, 12, 2), dtype=np.int8)
    gt[rng.random(gt.shape) < 0.08] = -1
    gt[3, 0, 0] = 2
    groups = {"P1": [0, 1, 2], "P2": [3, 4, 5], "P3": [6, 7, 8]}
    outgroup = [9, 10, 11]
    scalar = calculate_window_scalar(gt, groups, outgroup, min_mac=2)
    vector = calculate_window_vectorized(gt, groups, outgroup, min_mac=2)
    assert abs(scalar.abba - vector.abba) <= 1e-8
    assert abs(scalar.baba - vector.baba) <= 1e-8
    assert scalar.n_examined == vector.n_examined and scalar.n_pass == vector.n_pass
    assert scalar.n_informative == vector.n_informative and scalar.exclusions == vector.exclusions
    if scalar.d is not None: assert abs(scalar.d - vector.d) <= 1e-10

def test_exclusions_and_zero_denominator():
    gt = diploid([[0,0],[0,0],[0,0],[0,0]], [[0,2],[0,0],[0,0],[0,0]], [[0,0],[-1,-1],[0,0],[0,0]], [[0,0],[0,0],[0,0],[0,0]])
    result = calculate_window(gt, {"P1":[0],"P2":[1],"P3":[2]}, [3])
    assert result.d is None and result.exclusions["multiallelic"] == 1 and result.exclusions["insufficient_group_calls"] >= 1
    polymorphic = calculate_window(diploid([[0,0],[0,0],[1,1],[0,0],[1,1]]), {"P1":[0],"P2":[1],"P3":[2]}, [3,4])
    assert polymorphic.d is None and polymorphic.exclusions["outgroup_polymorphic"] == 1
    partial_outgroup = calculate_window(diploid([[0,0],[1,1],[1,1],[-1,0]]), {"P1":[0],"P2":[1],"P3":[2]}, [3])
    assert partial_outgroup.exclusions["outgroup_missing"] == 1

def test_ratio_of_sums_is_not_mean_of_sitewise_ratios():
    # One site has a 1:0 contribution and nine sites have 0:1 contributions. The
    # estimator is (1-9)/(1+9), rather than the mean of +1 and -1.
    abba_site = [[0, 0], [1, 1], [1, 1], [0, 0]]
    baba_site = [[1, 1], [0, 0], [1, 1], [0, 0]]
    gt = diploid(abba_site, *([baba_site] * 9))
    result = calculate_window(gt, {"P1": [0], "P2": [1], "P3": [2]}, [3])
    assert result.abba == 1.0 and result.baba == 9.0
    assert result.d == -0.8

def test_jackknife():
    result = jackknife([(10.0, 0.0, 10), (0.0, 10.0, 10), (5.0, 5.0, 10)])
    assert result["n_blocks"] == 3 and result["se"] is not None and result["pooled_d"] == 0.0
    filtered = jackknife([(10.0, 0.0, 10), (0.0, 0.0, 0), (0.0, 10.0, 10)])
    assert filtered["n_blocks"] == 2 and filtered["pooled_d"] == 0.0
    difference = jackknife_difference([(10.0, 0.0, 0.0, 10.0, 10, 10), (0.0, 10.0, 10.0, 0.0, 10, 10), (0.0, 0.0, 0.0, 0.0, 0, 0)])
    assert difference["n_blocks"] == 2 and difference["delta_d"] == 0.0

def test_jackknife_difference_handles_region_only_blocks_and_zero_delete():
    # The first block contributes only inside, the second only outside.  The
    # delete-one estimates have one undefined region each and therefore no
    # fabricated SE is returned.
    result = jackknife_difference([
        (4.0, 0.0, 0.0, 0.0, 4, 0),
        (0.0, 0.0, 0.0, 4.0, 0, 4),
    ])
    assert result["delta_d"] == 2.0
    assert result["n_blocks"] == 2
    assert result["se"] is None

def test_checkpoint_resume_key_simulation_has_no_repeatable_regions():
    windows = [
        {"window_id": "w1", "chrom": "2L", "start": "1", "end": "10"},
        {"window_id": "w2", "chrom": "2L", "start": "11", "end": "20"},
    ]
    designs = [{"test_id": "A", "P1": "G_s", "P2": "G_i", "P3": "C_s", "outgroup": "merus"}]
    all_keys = run_patterson_d.expected_result_keys(windows, designs)
    first = {k for k in all_keys if k[0] == "w1"}
    assert all_keys - first
    assert run_patterson_d.expected_result_keys(windows[:1], designs) - first == set()
    assert all_keys - all_keys == set()

def test_pilot_queries_are_exactly_eight_independent_50kb_windows():
    frozen = run_patterson_d.windows(True)
    regions = run_patterson_d.query_regions(frozen, True, 1_000_000)
    assert len(regions) == 8 and all(len(group) == 1 for group in regions)
    assert sum(int(group[0]["end"]) - int(group[0]["start"]) + 1 for group in regions) == 400_000
    assert [group[0]["window_id"] for group in regions] == [w["window_id"] for w in frozen]

def test_checkpoint_keys_skip_complete_regions_and_reject_bad_rows():
    windows = [{"window_id": "w1", "chrom": "2L", "start": "1", "end": "10"}, {"window_id": "w2", "chrom": "2L", "start": "11", "end": "20"}]
    designs = [{"test_id": "A", "P1": "G_s", "P2": "G_i", "P3": "C_s", "outgroup": "merus"}]
    complete = [{"window_id": w["window_id"], "chrom": w["chrom"], "start": w["start"], "end": w["end"], "test_id": "A", "filter_label": label, "P1": "G_s", "P2": "G_i", "P3": "C_s", "outgroup": "merus"} for w in windows[:1] for label in ("primary", "MAC2")]
    keys = run_patterson_d.validate_completed_rows(complete, windows, designs)
    assert keys == {("w1", "A", "primary"), ("w1", "A", "MAC2")}
    assert run_patterson_d.expected_result_keys([windows[0]], designs) - keys == set()
    assert run_patterson_d.expected_result_keys([windows[1]], designs) - keys
    bad = [dict(complete[0], start="2")]
    try:
        run_patterson_d.validate_completed_rows(bad, windows, designs)
    except ValueError as exc:
        assert "coordinate mismatch" in str(exc)
    else:
        raise AssertionError("incompatible checkpoint row was accepted")

def test_repository_root_and_preflight_output_are_isolated(tmp_path, monkeypatch):
    assert run_patterson_d.repository_root() == Path(__file__).resolve().parents[3]
    monkeypatch.setattr(run_patterson_d, "OUT", tmp_path)
    assert run_patterson_d.OUT / "bridges_preflight" != run_patterson_d.OUT / "pilot_results.tsv"

def test_full_report_finalization_does_not_overwrite_pilot_files(tmp_path, monkeypatch):
    monkeypatch.setattr(run_patterson_d, "OUT", tmp_path)
    frozen = {name: f"sentinel-{name}\n" for name in ("pilot_results.tsv", "pilot_report.md", "frozen_quartet_designs.tsv", "design_audit.md")}
    for name, value in frozen.items(): (tmp_path / name).write_text(value)
    (tmp_path / "inside_outside_summary.tsv").write_text("test_id\tregion\tpooled_D\tn_windows\tpooled_ABBA\tpooled_BABA\twindow_mean_D\twindow_median_D\nA_merus\tinside\t0.5\t1\t3\t1\t0.5\t0.5\n")
    (tmp_path / "jackknife_results.tsv").write_text("test_id\tregion\tblock_size_bp\tdelta_d\tse\tz\tci_low\tci_high\tn_blocks\nA_merus\tinside_minus_outside\t1000000\t0.5\t0.1\t5\t0.3\t0.7\t3\n")
    rows = [{"window_id": "w1", "eligibility": "valid", "query_seconds": "1", "compute_seconds": "2", "total_elapsed_seconds": "3", "peak_rss_mib": "10"}]
    run_patterson_d.write_full_report(rows)
    assert (tmp_path / "full_scan_report.md").exists()
    assert "Completed windows: 1" in (tmp_path / "full_scan_report.md").read_text()
    for name, value in frozen.items(): assert (tmp_path / name).read_text() == value

def test_full_output_completeness_is_enforced():
    windows = [{"window_id": "w1", "chrom": "2L", "start": "1", "end": "10"}]
    designs = [{"test_id": "A", "P1": "G_s", "P2": "G_i", "P3": "C_s", "outgroup": "merus"}]
    def row(label): return {"window_id": "w1", "chrom": "2L", "start": "1", "end": "10", "test_id": "A", "filter_label": label, "P1": "G_s", "P2": "G_i", "P3": "C_s", "outgroup": "merus"}
    try:
        run_patterson_d.validate_complete_output([row("primary")], windows, designs)
    except ValueError as exc:
        assert "incomplete full output" in str(exc)
    else:
        raise AssertionError("incomplete full output was accepted")
    run_patterson_d.validate_complete_output([row("primary"), row("MAC2")], windows, designs)

def test_slurm_scripts_propagate_python_and_ca_environment():
    full = (Path(__file__).resolve().parents[1] / "analysis" / "patterson_d" / "run_full_patterson_d.sbatch").read_text()
    preflight = (Path(__file__).resolve().parents[1] / "analysis" / "patterson_d" / "preflight_patterson_d.sbatch").read_text()
    for script in (full, preflight):
        assert 'PYTHON_BIN="${PYTHON_BIN:?set PYTHON_BIN' in script
        assert 'SSL_CERT_FILE' in script and 'REQUESTS_CA_BUNDLE' in script
        assert 'test -f "$MSRC_REPO/empirical/anopheles_2la/scripts/run_patterson_d.py"' in script
    assert "--preflight-one-chunk" in preflight
