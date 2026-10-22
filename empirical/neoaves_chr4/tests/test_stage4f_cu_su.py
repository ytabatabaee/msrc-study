"""Focused regression tests for the Stage 4F CU/SU extension."""

import importlib.util
from pathlib import Path

BASE = Path(__file__).parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


runner = load("stage4f_runner", BASE / "scripts/05_stage4e_run_aster.py")
analyzer = load("stage4f_analyzer", BASE / "scripts/05_stage4e_analyze_aster.py")
annotations = load("stage4e_annotations", BASE / "scripts/stage4e_annotations.py")

T0 = BASE / "results/stage4e/cluster/T0"
T0_EXTRACTED = BASE / "results/stage4e/cluster/T0_extracted"


def test_verified_length_commands_are_explicit():
    command = runner.astra_command(Path("astral4"), Path("in.tre"), Path("out.nwk"), 4,
                                    length="SULength", genelen=2500, sulengthtype=2)
    assert command[command.index("--length") + 1] == "SULength"
    assert command[command.index("--genelength") + 1] == "2500"
    assert command[command.index("--sulengthtype") + 1] == "2"


def test_fixed_topology_command_uses_official_constraint_flags():
    command = runner.astra_command(Path("astral4"), Path("in.tre"), Path("out.nwk"), 1,
                                   constraint=Path("species.nwk"))
    assert "-C" in command
    assert command[command.index("-c") + 1] == "species.nwk"


def test_both_mode_smoke_has_one_command_per_unit(tmp_path):
    binary = tmp_path / "astral4"
    binary.write_text("""#!/bin/sh
while [ "$#" -gt 0 ]; do
  case "$1" in
    -o) output="$2"; shift 2 ;;
    *) shift ;;
  esac
done
printf '((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);\\n' > "$output"
""")
    binary.chmod(0o755)
    report = runner.run_smoke_test(binary, tmp_path / "out", mode="both")
    assert report["smoke_test_passed"]
    assert len(report["smoke_test_commands"]) == 2
    assert (tmp_path / "out/aster_smoke_test/tiny_output.cu.nwk").exists()
    assert (tmp_path / "out/aster_smoke_test/tiny_output.su.nwk").exists()


def test_branch_length_unit_prefers_explicit_annotation():
    newick = "((A:9[&CULength=0.1,SULength=0.01],B:8[&CULength=0.2,SULength=0.02]):7[&CULength=0.3,SULength=0.03],(C:6[&CULength=0.4,SULength=0.04],D:5[&CULength=0.5,SULength=0.05]):4[&CULength=0.6,SULength=0.06]);"
    _, records = analyzer.tree_records(newick)
    terminal_a = records["A"]
    assert analyzer.number(terminal_a["CULength"]) == 0.1
    assert analyzer.number(terminal_a["SULength"]) == 0.01


def test_missing_empirical_outputs_are_not_coerced_to_zero(tmp_path, monkeypatch):
    monkeypatch.setattr(analyzer, "TREE_DIR", tmp_path)
    assert analyzer.load_trees("CU") == {}
    rows = analyzer.branch_tables({}, "CU")
    assert rows[0][0]["status"] == "PENDING_CLUSTER_OUTPUT"


def annotated_fixture():
    return ("((A:0.1,B:0.2)'[CULength=1.2;SULength=0.3;localPP=1;"
            "q1=0.9;q2=0.05;q3=0.05]':0.3,(C:0.4,D:0.5)"
            "'[CULength=1.2;SULength=0.3;localPP=1;q1=0.9;q2=0.05;q3=0.05]':0.3);")


def test_single_annotated_tree_extracts_internal_lengths_and_terminal_su(tmp_path):
    source = tmp_path / "tree.su.nwk"
    source.write_text(annotated_fixture())
    summary = annotations.extract_and_write(source, tmp_path / "out", "T0")
    assert summary["taxa"] == 4
    text = (tmp_path / "out/branch_lengths.tsv").read_text()
    assert "terminal_A\tterminal\tA\tA\tNA\t0.1" in text
    assert "internal_0001\tinternal\tA|B\tinternal_0001\t1.2\t0.3" in text
    assert (tmp_path / "out/tree.annotated.nwk").read_text() == source.read_text()


def test_single_run_command_requests_only_su_annotation():
    command = runner.astra_command(Path("astral4"), Path("in.tre"), Path("tree.annotated.nwk"), 32,
                                   length="SULength", genelen=1000, sulengthtype=1)
    assert command.count("--length") == 1
    assert command[command.index("--length") + 1] == "SULength"


def test_migration_is_idempotent_and_refuses_overwrite(tmp_path):
    source = tmp_path / "legacy.su.nwk"
    source.write_text(annotated_fixture())
    out = tmp_path / "out"
    first = annotations.extract_and_write(source, out, "T0", source_kind="extracted_from_existing_su_run")
    second = annotations.extract_and_write(source, out, "T0", source_kind="extracted_from_existing_su_run")
    assert first["annotated_sha256"] == second["annotated_sha256"]
    source.write_text(source.read_text().replace("0.3);", "0.31);"))
    try:
        annotations.extract_and_write(source, out, "T0")
    except FileExistsError:
        pass
    else:
        raise AssertionError("migration overwrote a differing annotated source")


def test_committed_t0_has_363_taxa_and_360_root_invariant_internal_splits():
    _, cu_taxa, cu_records = annotations.extract_records(
        (T0 / "tree.cu.nwk").read_text(), edge_annotation="CULength")
    _, su_taxa, su_records = annotations.extract_records(
        (T0 / "tree.su.nwk").read_text(), edge_annotation="SULength")
    cu = {r["split_key"] for r in cu_records if r["branch_type"] == "internal"}
    su = {r["split_key"] for r in su_records if r["branch_type"] == "internal"}
    assert len(cu_taxa) == len(su_taxa) == 363
    assert len(cu) == len(su) == 360
    assert cu == su


def test_committed_t0_extraction_matches_original_runs_and_has_no_terminal_cu():
    _, cu_taxa, cu_records = annotations.extract_records(
        (T0 / "tree.cu.nwk").read_text(), edge_annotation="CULength")
    _, su_taxa, su_records = annotations.extract_records(
        (T0 / "tree.su.nwk").read_text(), edge_annotation="SULength")
    _, ext_cu_taxa, ext_cu_records = annotations.extract_records(
        (T0_EXTRACTED / "tree.cu.nwk").read_text(), edge_annotation="CULength")
    _, ext_su_taxa, ext_su_records = annotations.extract_records(
        (T0_EXTRACTED / "tree.su.nwk").read_text(), edge_annotation="SULength")
    cu = {r["split_key"]: r for r in cu_records if r["branch_type"] == "internal"}
    su = {r["split_key"]: r for r in su_records if r["branch_type"] == "internal"}
    ext_cu = {r["split_key"]: r for r in ext_cu_records if r["branch_type"] == "internal"}
    ext_su = {r["split_key"]: r for r in ext_su_records if r["branch_type"] == "internal"}
    assert cu_taxa == su_taxa == ext_cu_taxa == ext_su_taxa
    assert all(ext_cu[k]["CU_length"] == cu[k]["CU_length"] for k in cu)
    assert all(ext_su[k]["SU_length"] == su[k]["SU_length"] for k in su)
    assert all(r["CU_length"] is None for r in ext_cu_records if r["branch_type"] == "terminal")
    assert len({r["taxon_label"] for r in ext_su_records if r["branch_type"] == "terminal"}) == 363


def test_committed_t0_focal_branches_are_exact_and_extraction_report_passes():
    rows = (BASE / "results/stage4e/t0_focal_branches.tsv").read_text().splitlines()
    assert rows[0].startswith("focal_relationship\t")
    assert {row.split("\t", 1)[0] for row in rows[1:]} == {"Columbea", "N61", "N62"}
    states = {row.split("\t")[0]: row.split("\t")[1] for row in rows[1:]}
    assert states == {"Columbea": "q2", "N61": "q1", "N62": "q1"}
    report = (BASE / "results/stage4e/t0_full363_validation.md").read_text()
    assert "Original CU/SU unrooted RF symmetric difference: **0**" in report
    assert "| taxa_363 | PASS |" in report
    assert "| extracted_cu_terminal_lengths_missing | PASS |" in report
    assert "No ASTRAL-IV subprocess was invoked" in report
