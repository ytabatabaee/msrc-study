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

