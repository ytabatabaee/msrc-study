#!/usr/bin/env python3
"""Dry-run-safe ASTER/ASTRAL4 runner for Stage 4E CU/SU jobs.

The upstream ASTER/ASTRAL-IV interface defines ``CULength`` as coalescent
units and ``SULength`` as CASTLES substitution-per-site units. Full-scale
execution is permitted only inside a SLURM job; local execution is limited to
the synthetic preflight smoke test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def astra_command(aster_bin: Path, input_path: Path, output_path: Path, threads: int,
                  length: str = "CULength", genelen: float = 1000.0,
                  sulengthtype: int = 1, constraint: Path | None = None) -> list[str]:
    if length not in {"CULength", "SULength"}:
        raise ValueError("length must be CULength or SULength")
    command = [str(aster_bin), "-i", str(input_path), "-o", str(output_path),
               "-u", "2", "-t", str(threads), "--length", length,
               "--genelength", str(genelen), "--sulengthtype", str(sulengthtype)]
    if constraint is not None:
        command.extend(["-C", "-c", str(constraint)])
    return command


def running_under_slurm() -> bool:
    return bool(os.environ.get("SLURM_JOB_ID"))


def recognized_cluster_hostname() -> bool:
    host = platform.node().lower()
    return any(token in host for token in ("bridges", "psc.edu", "login", "compute"))


def require_cluster_execution(execute: bool, allow_cluster_execution: bool) -> None:
    if not execute:
        return
    if running_under_slurm():
        return
    if allow_cluster_execution and recognized_cluster_hostname():
        return
    raise RuntimeError("Large Stage-4E ASTER analysis may only be executed inside a SLURM job.")


def help_text(aster_bin: Path) -> str:
    result = subprocess.run([str(aster_bin), "-h"], capture_output=True, text=True, check=False)
    return (result.stdout + result.stderr).strip()


def validate_four_taxon_newick(path: Path) -> None:
    import dendropy
    text = path.read_text().strip()
    if not text or not text.endswith(";"):
        raise RuntimeError(f"Invalid smoke-test output: {path}")
    tree = dendropy.Tree.get(data=text, schema="newick", rooting="force-unrooted", preserve_underscores=True)
    taxa = {node.taxon.label for node in tree.leaf_node_iter()}
    if taxa != {"A", "B", "C", "D"}:
        raise RuntimeError(f"Smoke-test output has unexpected taxa: {sorted(taxa)}")
    for edge in tree.postorder_edge_iter():
        if edge.length is not None and (edge.length < 0 or not float(edge.length) == float(edge.length)):
            raise RuntimeError(f"Nonfinite or negative branch length in {path}")


def run_smoke_test(aster_bin: Path, output_root: Path, mode: str = "topology_cu") -> dict[str, object]:
    smoke_dir = output_root / "aster_smoke_test"
    smoke_dir.mkdir(parents=True, exist_ok=True)
    tiny_input = smoke_dir / "tiny_input.tre"
    tiny_input.write_text(
        "((A,B),(C,D));\n((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6):0.7;\n((A:0.1,C:0.2):0.3,(B:0.4,D:0.5):0.6):0.7;\n"
        if mode == "both" else "((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n"
    )
    suffixes = {"topology_cu": ["cu"], "topology_su": ["su"], "both": ["cu", "su"]}[mode]
    stdout_path, stderr_path = smoke_dir / "stdout.log", smoke_dir / "stderr.log"
    commands, results = [], []
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        for suffix in suffixes:
            output = smoke_dir / f"tiny_output.{suffix}.nwk"
            length = "CULength" if suffix == "cu" else "SULength"
            cmd = astra_command(aster_bin, tiny_input, output, 1, length=length)
            commands.append(cmd)
            results.append(subprocess.run(cmd, stdout=stdout, stderr=stderr, check=False))
    passed, error = False, ""
    try:
        for suffix, result in zip(suffixes, results):
            if result.returncode != 0:
                raise RuntimeError(f"Smoke-test return code was {result.returncode}")
            validate_four_taxon_newick(smoke_dir / f"tiny_output.{suffix}.nwk")
        passed = True
    except Exception as exc:
        error = str(exc)
    return {"smoke_test_mode": mode,
            # Keep the historical single-command field stable for existing
            # preflight consumers; the complete command list is recorded below.
            "smoke_test_command": commands[0] if commands else [],
            "smoke_test_commands": commands,
            "smoke_test_returncode": max((r.returncode for r in results), default=1),
            "smoke_test_passed": passed, "smoke_test_tree": str(smoke_dir),
            "smoke_test_stdout": str(stdout_path), "smoke_test_stderr": str(stderr_path),
            "smoke_test_error": error}


def run_preflight(args: argparse.Namespace) -> None:
    if not args.aster_bin.exists() or not os.access(args.aster_bin, os.X_OK):
        raise RuntimeError(f"ASTER_BIN is missing or not executable: {args.aster_bin}")
    args.output_root.mkdir(parents=True, exist_ok=True)
    usage = shutil.disk_usage(args.output_root)
    rows = []
    with args.input_manifest.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            path = args.data_dir / Path(row["relative_path"]).name
            if not path.exists():
                raise RuntimeError(f"Missing input: {path}")
            if row["sha256"] and sha256(path) != row["sha256"]:
                raise RuntimeError(f"Input checksum mismatch: {path}")
            rows.append(row)
    smoke = run_smoke_test(args.aster_bin, args.output_root, mode="both")
    report = {"status": "preflight_passed" if smoke["smoke_test_passed"] else "preflight_failed",
              "large_inference_run": False, "aster_bin": str(args.aster_bin.resolve()),
              "aster_bin_sha256": sha256(args.aster_bin), "aster_help": help_text(args.aster_bin),
              "verified_length_modes": {"topology_cu": "CULength", "topology_su": "SULength"},
              "sulengthtype": 1, "default_genelength": 1000.0,
              "hostname": platform.node(), "date_utc": datetime.now(timezone.utc).isoformat(),
              "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
              "slurm_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
              "memory_request": os.environ.get("SLURM_MEM_PER_NODE") or os.environ.get("SLURM_MEM_PER_CPU"),
              "output_root": str(args.output_root), "output_free_bytes": usage.free,
              "n_inputs_checked": len(rows), **smoke}
    out = args.output_root / "stage4e_preflight_runtime.json"
    out.write_text(json.dumps(report, indent=2) + "\n")
    if not smoke["smoke_test_passed"]:
        raise RuntimeError(f"ASTER smoke test failed; see {smoke['smoke_test_stderr']}: {smoke['smoke_test_error']}")
    print(json.dumps({"status": "preflight_passed", "report": str(out)}, indent=2))


def run_one(args: argparse.Namespace, mode: str, output_path: Path) -> dict[str, object]:
    length = "CULength" if mode == "topology_cu" else "SULength"
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK") or args.threads)
    command = astra_command(args.aster_bin, args.input, output_path, threads, length=length,
                            genelen=args.genelength, sulengthtype=args.sulengthtype,
                            constraint=args.constraint)
    started = datetime.now(timezone.utc).isoformat()
    with (args.output_dir / f"{mode}.stdout.log").open("w") as stdout, (args.output_dir / f"{mode}.stderr.log").open("w") as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    if result.returncode == 0 and (not output_path.exists() or output_path.stat().st_size == 0):
        raise RuntimeError(f"{mode} completed without an output tree")
    return {"mode": mode, "length": length, "command": command, "returncode": result.returncode,
            "output": str(output_path), "started_utc": started}


def run_treatment(args: argparse.Namespace) -> None:
    require_cluster_execution(args.execute, args.allow_cluster_execution)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    mode = getattr(args, "mode", "topology_cu")
    genelen = getattr(args, "genelength", 1000.0)
    sulengthtype = getattr(args, "sulengthtype", 1)
    constraint = getattr(args, "constraint", None)
    args.genelength, args.sulengthtype, args.constraint = genelen, sulengthtype, constraint
    modes = ["topology_cu", "topology_su"] if mode == "both" else [mode]
    outputs = {"topology_cu": args.output_dir / "tree.cu.nwk", "topology_su": args.output_dir / "tree.su.nwk"}
    commands = [astra_command(args.aster_bin, args.input, outputs[m], int(os.environ.get("SLURM_CPUS_PER_TASK") or args.threads), length=("CULength" if m == "topology_cu" else "SULength"), genelen=genelen, sulengthtype=sulengthtype, constraint=constraint) for m in modes]
    (args.output_dir / "command.txt").write_text("\n".join(" ".join(c) for c in commands) + "\n")
    provenance = {"treatment": args.treatment, "dry_run": not args.execute, "command": commands,
                  "input": str(args.input), "input_sha256": sha256(args.input),
                  "aster_bin": str(args.aster_bin.resolve()),
                  "aster_bin_sha256": sha256(args.aster_bin) if args.aster_bin.exists() else None,
                  "hostname": platform.node(), "date_utc": datetime.now(timezone.utc).isoformat(),
                  "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                  "cpu_count": int(os.environ.get("SLURM_CPUS_PER_TASK") or args.threads),
                  "memory_request": os.environ.get("SLURM_MEM_PER_NODE") or os.environ.get("SLURM_MEM_PER_CPU"),
                  "mode": mode, "genelength": genelen, "sulengthtype": sulengthtype,
                  "constraint": str(constraint) if constraint else None}
    (args.output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    if not args.execute:
        (args.output_dir / "exit_code.txt").write_text("DRY_RUN\n")
        print("DRY RUN:\n" + "\n".join(" ".join(c) for c in commands))
        return
    runs = []
    for m in modes:
        run = run_one(args, m, outputs[m]); runs.append(run)
        if run["returncode"] != 0:
            (args.output_dir / "exit_code.txt").write_text(f"{run['returncode']}\n")
            raise SystemExit(run["returncode"])
    if mode == "topology_cu":
        shutil.copyfile(outputs["topology_cu"], args.output_dir / "tree.nwk")
    (args.output_dir / "exit_code.txt").write_text("0\n")
    provenance["runs"] = runs
    (args.output_dir / "run.json").write_text(json.dumps(provenance, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--execute", action="store_true", help="Actually launch ASTER/ASTRAL4; default is dry-run.")
    parser.add_argument("--allow-cluster-execution", action="store_true")
    parser.add_argument("--aster-bin", type=Path, required=True)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--treatment")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--mode", choices=["topology_cu", "topology_su", "both"], default="both")
    parser.add_argument("--genelength", type=float, default=1000.0)
    parser.add_argument("--sulengthtype", type=int, choices=[1, 2], default=1)
    parser.add_argument("--constraint", type=Path)
    parser.add_argument("--input-manifest", type=Path)
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument("--output-root", type=Path)
    args = parser.parse_args()
    if args.preflight:
        for required in ("input_manifest", "data_dir", "output_root"):
            if getattr(args, required) is None:
                parser.error(f"--preflight requires --{required.replace('_', '-')}")
        run_preflight(args)
        return
    if args.input is None or args.output_dir is None or not args.treatment:
        parser.error("--input, --output-dir, and --treatment are required outside --preflight")
    run_treatment(args)


if __name__ == "__main__":
    main()
