#!/usr/bin/env python3
"""Fail-closed Bridges2 orchestration for the Mai30 workflow."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path

from mai30_workflow import (
    EXPECTED_SIMULATOR_COMMIT, HERE, MANIFEST, backward_dir, forward_dir,
    history_path, output_root, read_gzip_json, rows, selected_row,
    sha256_file, task_manifest, validate_history_payload,
)

ADAPTERS = {
    "forward": "run_forward_histories.py", "benchmark_forward": "run_forward_histories.py",
    "persistence": "run_forward_persistence.py", "panel": "select_quartet_panel.py",
    "backward": "run_backward_genealogies.py", "benchmark_backward": "run_backward_genealogies.py",
    "analysis": "analyze_quartets.py", "validation": "validate_quartet_reversal.py",
    "astral": "run_optional_astral.py",
}
REQUIRED = {"stage", "task_id", "condition", "replicate", "scenario", "rho", "history_count", "seed_policy", "status"}


def check_manifest() -> None:
    data = rows()
    if not data or set(data[0]) != REQUIRED:
        raise SystemExit("manifest is empty or has an unexpected header")
    allowed = {"0", "0.0001", "0.001", "0.002", "0.01", "0.02", ""}
    for row in data:
        if row["rho"] not in allowed:
            raise SystemExit(f"rho outside prescribed grid: {row['rho']}")
    for stage in {row["stage"] for row in data}:
        ids = [int(row["task_id"]) for row in data if row["stage"] == stage]
        if ids != list(range(len(ids))):
            raise SystemExit(f"task ids for {stage} are not contiguous")


def live_compute_allocation() -> None:
    required = ("SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_JOB_PARTITION", "SLURM_CPUS_PER_TASK")
    missing = [name for name in required if not os.environ.get(name)]
    if missing:
        raise SystemExit("refusing computational stage outside Slurm: missing " + ", ".join(missing))
    scontrol = shutil.which("scontrol")
    if not scontrol:
        raise SystemExit("refusing computational stage: scontrol is unavailable")
    try:
        info = subprocess.check_output([scontrol, "show", "job", os.environ["SLURM_JOB_ID"]], text=True, stderr=subprocess.STDOUT)
        hosts = subprocess.check_output([scontrol, "show", "hostnames", os.environ["SLURM_JOB_NODELIST"]], text=True, stderr=subprocess.STDOUT).split()
    except (OSError, subprocess.CalledProcessError) as exc:
        raise SystemExit("refusing computational stage: Slurm could not verify the active job") from exc
    if "JobState=RUNNING" not in info:
        raise SystemExit("refusing computational stage: Slurm job is not RUNNING")
    host = socket.gethostname().split(".")[0]
    if host not in {item.split(".")[0] for item in hosts}:
        raise SystemExit(f"refusing computational stage: host {host} is outside allocated nodes")


def simulator_root() -> Path:
    value = os.environ.get("MSRC_SIM_ROOT")
    if not value:
        raise SystemExit("MSRC_SIM_ROOT must point to the pinned msrc-sim checkout")
    root = Path(value).expanduser().resolve()
    if not (root / ".git").exists():
        raise SystemExit(f"MSRC_SIM_ROOT is not a git checkout: {root}")
    actual = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    if actual != EXPECTED_SIMULATOR_COMMIT:
        raise SystemExit(f"msrc-sim commit mismatch: expected {EXPECTED_SIMULATOR_COMMIT}, found {actual}")
    return root


def preflight() -> None:
    check_manifest()
    if not shutil.which("git") or not shutil.which("sbatch"):
        raise SystemExit("git and sbatch are required")
    python = os.environ.get("MSRC_PYTHON")
    if not python or not Path(python).exists():
        raise SystemExit("MSRC_PYTHON must name the configured Bridges2 interpreter")
    root = simulator_root()
    missing = [name for name in ADAPTERS.values() if not (HERE / name).exists()]
    if missing:
        raise SystemExit("missing adapter entry points: " + ", ".join(sorted(set(missing))))
    print(json.dumps({"status": "preflight_ok", "simulator_root": str(root), "simulator_commit": EXPECTED_SIMULATOR_COMMIT}, indent=2))


def prerequisites(stage: str, row: dict[str, str], root: Path) -> None:
    if stage in {"backward", "benchmark_backward"}:
        required_stage = "forward" if stage == "backward" else "benchmark_forward"
        forward_row = next(r for r in rows() if r["stage"] == required_stage and r["condition"] == row["condition"] and r["replicate"] == row["replicate"] and r["scenario"] == row["scenario"])
        for history_id in range(int(row["history_count"])):
            path = history_path(forward_row, history_id, root)
            if not path.exists():
                raise SystemExit(f"missing forward checkpoint: {path}")
            validate_history_payload(read_gzip_json(path), forward_row, history_id)
    if stage in {"analysis", "validation", "astral"} and not (root / "tasks/panel/0/output_manifest.json").exists():
        raise SystemExit("validated quartet panel is required")
    if stage == "validation" and not (root / "tasks/analysis/0/output_manifest.json").exists():
        raise SystemExit("validated quartet analysis is required")
    if stage == "astral":
        for i in range(48):
            row_i = selected_row("backward", i)
            if not (backward_dir(row_i, root) / "genealogies.tsv").exists():
                raise SystemExit(f"validated backward collection is required for ASTRAL: task {i}")


def validate_output_manifest(stage: str, task_id: int, root: Path) -> dict:
    path = task_manifest(stage, task_id, root)
    if not path.exists():
        raise SystemExit(f"adapter did not write output manifest: {path}")
    try:
        value = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise SystemExit(f"invalid output manifest: {path}") from exc
    if value.get("status") != "complete" or value.get("stage") != stage or value.get("task_id") != task_id:
        raise SystemExit(f"output manifest identity/status invalid: {path}")
    if value.get("simulator_commit") != EXPECTED_SIMULATOR_COMMIT:
        raise SystemExit("output manifest simulator commit mismatch")
    if value.get("expected_count") != value.get("actual_count"):
        raise SystemExit("output manifest count validation failed")
    ids = value.get("ids")
    if ids is not None and len(ids) != len(set(ids)):
        raise SystemExit("output manifest contains duplicate identifiers")
    for item in value.get("files", []):
        file_path = root / item["path"]
        if not file_path.exists() or file_path.stat().st_size == 0:
            raise SystemExit(f"output manifest file missing or empty: {file_path}")
        if file_path.stat().st_size != item["size"] or sha256_file(file_path) != item["sha256"]:
            raise SystemExit(f"output manifest checksum mismatch: {file_path}")
    if not value.get("files"):
        raise SystemExit("output manifest has no files")
    return value


def marker_path(root: Path, stage: str, task_id: int) -> Path:
    return root / "markers" / f"{stage}_{task_id}.done.json"


def execute(stage: str, mode: str, task_id: int | None, root: Path) -> None:
    live_compute_allocation()
    preflight()
    if task_id is None:
        if len([r for r in rows() if r["stage"] == stage]) != 1:
            raise SystemExit(f"--task-id is required for {stage}")
        task_id = 0
    row = selected_row(stage, task_id)
    prerequisites(stage, row, root)
    marker = marker_path(root, stage, task_id)
    if marker.exists():
        validate_output_manifest(stage, task_id, root)
        print(f"already complete and revalidated: {marker}")
        return
    python = os.environ.get("MSRC_PYTHON", sys.executable)
    env = os.environ.copy()
    env.update({"MAI_STAGE": stage, "MAI_TASK_ID": str(task_id), "MAI_TASK_JSON": json.dumps(row), "MAI_OUTPUT_DIR": str(root), "MAI_BASE_SEED": "20261022"})
    started = time.monotonic()
    result = subprocess.run([python, str(HERE / ADAPTERS[stage]), "--task-id", str(task_id), "--output-dir", str(root)], cwd=str(HERE), env=env)
    if result.returncode:
        raise SystemExit(result.returncode)
    manifest = validate_output_manifest(stage, task_id, root)
    marker.parent.mkdir(parents=True, exist_ok=True)
    metadata = {"stage": stage, "mode": mode, "task_id": task_id, "manifest": str(MANIFEST), "simulator_commit": EXPECTED_SIMULATOR_COMMIT, "output_manifest": str(task_manifest(stage, task_id, root).relative_to(root)), "output_manifest_sha256": sha256_file(task_manifest(stage, task_id, root)), "wall_seconds": time.monotonic() - started, "actual_count": manifest["actual_count"]}
    temporary = marker.with_name(f".{marker.name}.tmp")
    temporary.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, marker)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("prepare-only", "preflight", "run", "analyze-existing"), default="prepare-only")
    parser.add_argument("--stage", choices=tuple(ADAPTERS), default="forward")
    parser.add_argument("--task-id", type=int)
    parser.add_argument("--output-dir", type=Path, default=output_root())
    args = parser.parse_args()
    check_manifest()
    if args.mode == "prepare-only":
        print(json.dumps({"status": "prepared", "manifest": str(MANIFEST), "rows": len(rows()), "simulator_imported": False}, indent=2))
    elif args.mode == "preflight":
        live_compute_allocation()
        preflight()
    else:
        execute(args.stage, args.mode, args.task_id, args.output_dir.resolve())


if __name__ == "__main__":
    main()
