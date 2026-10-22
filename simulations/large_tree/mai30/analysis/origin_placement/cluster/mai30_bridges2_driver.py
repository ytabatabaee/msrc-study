#!/usr/bin/env python3
"""Fail-closed Bridges2 orchestration for Mai30 origin-placement stages.

This module has no simulator imports.  The simulator-specific entry point is
provided explicitly through MAI_<STAGE>_COMMAND after the optimized checkout
has been installed on Bridges2 and verified at the required commit.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import resource
from pathlib import Path

EXPECTED_SIMULATOR_COMMIT = "bfba58d627f6f467af649cdb0b37f969349e0aab"
HERE = Path(__file__).resolve().parent
MAI = HERE.parents[2]
MANIFEST = HERE / "mai30_bridges2_manifest.tsv"
ALLOWED_RHOS = {"0", "0.0001", "0.001", "0.002", "0.01", "0.02", ""}


def rows() -> list[dict[str, str]]:
    with MANIFEST.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def selected_row(stage: str, task_id: int | None) -> dict[str, str] | None:
    matches = [row for row in rows() if row["stage"] == stage]
    if task_id is None:
        if len(matches) == 1:
            return matches[0]
        return None
    if task_id < 0 or task_id >= len(matches):
        raise SystemExit(f"task id {task_id} is outside {stage} range 0..{len(matches)-1}")
    return matches[task_id]


def check_manifest() -> None:
    required = {"stage", "task_id", "condition", "replicate", "scenario", "rho", "history_count", "seed_policy", "status"}
    data = rows()
    if not data or set(data[0]) != required:
        raise SystemExit("manifest is empty or has an unexpected header")
    for row in data:
        if row["rho"] not in ALLOWED_RHOS:
            raise SystemExit(f"rho outside prescribed grid: {row['rho']}")
    for stage in {row["stage"] for row in data}:
        ids = [int(row["task_id"]) for row in data if row["stage"] == stage]
        if ids != list(range(len(ids))):
            raise SystemExit(f"task ids for {stage} are not contiguous")


def require_slurm() -> None:
    missing = [name for name in ("SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_CPUS_PER_TASK") if not os.environ.get(name)]
    if missing:
        raise SystemExit("refusing computational stage outside a Slurm allocation; missing " + ", ".join(missing))
    if not os.environ.get("SLURM_JOB_PARTITION"):
        raise SystemExit("refusing computational stage without SLURM_JOB_PARTITION")


def simulator_root() -> Path:
    value = os.environ.get("MSRC_SIM_ROOT")
    if not value:
        raise SystemExit("MSRC_SIM_ROOT must point to the Bridges2 msrc-sim checkout")
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
        raise SystemExit("git and sbatch are required on Bridges2")
    if not os.environ.get("BRIDGES2_ACCOUNT") and not os.environ.get("SLURM_JOB_ACCOUNT"):
        raise SystemExit("set BRIDGES2_ACCOUNT; no account is encoded in this repository")
    if not os.environ.get("BRIDGES2_PARTITION") and not os.environ.get("SLURM_JOB_PARTITION"):
        raise SystemExit("set BRIDGES2_PARTITION; no partition is encoded in this repository")
    python = os.environ.get("MSRC_PYTHON")
    if not python or not Path(python).exists():
        raise SystemExit("MSRC_PYTHON must name the configured interpreter")
    root = simulator_root()
    print(json.dumps({"status": "preflight_ok", "simulator_root": str(root), "simulator_commit": EXPECTED_SIMULATOR_COMMIT}, indent=2))


def marker_path(out_dir: Path, stage: str, row: dict[str, str]) -> Path:
    key = f"{stage}_{row['task_id']}"
    return out_dir / "markers" / f"{key}.done.json"


def command_for(stage: str) -> str:
    key = "MAI_" + stage.upper() + "_COMMAND"
    value = os.environ.get(key)
    if not value:
        raise SystemExit(f"{key} is required; refusing to guess an unavailable simulator API")
    return value


def execute(stage: str, mode: str, task_id: int | None, out_dir: Path) -> None:
    require_slurm()
    preflight()
    row = selected_row(stage, task_id)
    if row is None:
        raise SystemExit(f"exactly one {stage} manifest row is required when --task-id is omitted")
    marker = marker_path(out_dir, stage, row)
    marker.parent.mkdir(parents=True, exist_ok=True)
    metadata = {"stage": stage, "mode": mode, "task_id": row["task_id"], "manifest": str(MANIFEST), "simulator_commit": EXPECTED_SIMULATOR_COMMIT}
    if marker.exists():
        saved = json.loads(marker.read_text())
        if saved == metadata:
            print(f"already complete: {marker}")
            return
        raise SystemExit(f"completion marker metadata mismatch: {marker}")
    command = command_for(stage)
    env = os.environ.copy()
    env.update({"MAI_STAGE": stage, "MAI_TASK_ID": row["task_id"], "MAI_TASK_JSON": json.dumps(row), "MAI_OUTPUT_DIR": str(out_dir), "MAI_BASE_SEED": "20261022"})
    print("executing configured stage command:", " ".join(shlex.split(command)))
    started = time.monotonic()
    result = subprocess.run(command, shell=True, cwd=str(MAI), env=env)
    if result.returncode:
        raise SystemExit(result.returncode)
    metric = {**metadata, "wall_seconds": time.monotonic() - started,
              "child_max_rss_kb": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
              "completed_history_count": row["history_count"]}
    fd, temporary = tempfile.mkstemp(prefix=f".{marker.name}.", dir=str(marker.parent), text=True)
    os.close(fd)
    try:
        Path(temporary).write_text(json.dumps(metadata, indent=2) + "\n")
        os.replace(temporary, marker)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    metric_path = marker.with_suffix(".metrics.json")
    metric_path.write_text(json.dumps(metric, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("prepare-only", "preflight", "run", "analyze-existing"), default="prepare-only")
    parser.add_argument("--stage", choices=("forward", "backward", "panel", "analysis", "validation", "astral", "benchmark_forward", "benchmark_backward"), default="forward")
    parser.add_argument("--task-id", type=int)
    parser.add_argument("--output-dir", type=Path, default=HERE / "bridges2_outputs")
    args = parser.parse_args()
    if args.mode == "prepare-only":
        check_manifest()
        print(json.dumps({"status": "prepared", "manifest": str(MANIFEST), "rows": len(rows()), "simulator_imported": False}, indent=2))
    elif args.mode == "preflight":
        require_slurm()
        preflight()
    else:
        execute(args.stage, args.mode, args.task_id, args.output_dir)


if __name__ == "__main__":
    main()
