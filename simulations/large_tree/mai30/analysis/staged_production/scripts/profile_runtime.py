#!/usr/bin/env python3
"""Benchmark one authoritative MSRC history without changing its process.

The script is intentionally opt-in because the slow 0.15 condition can take
many minutes per history.  It records phase timings from the existing
collection runner and writes no completion marker unless the runner validates
the history.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--task-id", type=int, required=True)
    p.add_argument("--study-root", type=Path, default=Path(__file__).resolve().parents[3])
    p.add_argument("--output-root", type=Path)
    p.add_argument("--histories", type=int, default=1)
    args = p.parse_args()
    root = args.study_root.resolve()
    command = [sys.executable, str(root / "analysis/scripts/run_production_collection.py"),
               "--manifest", str(args.manifest), "--task-id", str(args.task_id),
               "--histories", str(args.histories)]
    if args.output_root:
        command += ["--output-root", str(args.output_root)]
    start = time.perf_counter()
    completed = subprocess.run(command, cwd=root, text=True, capture_output=True, check=False)
    elapsed = time.perf_counter() - start
    record = {"task_id": args.task_id, "histories": args.histories, "elapsed_seconds": elapsed,
              "returncode": completed.returncode, "stdout_tail": completed.stdout[-2000:],
              "stderr_tail": completed.stderr[-2000:]}
    print(json.dumps(record, indent=2))
    if completed.returncode:
        raise SystemExit(completed.returncode)


if __name__ == "__main__":
    main()
