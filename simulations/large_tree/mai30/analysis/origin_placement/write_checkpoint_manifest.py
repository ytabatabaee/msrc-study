#!/usr/bin/env python3
"""Build an aggregate manifest after interrupted or parallel checkpoint runs."""
from __future__ import annotations
import argparse, csv, gzip, json
from pathlib import Path

def main():
    p = argparse.ArgumentParser(); p.add_argument("--run-dir", type=Path, required=True); args = p.parse_args()
    rows = []
    for directory in sorted(x for x in args.run_dir.iterdir() if x.is_dir()):
        files = sorted(directory.glob("history_*.json.gz"))
        if not files: continue
        with gzip.open(files[0], "rt", encoding="utf-8") as h: first = json.load(h)
        rows.append({"condition": first["condition"], "replicate": first["replicate"], "scenario": first["scenario"],
                     "histories_requested": 100, "histories_completed": len(files),
                     "rho_count": len(first["rho_grid"]), "simulator_commit": first["simulator_commit"],
                     "status": "complete" if len(files) == 100 else "partial", "checkpoint_dir": str(directory)})
    out = args.run_dir / "forward_histories_manifest.tsv"
    with out.open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    manifest = {"status": "complete" if rows and all(r["status"] == "complete" for r in rows) else "partial",
                "histories_requested": 100, "rho_grid": [0.0, 0.0001, 0.001, 0.002, 0.01, 0.02],
                "forward_histories_reused_across_rho": True, "genealogies_per_history_per_rho": 1,
                "simulator_commit": rows[0]["simulator_commit"], "checkpoint_count": sum(r["histories_completed"] for r in rows)}
    (args.run_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
if __name__ == "__main__": main()
