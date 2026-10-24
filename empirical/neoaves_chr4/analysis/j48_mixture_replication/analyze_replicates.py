#!/usr/bin/env python3
"""Postprocess returned replicate trees; missing results remain pending."""
from __future__ import annotations
import argparse, csv, importlib.util, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / "empirical" / "neoaves_chr4"
SCRIPT = BASE / "scripts" / "05_stage4e_analyze_j48_mixtures.py"
spec = importlib.util.spec_from_file_location("j48_original", SCRIPT); mod = importlib.util.module_from_spec(spec); sys.modules[spec.name] = mod; spec.loader.exec_module(mod)

def read_tsv(path):
    with path.open(newline="") as f: return list(csv.DictReader(f, delimiter="\t"))

def load_records(path):
    return mod.parse_tree(path, "CULength")

def analyze(result_root: Path, manifest: Path) -> tuple[list[dict], list[dict]]:
    original = {r["treatment"]: r for r in mod.read_tsv(mod.RESULTS / "j48_mixture_topology_summary.tsv")}
    t0_taxa, t0 = mod.parse_tree(mod.CLUSTER / mod.T0 / "tree.cu.nwk", "CULength")
    _, non = mod.parse_tree(mod.CLUSTER / mod.NON / "tree.cu.nwk", "CULength")
    rows, summary = [], []
    t0_splits, non_splits = set(t0), set(non); groups = mod.load_frozen_groups()
    for row in read_tsv(manifest):
        d = result_root / Path(row["result_relative_path"]); cu_path = d / "tree.cu.nwk"
        if not cu_path.exists():
            rows.append({"mixture_id": row["mixture_id"], "epsilon": row["epsilon_requested"], "replicate": row["replicate"], "status": "pending"}); continue
        taxa, rec = load_records(cu_path); splits = set(rec)
        focal = mod.focal_states(rec, taxa, groups)
        rows.append({"mixture_id": row["mixture_id"], "epsilon": row["epsilon_requested"], "replicate": row["replicate"], "status": "complete", "rf_vs_T0": len(splits ^ t0_splits), "rf_vs_NONCHR4": len(splits ^ non_splits), "topology_identity": mod.topology_identity(splits, t0_splits, non_splits), "equals_T0": str(splits == t0_splits).lower(), "equals_NONCHR4": str(splits == non_splits).lower(), "novel_splits": ";".join(sorted(splits - t0_splits - non_splits)), "n_novel_splits": len(splits - t0_splits - non_splits), "focal_states": json.dumps({k: v.get("state") for k,v in focal.items()}, sort_keys=True), "focal_splits": json.dumps({k: v.get("split_key", "") for k,v in focal.items()}, sort_keys=True)})
    for epsilon in sorted({r["epsilon"] for r in rows}, key=float):
        complete = [r for r in rows if r["epsilon"] == epsilon and r["status"] == "complete"]
        n = len(complete)
        for label, key in (("T0", "equals_T0"), ("NONCHR4", "equals_NONCHR4")):
            summary.append({"epsilon": epsilon, "category": "complete_topology", "state": label, "n_complete": n, "n_recoveries": sum(r[key] == "true" for r in complete), "fraction": (sum(r[key] == "true" for r in complete) / n if n else "")})
        for relationship in ("N61", "N62", "Columbea"):
            states = [json.loads(r["focal_states"]).get(relationship, "") for r in complete]
            for state in sorted(set(states)):
                summary.append({"epsilon": epsilon, "category": relationship, "state": state, "n_complete": n, "n_recoveries": states.count(state), "fraction": (states.count(state) / n if n else "")})
    return rows, summary

def main():
    p = argparse.ArgumentParser(); p.add_argument("--result-root", type=Path, required=True); p.add_argument("--manifest", type=Path, required=True); p.add_argument("--output", type=Path, required=True); args = p.parse_args(); rows, summary = analyze(args.result_root, args.manifest); args.output.mkdir(parents=True, exist_ok=True)
    for name, data in (("replicate_topology_rows.tsv", rows), ("replicate_recovery_summary.tsv", summary)):
        fields = sorted({k for r in data for k in r}) if data else ["status"]
        with (args.output / name).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(data)
if __name__ == "__main__": main()
