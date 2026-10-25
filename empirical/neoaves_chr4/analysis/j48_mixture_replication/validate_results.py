#!/usr/bin/env python3
"""Validate returned J48 replicate outputs without launching ASTRAL."""
from __future__ import annotations
import argparse, csv, hashlib, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / "empirical" / "neoaves_chr4"
sys.path.insert(0, str(BASE / "scripts"))
from stage4e_annotations import extract_records  # noqa: E402

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""): h.update(b)
    return h.hexdigest()

def rows(path):
    with path.open(newline="") as f: return list(csv.DictReader(f, delimiter="\t"))

def validate(root: Path, manifest: Path, only: str | None = None) -> list[dict[str, object]]:
    manifest_rows = rows(manifest)
    if only is None:
        if len(manifest_rows) != 25 or len({r["mixture_id"] for r in manifest_rows}) != 25:
            raise ValueError(f"expected 25 unique manifest rows, observed {len(manifest_rows)}")
        if sorted({float(r["epsilon_requested"]) for r in manifest_rows}) != [0.02, 0.05, 0.075, 0.1, 0.15]:
            raise ValueError("manifest epsilon grid is not the five-replicate grid")
        if len({r["seed"] for r in manifest_rows}) != 25:
            raise ValueError("replicate seeds are not unique")
    result = []
    for row in manifest_rows:
        if only and row["mixture_id"] != only: continue
        d = root / Path(row["result_relative_path"])
        if not d.is_dir(): raise ValueError(f"missing result directory: {d}")
        if (d / "exit_code.txt").read_text().strip() != "0": raise ValueError(f"nonzero exit: {d}")
        cu_tree = d / "tree.cu.nwk"; su_tree = d / "tree.su.nwk"
        cu, taxa_cu, rec_cu = extract_records(cu_tree.read_text(), edge_annotation="CULength")
        su, taxa_su, rec_su = extract_records(su_tree.read_text(), edge_annotation="SULength")
        if taxa_cu != taxa_su or len(taxa_cu) != 48 or set(r["split_key"] for r in rec_cu if r["branch_type"] == "internal") != set(r["split_key"] for r in rec_su if r["branch_type"] == "internal"):
            raise ValueError(f"taxa or topology mismatch: {d}")
        run = json.loads((d / "run.json").read_text()); prov = json.loads((d / "provenance.json").read_text())
        expected = run.get("output_sha256", {})
        if expected.get("tree.cu.nwk") != sha256(cu_tree) or expected.get("tree.su.nwk") != sha256(su_tree): raise ValueError(f"checksum mismatch: {d}")
        if prov.get("output_sha256") and prov["output_sha256"] != expected: raise ValueError(f"provenance mismatch: {d}")
        input_sidecar = root / Path(row["input_relative_path"]).with_suffix(".json")
        if not input_sidecar.exists(): raise ValueError(f"missing input provenance sidecar: {input_sidecar}")
        input_meta = json.loads(input_sidecar.read_text())
        if input_meta.get("input_sha256") != prov.get("input_sha256"): raise ValueError(f"input checksum mismatch: {d}")
        if str(input_meta.get("seed")) != str(row.get("seed")) or str(input_meta.get("epsilon_requested")) != str(row.get("epsilon_requested")): raise ValueError(f"input manifest metadata mismatch: {d}")
        if int(input_meta.get("total_records", -1)) != 62945: raise ValueError(f"record-count mismatch: {d}")
        if int(input_meta.get("n_chr4_records", -1)) != int(row["n_chr4_records"]): raise ValueError(f"chr4-count mismatch: {d}")
        if int(input_meta.get("n_background_records", -1)) != int(row["n_background_records"]): raise ValueError(f"background-count mismatch: {d}")
        if not run.get("runs") or any(int(x.get("returncode", 1)) != 0 for x in run["runs"]): raise ValueError(f"ASTRAL run record is not successful: {d}")
        command_text = json.dumps(run.get("command", []))
        if run.get("mode") != "topology_su" or '"--length", "SULength"' not in command_text: raise ValueError(f"not a single topology_su/SULength run: {d}")
        result.append({"mixture_id": row["mixture_id"], "epsilon": row["epsilon_requested"], "replicate": row["replicate"], "seed": row["seed"], "realized_epsilon": row["realized_epsilon"], "total_records": row["total_records"], "n_chr4_records": row["n_chr4_records"], "n_background_records": row["n_background_records"], "n_taxa": len(taxa_cu), "n_internal_splits": len([r for r in rec_cu if r["branch_type"] == "internal"]), "cu_su_topology_agree": "true", "exit_code": 0, "tree_cu_sha256": sha256(cu_tree), "tree_su_sha256": sha256(su_tree), "input_sha256": prov.get("input_sha256", ""), "astral_binary_sha256": prov.get("aster_bin_sha256", ""), "mode": run.get("mode", ""), "checksum_valid": "true", "metadata_valid": "true"})
    return result

def main():
    p = argparse.ArgumentParser(); p.add_argument("--root", type=Path, required=True); p.add_argument("--manifest", type=Path, required=True); p.add_argument("--only"); p.add_argument("--output", type=Path); args = p.parse_args()
    result = validate(args.root, args.manifest, args.only)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fields = list(result[0]) if result else ["status"]
        with args.output.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(result)
    print(json.dumps({"validated": len(result), "output": str(args.output) if args.output else None}, indent=2))
if __name__ == "__main__": main()
