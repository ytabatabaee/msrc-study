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
    result = []
    for row in rows(manifest):
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
        if input_meta.get("seed") != row.get("seed") or input_meta.get("epsilon_requested") != row.get("epsilon_requested"): raise ValueError(f"input manifest metadata mismatch: {d}")
        command_text = json.dumps(run.get("command", []))
        if run.get("mode") != "topology_su" or '"--length", "SULength"' not in command_text: raise ValueError(f"not a single topology_su/SULength run: {d}")
        result.append({"mixture_id": row["mixture_id"], "epsilon": row["epsilon_requested"], "replicate": row["replicate"], "n_taxa": len(taxa_cu), "n_internal_splits": len([r for r in rec_cu if r["branch_type"] == "internal"]), "cu_su_topology_agree": "true", "exit_code": 0, "tree_cu_sha256": sha256(cu_tree), "tree_su_sha256": sha256(su_tree), "input_sha256": prov.get("input_sha256", ""), "astral_binary_sha256": prov.get("aster_bin_sha256", ""), "mode": run.get("mode", "")})
    return result

def main():
    p = argparse.ArgumentParser(); p.add_argument("--root", type=Path, required=True); p.add_argument("--manifest", type=Path, required=True); p.add_argument("--only"); args = p.parse_args()
    print(json.dumps(validate(args.root, args.manifest, args.only), indent=2))
if __name__ == "__main__": main()
