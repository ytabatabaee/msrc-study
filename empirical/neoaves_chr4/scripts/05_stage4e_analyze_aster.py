#!/usr/bin/env python3
"""Analyze collected Stage 4E ASTRAL-IV topology, CU, and SU trees.

No missing output is converted to zero.  CU and SU trees are loaded from their
explicit suffixes and branch matching is by normalized unrooted bipartition
for internal branches or exact taxon label for terminal branches.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from stage4d_core import ROOT, RESULTS, rf_comparison, focal_split_states, frozen_focal_groups, read_tsv

BASE = ROOT / "empirical/neoaves_chr4"
TREE_DIR = BASE / "results/stage4e/trees"
OUT = BASE / "results/stage4e"


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def manifest_counts() -> dict[str, dict[str, str]]:
    path = RESULTS / "stage4e_input_manifest.tsv"
    return {row["treatment"]: row for row in read_tsv(path)} if path.exists() else {}


def parse_label(label: object) -> dict[str, str]:
    text = "" if label is None else str(label)
    out = {key: "" for key in ("localPP", "q1", "q2", "q3", "CULength", "SULength")}
    for key in out:
        match = re.search(rf"{key}=([^,;:)\]]+)", text)
        if match:
            out[key] = match.group(1)
    if text and not out["localPP"]:
        try:
            float(text); out["localPP"] = text
        except ValueError:
            pass
    return out


def number(value: object) -> float | None:
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def tree_records(newick: str) -> tuple[set[str], dict[str, dict]]:
    import dendropy
    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-unrooted", preserve_underscores=True)
    taxa = {node.taxon.label for node in tree.leaf_node_iter()}
    if len(taxa) != len(list(tree.leaf_node_iter())):
        raise ValueError("duplicate taxon labels")
    records = {}
    for node in tree.postorder_node_iter():
        if node.parent_node is None:
            continue
        annotation_values = {str(annotation.name): str(annotation.value)
                             for annotation in node.annotations}
        parsed_label = parse_label(node.label)
        for key in ("localPP", "q1", "q2", "q3", "CULength", "SULength"):
            if not parsed_label[key] and key in annotation_values:
                parsed_label[key] = annotation_values[key]
        side = frozenset(n.taxon.label for n in node.leaf_iter())
        other = frozenset(taxa - side)
        branch_type = "terminal" if len(side) == 1 or len(other) == 1 else "internal"
        if branch_type == "internal":
            key = "|".join(sorted(min(tuple(sorted(side)), tuple(sorted(other)))))
        else:
            key = "|".join(sorted(side if len(side) == 1 else other))
        records[key] = {"branch_id": key, "branch_type": branch_type,
                        "length": node.edge.length, "label": node.label,
                        **parsed_label}
    return taxa, records


def load_trees(unit: str) -> dict[str, str]:
    suffix = "cu" if unit == "CU" else "su"
    trees = {}
    for path in sorted(TREE_DIR.glob(f"*.{suffix}.nwk")):
        treatment = path.name[: -len(f".{suffix}.nwk")]
        trees[treatment] = path.read_text()
    # Legacy single-output CU trees are accepted only when no explicit CU tree
    # exists for that treatment, preserving old Stage 4E compatibility.
    for path in sorted(TREE_DIR.glob("*.nwk")):
        if path.name.endswith((".cu.nwk", ".su.nwk")):
            continue
        trees.setdefault(path.stem, path.read_text()) if unit == "CU" else None
    return trees


def topology_rows(trees: dict[str, str]) -> list[dict]:
    counts = manifest_counts(); t0 = trees.get("T0"); non = trees.get("NONCHR4"); rows = []
    for treatment, newick in sorted(trees.items()):
        row = {"treatment": treatment, "n_loci": counts.get(treatment, {}).get("n_loci", ""), "n_chr4_loci": counts.get(treatment, {}).get("n_chr4_loci", "")}
        if t0:
            c = rf_comparison(newick, t0); row.update(rf_to_T0=c["rf"], nrf_to_T0=c["normalized_rf"], same_as_T0=c["rf"] == 0)
        else:
            row.update(rf_to_T0="PENDING_CLUSTER_OUTPUT", nrf_to_T0="PENDING_CLUSTER_OUTPUT", same_as_T0="")
        if non:
            c = rf_comparison(newick, non); row.update(rf_to_NONCHR4=c["rf"], nrf_to_NONCHR4=c["normalized_rf"], same_as_NONCHR4=c["rf"] == 0)
        else:
            row.update(rf_to_NONCHR4="PENDING_CLUSTER_OUTPUT", nrf_to_NONCHR4="PENDING_CLUSTER_OUTPUT", same_as_NONCHR4="")
        rows.append(row)
    return rows


def changed_splits(trees: dict[str, str]) -> list[dict]:
    if "T0" not in trees:
        return [{"treatment": "PENDING_CLUSTER_OUTPUT", "splits_lost_from_T0": "PENDING_CLUSTER_OUTPUT", "splits_gained_vs_T0": "PENDING_CLUSTER_OUTPUT", "splits_matching_NONCHR4": "PENDING_CLUSTER_OUTPUT"}]
    _, base = tree_records(trees["T0"]); base = {k for k,v in base.items() if v["branch_type"] == "internal"}; non = set()
    if "NONCHR4" in trees:
        _, nonrec = tree_records(trees["NONCHR4"]); non = {k for k,v in nonrec.items() if v["branch_type"] == "internal"}
    out=[]
    for treatment,newick in sorted(trees.items()):
        _, rec = tree_records(newick); keys={k for k,v in rec.items() if v["branch_type"] == "internal"}
        out.append({"treatment":treatment,"splits_lost_from_T0":len(base-keys),"splits_gained_vs_T0":len(keys-base),"splits_matching_NONCHR4":len(keys&non) if non else "PENDING_CLUSTER_OUTPUT"})
    return out


def focal_rows(trees: dict[str, str], unit: str) -> list[dict]:
    if not trees or "T0" not in trees:
        return [{"treatment":"PENDING_CLUSTER_OUTPUT","unit":unit,"clade":"PENDING_CLUSTER_OUTPUT","state":"PENDING_CLUSTER_OUTPUT"}]
    groups=frozen_focal_groups(); base={r["clade"]:r for r in focal_split_states(trees["T0"],groups)}; rows=[]
    for treatment,newick in sorted(trees.items()):
        for row in focal_split_states(newick,groups):
            parsed = parse_label(row.get("branch_support")); row.update(treatment=treatment,unit=unit,CULength=parsed["CULength"],SULength=parsed["SULength"],matches_T0=row["state"]==base.get(row["clade"],{}).get("state")); rows.append(row)
    return rows


def branch_tables(trees: dict[str, str], unit: str) -> tuple[list[dict], list[dict]]:
    if not trees:
        return ([{"unit":unit,"status":"PENDING_CLUSTER_OUTPUT"}], [{"unit":unit,"status":"PENDING_CLUSTER_OUTPUT"}])
    base_records = tree_records(trees["T0"])[1] if "T0" in trees else {}
    lengths=[]; changes=[]
    for treatment,newick in sorted(trees.items()):
        _, rec = tree_records(newick)
        for branch_id,data in sorted(rec.items()):
            annotation_key = "CULength" if unit == "CU" else "SULength"
            annotated = number(data.get(annotation_key))
            value = annotated if annotated is not None else number(data["length"])
            if value is not None and (value < 0 or not math.isfinite(value)):
                raise ValueError(f"invalid {unit} branch length: {treatment} {branch_id}")
            lengths.append({"unit":unit,"treatment":treatment,"branch_id":branch_id,"branch_type":data["branch_type"],"length":value if value is not None else "","label":data["label"] or ""})
            if treatment == "T0" or branch_id not in base_records:
                continue
            base_annotated = number(base_records[branch_id].get(annotation_key))
            base = base_annotated if base_annotated is not None else number(base_records[branch_id]["length"])
            delta = value-base if value is not None and base is not None else ""
            ratio = value/base if value is not None and base is not None and abs(base)>1e-12 else ""
            changes.append({"unit":unit,"treatment":treatment,"branch_id":branch_id,"branch_type":data["branch_type"],"baseline_length":base if base is not None else "","treatment_length":value if value is not None else "","delta":delta,"ratio":ratio,"percent_change":(100*(ratio-1) if ratio != "" else ""),"baseline_support":base_records[branch_id]["label"] or "","treatment_support":data["label"] or ""})
    return lengths, changes


def random_rows() -> list[dict]:
    path=RESULTS/"stage4e_random_matched_manifest.tsv"
    if not path.exists(): return [{"regime":"ALL","status":"PENDING_CLUSTER_OUTPUT"}]
    return [{"regime":r["regime"],"status":"PENDING_CLUSTER_OUTPUT","replicate":r["replicate"]} for r in read_tsv(path)]


def write_pending_docs():
    (OUT/"stage4e_cu_su_validation.md").write_text("""# Stage 4E CU/SU validation\n\nThe verified ASTER source interface uses `--length CULength` for coalescent-unit branch lengths and `--length SULength` for CASTLES substitution-per-site branch lengths. SU estimation uses `--genelength` (default 1000) and `--sulengthtype` (1: species-tree branch length; 2: mean gene-tree branch length). The local smoke test is generated from a four-taxon input with nonzero gene-tree branches and validates both output modes.\n\nNo 363-species ASTER output is present in the repository at collection time; all empirical treatment and random-control comparisons therefore remain pending. Missing CU/SU estimates are retained as missing and never converted to zero.\n""")
    (OUT/"stage4e_cu_su_methods.md").write_text("""# Methods\n\nPrimary treatments retain the frozen Stage 4E input manifests. The runner writes separate `tree.cu.nwk` and `tree.su.nwk` outputs for `--mode both`, records the exact commands and binary checksum, and requires SLURM for full-scale execution. Internal branches are matched by normalized unrooted bipartition; terminal branches by exact taxon label. Relative changes are omitted when the baseline is zero or near zero.\n""")
    (OUT/"stage4e_cu_su_results.md").write_text("# Results\n\nFull-scale Stage 4E CU/SU results are pending returned and validated Bridges-2 outputs. The collector reports each treatment and unit independently.\n")


def write_manifest(cu: dict[str, str], su: dict[str, str]) -> None:
    payload = {
        "stage": "neoaves_chr4_stage4f",
        "status": "complete" if cu and su else "pending_cluster_output",
        "units": {"CU": "coalescent units", "SU": "CASTLES substitution-per-site units"},
        "runner_modes": ["topology_cu", "topology_su", "both"],
        "fixed_topology": {"supported_command": "-C -c <species-tree>", "status": "pending_validated_T0"},
        "cu_treatments": sorted(cu),
        "su_treatments": sorted(su),
        "missing_empirical_outputs_are_not_zero": True,
    }
    (OUT / "stage4e_cu_su_manifest.json").write_text(json.dumps(payload, indent=2) + "\n")


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--allow-pending",action="store_true"); args=parser.parse_args()
    cu=load_trees("CU"); su=load_trees("SU")
    if not cu and not su and not args.allow_pending: raise SystemExit("PENDING_CLUSTER_OUTPUT: no completed Stage 4E trees are available")
    tops=topology_rows(cu); write_tsv(OUT/"stage4e_topology_comparison.tsv",tops,["treatment","n_loci","n_chr4_loci","rf_to_T0","nrf_to_T0","rf_to_NONCHR4","nrf_to_NONCHR4","same_as_T0","same_as_NONCHR4"])
    write_tsv(OUT/"stage4e_tree_comparison.tsv",tops,["treatment","n_loci","n_chr4_loci","rf_to_T0","nrf_to_T0","rf_to_NONCHR4","nrf_to_NONCHR4","same_as_T0","same_as_NONCHR4"])
    write_tsv(OUT/"stage4e_changed_splits.tsv",changed_splits(cu),["treatment","splits_lost_from_T0","splits_gained_vs_T0","splits_matching_NONCHR4"])
    write_tsv(OUT/"stage4e_focal_relationships.tsv",focal_rows(cu,"CU")+focal_rows(su,"SU"),["treatment","unit","clade","state","branch_support","CULength","SULength","matches_T0","n_C1","n_C2","n_S","n_O"])
    cu_lengths,cu_changes=branch_tables(cu,"CU"); su_lengths,su_changes=branch_tables(su,"SU")
    write_tsv(OUT/"stage4e_branch_lengths_CU.tsv",cu_lengths,["unit","treatment","branch_id","branch_type","length","label"])
    write_tsv(OUT/"stage4e_branch_lengths_SU.tsv",su_lengths,["unit","treatment","branch_id","branch_type","length","label"])
    write_tsv(OUT/"stage4e_branch_length_changes_CU.tsv",cu_changes,["unit","treatment","branch_id","branch_type","baseline_length","treatment_length","delta","ratio","percent_change","baseline_support","treatment_support"])
    write_tsv(OUT/"stage4e_branch_length_changes_SU.tsv",su_changes,["unit","treatment","branch_id","branch_type","baseline_length","treatment_length","delta","ratio","percent_change","baseline_support","treatment_support"])
    write_tsv(OUT/"stage4e_branch_support.tsv",[],["treatment","split","localPP","q1","q2","q3","CULength"])
    write_tsv(OUT/"stage4e_branch_length_sensitivity.tsv",cu_changes,["unit","treatment","branch_id","branch_type","baseline_length","treatment_length","delta","ratio","percent_change","baseline_support","treatment_support"])
    write_tsv(OUT/"stage4e_fixed_topology_comparison.tsv",[{"status":"PENDING_CLUSTER_OUTPUT","note":"Run ASTRAL-IV with -C -c only after validated T0 topology is returned."}],["status","note"])
    write_tsv(OUT/"stage4e_random_matched_comparison.tsv",random_rows(),["regime","status","replicate"])
    write_tsv(OUT/"stage4e_branch_length_outliers.tsv",[],["unit","treatment","branch_id","branch_type","length","reason"])
    write_tsv(OUT/"stage4e_treatment_summary.tsv",[{"treatment":t,"cu_status":"complete" if t in cu else "pending","su_status":"complete" if t in su else "pending"} for t in sorted(set(cu)|set(su)|set(manifest_counts()))],["treatment","cu_status","su_status"])
    write_pending_docs()
    write_manifest(cu, su)
    print(json.dumps({"cu_treatments":sorted(cu),"su_treatments":sorted(su),"status":"complete" if cu and su else "pending_cluster_output"},indent=2))


if __name__ == "__main__": main()
