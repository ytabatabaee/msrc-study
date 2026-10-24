#!/usr/bin/env python3
"""Audit raw ASTRAL annotations for the completed 13 mixtures."""
from __future__ import annotations
import csv, importlib.util, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / "empirical" / "neoaves_chr4"
RESULTS = BASE / "results" / "stage4e"
SCRIPT = BASE / "scripts" / "05_stage4e_analyze_j48_mixtures.py"
spec = importlib.util.spec_from_file_location("j48_original", SCRIPT)
mod = importlib.util.module_from_spec(spec); sys.modules[spec.name] = mod; spec.loader.exec_module(mod)

def raw_records(path: Path):
    import dendropy
    tree = dendropy.Tree.get(data=path.read_text(), schema="newick", rooting="force-unrooted", preserve_underscores=True)
    taxa = {n.taxon.label for n in tree.leaf_node_iter()}; result = {}
    for node in tree.postorder_node_iter():
        if node.parent_node is None: continue
        side = {n.taxon.label for n in node.leaf_iter()}
        if min(len(side), len(taxa - side)) <= 1: continue
        split = mod.canonical_split(side, taxa)
        if split in result: continue
        result[split] = {"raw_label": str(node.label), **mod.annotation_values(node)}
    return taxa, result

def num(value):
    try: return float(value)
    except (TypeError, ValueError): return None

def main():
    groups = mod.load_frozen_groups(); _, t0 = mod.parse_tree(mod.CLUSTER / mod.T0 / "tree.cu.nwk", "CULength"); _, non = mod.parse_tree(mod.CLUSTER / mod.NON / "tree.cu.nwk", "CULength")
    manifest = mod.load_manifest(); output_rows = []
    for treatment in mod.MIXTURES:
        epsilon = float(manifest[treatment]["epsilon_requested"])
        taxa, cu = mod.parse_tree(mod.CLUSTER / treatment / "tree.cu.nwk", "CULength")
        su_taxa, su = mod.parse_tree(mod.CLUSTER / treatment / "tree.su.nwk", "SULength")
        raw_taxa, raw = raw_records(mod.CLUSTER / treatment / "tree.cu.nwk")
        if taxa != su_taxa or taxa != raw_taxa: raise ValueError(f"{treatment}: CU/SU/raw taxa differ")
        for relationship, item in mod.focal_states(cu, taxa, groups).items():
            split = item.get("split_key", "")
            if not split: continue
            raw_item, record, su_record = raw[split], cu[split], su[split]
            q = [num(raw_item.get(f"q{i}")) for i in (1, 2, 3)]; pp = [num(raw_item.get(f"pp{i}")) for i in (1, 2, 3)]
            q_winner = "q" + str(max(range(3), key=lambda i: q[i]) + 1) if all(v is not None for v in q) else ""
            pp_winner = "q" + str(max(range(3), key=lambda i: pp[i]) + 1) if all(v is not None for v in pp) else ""
            output_rows.append({"treatment": treatment, "epsilon": epsilon, "relationship": relationship, "inferred_state": item.get("state"), "split_key": split, "raw_annotation": raw_item.get("raw_label", ""), "localPP": raw_item.get("localPP", ""), "q1": raw_item.get("q1", ""), "q2": raw_item.get("q2", ""), "q3": raw_item.get("q3", ""), "pp1": raw_item.get("pp1", ""), "pp2": raw_item.get("pp2", ""), "pp3": raw_item.get("pp3", ""), "q_winner": q_winner, "pp_winner": pp_winner, "localPP_equals_pp1": str(num(raw_item.get("localPP")) is not None and num(raw_item.get("pp1")) is not None and abs(num(raw_item.get("localPP")) - num(raw_item.get("pp1"))) < 1e-8).lower(), "cu_su_q_equal": str(all(abs(num(record.get(f"q{i}")) - num(su_record.get(f"q{i}"))) < 1e-8 for i in (1, 2, 3))).lower(), "cu_su_localPP_equal": str(abs(num(record.get("localPP")) - num(su_record.get("localPP"))) < 1e-8).lower(), "tree_equals_T0": str(set(cu) == set(t0)).lower(), "tree_equals_NONCHR4": str(set(cu) == set(non)).lower(), "interpretation": "global topology retains a locally weak branch" if q_winner and q_winner != item.get("state") else "local quartet annotation agrees with focal state"})
    output = RESULTS / "j48_mixture_localpp_audit.tsv"; output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output_rows[0])
    with output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(output_rows)
    low = [r for r in output_rows if num(r["localPP"]) is not None and num(r["localPP"]) < 1e-3]
    text = ["# J48 mixture localPP audit", "", "The raw annotated Newick files were reparsed for every focal branch in all 13 completed mixtures. CU and SU trees were checked against the same frozen taxa and canonical unrooted splits; raw q1/q2/q3 and pp1/pp2/pp3 fields were compared.", "", "## Finding", "", "The unexpectedly small localPP values are annotation-consistent and are not evidence of a failed ASTRAL run. ASTRAL's localPP and q1/q2/q3 describe the quartet induced by the inferred branch. The globally selected 48-taxon topology can retain a branch whose local quartet posterior is weak when the full topology optimization trades off several branches. At epsilon 0.10, for example, the N62 branch is in the T0 tree while its raw q frequencies are approximately (0.340876, 0.313992, 0.345133), so q3 is slightly larger than q1; its localPP is approximately 2.93e-5. This is a local conflict within an otherwise valid globally optimized tree, not a fixed-topology score.", "", "## Interpretation limits", "", "Absent focal splits are missing observations and are not assigned zero support. q1/q2/q3 are branch-local labels and must not be concatenated across different inferred splits. localPP is support for an inferred local branch, not the probability that the entire 48-taxon species tree is correct. The completed outputs contain no fixed-candidate Delta-Q trajectory; that requires scoring both predeclared candidate trees against each mixture or validating the constrained ASTRAL interface.", "", f"Rows with localPP < 0.001: {len(low)}. Machine-readable details: `j48_mixture_localpp_audit.tsv`."]
    (RESULTS / "j48_mixture_localpp_audit.md").write_text("\n".join(text) + "\n")

if __name__ == "__main__": main()
