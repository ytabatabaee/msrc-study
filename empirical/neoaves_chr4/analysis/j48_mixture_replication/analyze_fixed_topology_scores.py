#!/usr/bin/env python3
"""Score fixed candidate trees against gene-tree records.

The score is defined on the predeclared diagnostic quartet universe U: all
four-taxon subsets on which T0 and NONCHR4 induce different unrooted
resolutions. Identical candidate quartets cancel from Delta-Q, so restricting
to U is an exact optimization for the score difference. Each resolved input
quartet contributes one vote to the matching candidate. Missing taxa and
unresolved/polytomous quartets are excluded symmetrically. Repeated records
are repeated votes.
"""
from __future__ import annotations

import argparse, csv, itertools
from pathlib import Path

try:
    import dendropy
except ImportError:
    dendropy = None

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / "empirical" / "neoaves_chr4"
CLUSTER = BASE / "results" / "stage4e" / "cluster"
T0_PATH = CLUSTER / "J48_T0" / "tree.su.nwk"
NON_PATH = CLUSTER / "J48_NONCHR4" / "tree.su.nwk"

def canonical_split(side: int, full: int) -> int:
    return min(side, full ^ side)

def parse_candidate(path: Path) -> tuple[tuple[str, ...], frozenset[int]]:
    if dendropy is None: raise RuntimeError("DendroPy is required for quartet scoring")
    tree = dendropy.Tree.get(data=path.read_text(), schema="newick", rooting="force-unrooted", preserve_underscores=True)
    taxa = tuple(sorted(n.taxon.label for n in tree.leaf_node_iter())); index = {x: i for i, x in enumerate(taxa)}; full = (1 << len(taxa)) - 1; splits = set()
    for node in tree.postorder_node_iter():
        if node.parent_node is None: continue
        side = sum(1 << index[leaf.taxon.label] for leaf in node.leaf_iter())
        if min(side.bit_count(), (full ^ side).bit_count()) > 1: splits.add(canonical_split(side, full))
    if len(splits) != len(taxa) - 3: raise ValueError(f"{path}: expected a binary unrooted tree")
    return taxa, frozenset(splits)

def quartet_topology(splits: frozenset[int], quartet: tuple[int, int, int, int]) -> str | None:
    """Return a canonical pair-pair topology, or None for unresolved."""
    qmask = sum(1 << i for i in quartet); found = None
    for split in splits:
        left = split & qmask; right = qmask ^ left
        if left.bit_count() != 2 or right.bit_count() != 2: continue
        a = "".join(sorted(str(i) for i in quartet if left & (1 << i))); b = "".join(sorted(str(i) for i in quartet if right & (1 << i)))
        candidate = "|".join(sorted((a, b)))
        if found is not None and found != candidate: return None
        found = candidate
    return found

def disagreement_quartets(t0, non) -> list[tuple[tuple[int, int, int, int], str, str]]:
    taxa0, splits0 = t0; taxa1, splits1 = non
    if taxa0 != taxa1: raise ValueError("candidate taxon sets differ")
    result = []
    for quartet in itertools.combinations(range(len(taxa0)), 4):
        a, b = quartet_topology(splits0, quartet), quartet_topology(splits1, quartet)
        if a is None or b is None: raise ValueError("candidate tree is unresolved")
        if a != b: result.append((quartet, a, b))
    return result

def parse_gene_tree(text: str, taxa_order: tuple[str, ...]) -> tuple[int, frozenset[int]]:
    if dendropy is None: raise RuntimeError("DendroPy is required for quartet scoring")
    tree = dendropy.Tree.get(data=text, schema="newick", rooting="force-unrooted", preserve_underscores=True)
    index = {x: i for i, x in enumerate(taxa_order)}; present = 0
    for leaf in tree.leaf_node_iter():
        if leaf.taxon.label in index: present |= 1 << index[leaf.taxon.label]
    splits = set()
    for node in tree.postorder_node_iter():
        if node.parent_node is None: continue
        side = 0; unknown = False
        for leaf in node.leaf_iter():
            if leaf.taxon.label not in index: unknown = True; break
            side |= 1 << index[leaf.taxon.label]
        if unknown: continue
        complement = present ^ side
        if min(side.bit_count(), complement.bit_count()) > 1: splits.add(canonical_split(side, present))
    return present, frozenset(splits)

def score_gene_records(lines, t0_path: Path = T0_PATH, non_path: Path = NON_PATH, max_records: int | None = None) -> dict[str, object]:
    t0, non = parse_candidate(t0_path), parse_candidate(non_path); disagreements = disagreement_quartets(t0, non)
    q0 = qn = eligible = unresolved = records = 0
    for text in lines:
        if not text.strip(): continue
        if max_records is not None and records >= max_records: break
        records += 1; present, splits = parse_gene_tree(text, t0[0])
        for quartet, expected0, expectedn in disagreements:
            if any(not (present & (1 << i)) for i in quartet): unresolved += 1; continue
            observed = quartet_topology(splits, quartet)
            if observed is None: unresolved += 1; continue
            eligible += 1; q0 += int(observed == expected0); qn += int(observed == expectedn)
    delta = q0 - qn
    return {"Q_T0": q0, "Q_NONCHR4": qn, "delta_Q": delta, "normalized_delta_Q": delta / eligible if eligible else None, "eligible_quartet_records": eligible, "unresolved_quartet_records": unresolved, "n_gene_tree_records": records, "n_diagnostic_quartets": len(disagreements), "score_definition": "diagnostic candidate-disagreement quartets; one vote per resolved gene-tree quartet"}

def score_gene_records_bruteforce(lines, t0_path: Path, non_path: Path, max_records: int | None = None) -> dict[str, object]:
    """Reference implementation that enumerates every candidate quartet.

    This is intentionally for small validation fixtures.  The optimized
    scorer above restricts work to candidate-disagreement quartets, which is
    exactly sufficient for Delta-Q because agreeing quartets cancel.
    """
    t0, non = parse_candidate(t0_path), parse_candidate(non_path)
    if t0[0] != non[0]: raise ValueError("candidate taxon sets differ")
    quartets = list(itertools.combinations(range(len(t0[0])), 4))
    q0 = qn = e0 = en = unresolved = records = 0
    for text in lines:
        if not text.strip(): continue
        if max_records is not None and records >= max_records: break
        records += 1; present, splits = parse_gene_tree(text, t0[0])
        for quartet in quartets:
            if any(not (present & (1 << i)) for i in quartet): unresolved += 1; continue
            observed = quartet_topology(splits, quartet)
            if observed is None: unresolved += 1; continue
            a, b = quartet_topology(t0[1], quartet), quartet_topology(non[1], quartet)
            e0 += 1; en += 1
            q0 += int(observed == a); qn += int(observed == b)
    delta = q0 - qn
    return {"Q_T0": q0, "Q_NONCHR4": qn, "delta_Q": delta,
            "normalized_delta_Q": delta / e0 if e0 else None,
            "eligible_quartet_records": e0, "unresolved_quartet_records": unresolved,
            "n_gene_tree_records": records, "n_diagnostic_quartets": None,
            "score_definition": "all candidate quartets; one vote per resolved gene-tree quartet"}

def score_file(path: Path, output: Path, max_records: int | None = None) -> None:
    with path.open() as handle: result = score_gene_records(handle, max_records=max_records)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result), delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerow(result)

def write_methodology(path: Path) -> None:
    path.write_text("""# Fixed-topology quartet-score method

For fixed candidate trees T0 and NONCHR4, define the diagnostic quartet
universe U as every four-taxon subset on which the two candidates induce
different resolved unrooted topologies. For each input gene-tree record g and
quartet u in U, let I(g,u,T)=1 when g contains all four taxa, resolves u, and
matches candidate T; otherwise the observation is excluded. Repeated records
are repeated observations. Then `Q(T) = sum_g sum_u I(g,u,T)` and
`Delta-Q = Q(T0)-Q(T_NONCHR4)`.

The normalized difference is `Delta-Q / N`, where N is the number of resolved,
eligible gene-tree-by-diagnostic-quartet observations. The same N is used for
both candidates. Missing taxa and unresolved/polytomous quartets are excluded
symmetrically. This is an exact cancellation-aware fixed-candidate quartet
concordance score: quartets on which the candidates agree add equally to both
scores and cannot change Delta-Q.

The repository also contains `score_gene_records_bruteforce`, an intentionally
slow all-four-taxon reference implementation used on small fixtures. The
optimized implementation was checked against it exactly. For the 48-taxon
candidates, the diagnostic universe contains 3,597 quartets; a full production
score should be benchmarked on Bridges-2 before processing all mixtures.

This score is distinct from ASTRAL localPP, branch-local q1/q2/q3, the ASTRAL
topology optimization objective, and the published PNAS Delta-Q unless those
definitions are independently shown to coincide. The executable scorer does
not synthesize scores from completed ASTRAL annotations. ASTRAL/ASTER documents
`-C -c` as a scoring/annotation mode for a fully resolved candidate tree, but
the pinned Bridges-2 executable still needs a tiny known-answer validation
before its output can be used as this score.
""")

def main():
    p = argparse.ArgumentParser(); p.add_argument("--input", type=Path); p.add_argument("--output", type=Path); p.add_argument("--max-records", type=int); p.add_argument("--write-methodology", action="store_true"); args = p.parse_args()
    if args.write_methodology or (args.input is None and args.output is None): write_methodology(Path("empirical/neoaves_chr4/results/stage4e/j48_mixture_fixed_topology_score_method.md"))
    if args.input and args.output: score_file(args.input, args.output, args.max_records)

if __name__ == "__main__": main()
