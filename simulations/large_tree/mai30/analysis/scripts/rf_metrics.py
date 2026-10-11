"""Strict unrooted Robinson--Foulds metrics for Newick trees.

Branch lengths and rooting are ignored.  The parser nevertheless preserves
the rooted representation long enough to enumerate every edge bipartition.
For unresolved trees the normalization denominator is the sum of the numbers
of nontrivial splits in the two input trees.  For two binary n-taxon trees it
reduces to 2(n-3).
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass


class TreeValidationError(ValueError):
    pass


@dataclass(frozen=True)
class Tree:
    taxa: frozenset[str]
    splits: frozenset[frozenset[str]]
    leaf_count: int
    resolved: bool


@dataclass(frozen=True)
class RFResult:
    distance: int
    normalized: float
    denominator: int
    n_splits_a: int
    n_splits_b: int
    resolved_a: bool
    resolved_b: bool


_TOKEN = re.compile(r"[(),:;]|[^(),:;\s]+")


class _Node:
    def __init__(self, label: str | None = None, children: list["_Node"] | None = None):
        self.label = label
        self.children = children or []


def parse_newick(text: str) -> Tree:
    tokens = _TOKEN.findall(text.strip())
    if not tokens:
        raise TreeValidationError("empty Newick record")
    pos = 0

    def subtree() -> _Node:
        nonlocal pos
        if pos >= len(tokens):
            raise TreeValidationError("unexpected end of Newick record")
        if tokens[pos] == "(":
            pos += 1
            children = [subtree()]
            while pos < len(tokens) and tokens[pos] == ",":
                pos += 1
                children.append(subtree())
            if pos >= len(tokens) or tokens[pos] != ")":
                raise TreeValidationError("expected ')' in Newick record")
            pos += 1
            label = None
            if pos < len(tokens) and tokens[pos] not in ",():;":
                label = tokens[pos]
                pos += 1
            if pos < len(tokens) and tokens[pos] == ":":
                pos += 1
                if pos >= len(tokens) or tokens[pos] in ",():;":
                    raise TreeValidationError("missing internal branch length")
                try:
                    length = float(tokens[pos])
                except ValueError as exc:
                    raise TreeValidationError("invalid internal branch length") from exc
                if not math.isfinite(length) or length < 0:
                    raise TreeValidationError("negative or non-finite branch length")
                pos += 1
            if len(children) < 2:
                raise TreeValidationError("internal node has fewer than two children")
            return _Node(label, children)

        label = tokens[pos]
        if label in ",():;":
            raise TreeValidationError("missing leaf label")
        pos += 1
        if pos < len(tokens) and tokens[pos] == ":":
            pos += 1
            if pos >= len(tokens) or tokens[pos] in ",():;":
                raise TreeValidationError("missing leaf branch length")
            try:
                length = float(tokens[pos])
            except ValueError as exc:
                raise TreeValidationError("invalid leaf branch length") from exc
            if not math.isfinite(length) or length < 0:
                raise TreeValidationError("negative or non-finite branch length")
            pos += 1
        return _Node(label=label)

    root = subtree()
    if pos < len(tokens) and tokens[pos] == ";":
        pos += 1
    if pos != len(tokens):
        raise TreeValidationError(f"unparsed Newick tokens: {tokens[pos:]}")

    leaves: list[str] = []
    splits: set[frozenset[str]] = set()

    def canonical(side: frozenset[str], all_taxa: frozenset[str]) -> frozenset[str] | None:
        other = all_taxa - side
        if len(side) < 2 or len(other) < 2:
            return None
        a = tuple(sorted(side))
        b = tuple(sorted(other))
        return side if a <= b else frozenset(other)

    def visit(node: _Node, all_taxa: frozenset[str] | None = None) -> frozenset[str]:
        if not node.children:
            leaves.append(node.label or "")
            return frozenset([node.label or ""])
        taxa = frozenset().union(*(visit(child, all_taxa) for child in node.children))
        if all_taxa is not None:
            for child in node.children:
                child_taxa = leaf_taxa(child)
                split = canonical(child_taxa, all_taxa)
                if split is not None:
                    splits.add(split)
        return taxa

    def leaf_taxa(node: _Node) -> frozenset[str]:
        if not node.children:
            return frozenset([node.label or ""])
        return frozenset().union(*(leaf_taxa(child) for child in node.children))

    all_taxa = leaf_taxa(root)
    visit(root, all_taxa)
    if "" in leaves:
        raise TreeValidationError("empty taxon label")
    if len(leaves) != len(set(leaves)):
        raise TreeValidationError("duplicate taxon label")
    if len(leaves) < 4:
        raise TreeValidationError("RF comparison requires at least four taxa")
    expected_binary_splits = len(leaves) - 3
    resolved = len(splits) == expected_binary_splits
    return Tree(frozenset(leaves), frozenset(splits), len(leaves), resolved)


def rf_distance(a: Tree, b: Tree) -> RFResult:
    if a.taxa != b.taxa:
        missing_a = sorted(b.taxa - a.taxa)
        missing_b = sorted(a.taxa - b.taxa)
        raise TreeValidationError(f"taxon sets differ; missing_from_a={missing_a}; missing_from_b={missing_b}")
    distance = len(a.splits ^ b.splits)
    denominator = len(a.splits) + len(b.splits)
    normalized = 0.0 if denominator == 0 else distance / denominator
    if not 0.0 <= normalized <= 1.0:
        raise AssertionError(f"normalized RF outside [0,1]: {normalized}")
    return RFResult(distance, normalized, denominator, len(a.splits), len(b.splits), a.resolved, b.resolved)


def read_tree_file(path: str, expected_count: int | None = None) -> list[Tree]:
    from pathlib import Path
    trees: list[Tree] = []
    with Path(path).open() as handle:
        for line_no, line in enumerate(handle, 1):
            if line.strip():
                try:
                    trees.append(parse_newick(line))
                except TreeValidationError as exc:
                    raise TreeValidationError(f"{path}:{line_no}: {exc}") from exc
    if not trees:
        raise TreeValidationError(f"{path}: no Newick records")
    taxa = trees[0].taxa
    for index, tree in enumerate(trees, 1):
        if tree.taxa != taxa:
            raise TreeValidationError(f"{path}: tree {index} has a different taxon set")
    if expected_count is not None and len(trees) != expected_count:
        raise TreeValidationError(f"{path}: expected {expected_count} trees, found {len(trees)}")
    return trees
