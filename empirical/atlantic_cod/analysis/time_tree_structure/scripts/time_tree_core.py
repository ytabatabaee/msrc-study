"""Rooted ultrametric SNAPP time-tree feature and triplet calculations."""
from __future__ import annotations

import itertools
import math
from Bio import Phylo


def extract_mrca_ages(newick, expected_taxa=None, ultrametric_rtol=1e-6, ultrametric_atol=1e-10):
    """Return pairwise MRCA ages backward from present and strict validity audit.

    Root stems are excluded by using root-to-node path lengths (Bio.Phylo's
    distance from the root excludes the root clade's own branch length).
    Missing or negative non-root branches invalidate the tree; zero branches
    are retained and can induce unresolved triplets.
    """
    try:
        tree = Phylo.read(__import__('io').StringIO(newick), 'newick')
    except Exception as exc:
        return {'valid': False, 'issue': f'newick_parse_error:{type(exc).__name__}', 'taxa': (), 'root_age': math.nan, 'ultrametric_spread': math.nan, 'pair_ages': {}}
    root = tree.root
    taxa = tuple(sorted(t.name for t in tree.get_terminals() if t.name))
    issues = []
    if len(taxa) != len(tree.get_terminals()) or len(set(taxa)) != len(taxa):
        issues.append('missing_or_duplicate_taxon_labels')
    if expected_taxa is not None and set(taxa) != set(expected_taxa):
        issues.append('inconsistent_taxon_set')
    for clade in tree.find_clades(order='preorder'):
        if clade is root:
            continue
        x = clade.branch_length
        if x is None:
            issues.append('missing_branch_length')
        elif not math.isfinite(float(x)):
            issues.append('nonfinite_branch_length')
        elif float(x) < 0:
            issues.append('negative_branch_length')
    if issues:
        return {'valid': False, 'issue': ';'.join(sorted(set(issues))), 'taxa': taxa, 'root_age': math.nan, 'ultrametric_spread': math.nan, 'pair_ages': {}, 'tree': tree}
    tips = tree.get_terminals()
    if len(tips) < 3:
        return {'valid': False, 'issue': 'fewer_than_three_tips', 'taxa': taxa, 'root_age': math.nan, 'ultrametric_spread': math.nan, 'pair_ages': {}, 'tree': tree}
    tip_depth = {t.name: float(tree.distance(root, t)) for t in tips}
    root_age = max(tip_depth.values())
    spread = max(tip_depth.values()) - min(tip_depth.values())
    tol = max(ultrametric_atol, abs(root_age) * ultrametric_rtol)
    if not math.isfinite(root_age) or root_age <= 0:
        issues.append('zero_or_invalid_root_age')
    if spread > tol:
        issues.append('non_ultrametric')
    pair_ages = {}
    if not issues:
        by_name = {t.name: t for t in tips}
        for a, b in itertools.combinations(taxa, 2):
            ancestor = tree.common_ancestor(by_name[a], by_name[b])
            age = root_age - float(tree.distance(root, ancestor))
            if age < -tol:
                issues.append('negative_mrca_age')
                break
            pair_ages[(a, b)] = max(0.0, age) if age >= -tol else age
    return {'valid': not issues, 'issue': ';'.join(sorted(set(issues))), 'taxa': taxa, 'root_age': root_age, 'ultrametric_spread': spread, 'pair_ages': pair_ages if not issues else {}, 'tree': tree}


def triplet_features(taxa, pair_ages, root_age, unresolved_rtol=1e-6, unresolved_atol=1e-10):
    """Return rooted topology, separation, and all three signed Q contrasts."""
    out = []
    tol = max(unresolved_atol, abs(root_age) * unresolved_rtol)
    for a, b, c in itertools.combinations(sorted(taxa), 3):
        pairs = {tuple(sorted((a, b))): pair_ages[tuple(sorted((a, b)))],
                 tuple(sorted((a, c))): pair_ages[tuple(sorted((a, c)))],
                 tuple(sorted((b, c))): pair_ages[tuple(sorted((b, c)))]}
        # Closest pair has the youngest MRCA in an ultrametric rooted tree.
        order = sorted(pairs.items(), key=lambda x: (x[1], x[0]))
        sep = max(pairs.values()) - min(pairs.values())
        resolved = sep > tol and (len(order) < 2 or order[1][1] - order[0][1] > tol)
        topology = f'{order[0][0][0]}+{order[0][0][1]}|{next(x for x in (a,b,c) if x not in order[0][0])}' if resolved else 'unresolved'
        q = {}; separation_by_resolution = {}; resolution_topologies={}
        for index,(x, y, z) in enumerate(((a,b,c),(a,c,b),(b,c,a)),1):
            split=f'{x}+{y}|{z}';resolution=f'r{index}'
            resolution_topologies[resolution]=split
            q[resolution] = (pairs[tuple(sorted((x,z)))] + pairs[tuple(sorted((y,z)))]) / 2 - pairs[tuple(sorted((x,y)))]
            separation_by_resolution[resolution] = max(pairs.values()) - pairs[tuple(sorted((x,y)))]
        out.append({'triplet': f'{a},{b},{c}', 'taxa': (a,b,c), 'topology': topology,
                    'temporal_separation': sep if resolved else 0.0,
                    'normalized_temporal_separation': sep/root_age if resolved and root_age > 0 else 0.0,
                    'unresolved': not resolved,
                    'resolution_topologies': resolution_topologies,
                    'q_contrasts': q,
                    'separation_by_resolution': separation_by_resolution,
                    'normalized_separation_by_resolution': {k:v/root_age for k,v in separation_by_resolution.items()} if root_age > 0 else {k:math.nan for k in separation_by_resolution},
                    'normalized_q_contrasts': {k:v/root_age for k,v in q.items()} if root_age > 0 else {k:math.nan for k in q}})
    return out
