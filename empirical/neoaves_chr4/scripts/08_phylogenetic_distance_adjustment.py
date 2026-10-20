#!/usr/bin/env python3
"""Neoaves chromosome 4 Stage 8: pair-specific phylogenetic adjustment.

This stage preserves the Stage 7 exploratory raw distance contrasts and adds
pair-specific residual contrasts based on outside-region ML patristic
distances.  Structural states are independent, multistate synteny labels;
they are not inferred from the local trees.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical/neoaves_chr4"
DATA = ROOT / "data/neoaves_chr4/processed"
RES = EMP / "results/stage8_phylogenetic_distance_adjustment"
FIG = EMP / "figures"
TREES = EMP / "external/stiller2024/63430.named.gene.trees.gz"
MANIFEST = EMP / "results/stage4d_locus_manifest.tsv"
LOCUS = DATA / "locus_table_chr4.tsv"
STATES = DATA / "canonical_arrangement_states.tsv"
REGIONS = DATA / "canonical_structural_regions.tsv"
SUPPORT = EMP / "results/stage6_genomic_tracks/window_quartet_support.tsv"
STAGE7 = EMP / "results/stage7_arrangement_divergence"
ROLES = DATA / "stage3b_v2_quartet_role_mapping.tsv"
J48 = EMP / "data/stage4d/jarvis48_taxa.tsv"

SPECIES = {
    "Stork": "Ciconia_maguari",
    "Flamingo": "Phoenicopterus_ruber",
    "Dove": "Columba_livia",
    "Sandgrouse": "Pterocles_gutturalis",
    "Turaco": "Tauraco_erythrolophus",
    "Cuckoo": "Cuculus_canorus",
}
TREE_ORDER = {tip: i for i, tip in enumerate(SPECIES.values())}
CLADES = ["Columbea", "N61", "N62"]
FLOAT_FIELDS = {
    "start", "end", "midpoint", "baseline_median", "baseline_q25",
    "baseline_q75", "distance", "residual", "relative_residual",
    "mean_same_residual", "mean_different_residual", "A_adjusted",
    "A_balanced", "mean_same_relative", "mean_different_relative",
    "A_relative_residual", "inside_mean", "outside_mean", "delta_inside_minus_outside",
}


def read(path: Path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})


def sha(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def finite(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def fmean(vals):
    vals = [float(x) for x in vals if finite(x)]
    return mean(vals) if vals else math.nan


def fmt(x):
    return "" if x is None or not finite(x) else f"{float(x):.12g}"


def parse_newick_distances(text, wanted):
    parent, edge, label, children = [], [], [], []
    i = 0

    def new(par=-1):
        k = len(parent)
        parent.append(par); edge.append(0.0); label.append(""); children.append([])
        return k

    def name():
        nonlocal i
        j = i
        while i < len(text) and text[i] not in ",():;":
            i += 1
        return text[j:i].strip()

    def length(k):
        nonlocal i
        if i < len(text) and text[i] == ":":
            i += 1; j = i
            while i < len(text) and text[i] not in ",();":
                i += 1
            try:
                edge[k] = float(text[j:i])
            except ValueError:
                edge[k] = math.nan

    def parse(par=-1):
        nonlocal i
        k = new(par)
        if i < len(text) and text[i] == "(":
            i += 1
            c = parse(k); children[k].append(c)
            while i < len(text) and text[i] == ",":
                i += 1; c = parse(k); children[k].append(c)
            if i >= len(text) or text[i] != ")":
                raise ValueError("unbalanced Newick")
            i += 1
            if i < len(text) and text[i] not in ":,();":
                label[k] = name()
        else:
            label[k] = name()
        length(k)
        return k

    parse()
    depth = [0.0] * len(parent)
    for k in range(1, len(parent)):
        depth[k] = depth[parent[k]] + edge[k] if math.isfinite(edge[k]) else math.nan
    leaves = {label[k]: k for k in range(len(parent)) if not children[k] and label[k]}
    present = [x for x in wanted if x in leaves]
    paths = {}
    for x in present:
        k = leaves[x]; d = {}
        while k >= 0:
            d[k] = depth[k]; k = parent[k]
        paths[x] = d
    out = {}
    for ai, a in enumerate(present):
        for b in present[ai + 1:]:
            common = set(paths[a]) & set(paths[b])
            if not common:
                continue
            lca = max(common, key=lambda k: depth[k])
            value = depth[leaves[a]] + depth[leaves[b]] - 2 * depth[lca]
            if math.isfinite(value):
                out[(a, b)] = value
    return out, len(leaves)


def pair_key(a, b):
    # The Newick parser stores pairs in the frozen SPECIES order.  Normalize
    # lookup order to that same order rather than alphabetical order.
    return (a, b) if TREE_ORDER.get(a, 10**9) <= TREE_ORDER.get(b, 10**9) else (b, a)


def state_pairs(state_by_species, species_names=None):
    names = species_names or list(SPECIES)
    names = [x for x in names if x in state_by_species]
    all_pairs = [(a, b) for ia, a in enumerate(names) for b in names[ia + 1:]]
    same = [p for p in all_pairs if state_by_species[p[0]]["normalized_arrangement_state"] == state_by_species[p[1]]["normalized_arrangement_state"]]
    diff = [p for p in all_pairs if state_by_species[p[0]]["normalized_arrangement_state"] != state_by_species[p[1]]["normalized_arrangement_state"]]
    return all_pairs, same, diff


def balanced_mean(values_by_category):
    vals = [fmean(v) for v in values_by_category.values() if v]
    return fmean(vals)


def main():
    RES.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    locus = read(LOCUS)
    unique = {r["Gene"]: r for r in locus}
    manifest_rows = [r for r in read(MANIFEST) if r["chromosome"] == "chr4"]
    expected_to_gene = {
        f"chr4_{r['start']}_{r['end']}_1k_start{r['locus_start_within_window']}": r["original_gene_id"]
        for r in manifest_rows
    }
    regions = read(REGIONS)
    states = read(STATES)
    state_by_region = defaultdict(dict)
    for r in states:
        state_by_region[r["canonical_region_id"]][r["species"]] = r
    support = read(SUPPORT)

    mapping_rows = []
    local_labels = set()
    # The local tip universe is collected while parsing; exact equality is required.
    for common, tip in SPECIES.items():
        mapping_rows.append({
            "structural_species": common, "local_tree_tip": tip,
            "exact_local_tip_name": tip, "mapping_status": "exact_frozen_stage7_mapping",
            "independent_state_source": "canonical_arrangement_states.tsv",
            "state_inferred_from_local_tree": "False",
            "note": "Ciconia is outside the exact 47-taxon Jarvis mapping but present in local trees" if common == "Stork" else "",
        })
    write(RES / "expanded_structural_state_mapping.tsv", mapping_rows, list(mapping_rows[0]))

    # Audit all structural tables for possible expansion. Only exact tip names and
    # independently supplied states are admissible; no phylogenetic imputation occurs.
    all_state_species = sorted({r["species"] for r in states})
    extra_rows = []
    for sp in all_state_species:
        extra_rows.append({"structural_species": sp, "local_tree_exact_tip": "", "independent_state_table": "canonical_arrangement_states.tsv", "expansion_status": "not_in_stage7_frozen_six_species_scope", "reason": "No exact local-tree tip/state validation was available for an additional species; not imputed"})
    (RES / "structural_taxon_expansion_audit.md").write_text(
        "# Structural taxon expansion audit\n\n"
        f"The frozen canonical state table contains {len(all_state_species)} species labels across {len(regions)} regions, but the Stage 7 exact compatibility audit supports six exact local-tree mappings: "
        + ", ".join(f"{k}={v}" for k, v in SPECIES.items()) + ".\n\n"
        "No additional species were added. The local trees contain many tips, but assigning a state to an additional tip would require an independently supported structural classification and an exact label mapping. Relatedness, local topology, or a partial name match was not used as an assignment rule. Ciconia is retained because its exact local tip is present even though it is absent from the 47-taxon Jarvis mapping.\n\n"
        "The canonical labels are region-specific multistate synteny/anchor states, with missing and uncertain assignments. They are independent of the local ML trees used for distances and QQS. The six-taxon limitation means all adjusted contrasts remain exploratory and cannot be treated as independent population-level replication.\n"
    )
    write(RES / "structural_taxon_expansion_audit.tsv", extra_rows, list(extra_rows[0]) if extra_rows else ["structural_species", "local_tree_exact_tip", "independent_state_table", "expansion_status", "reason"])

    # Eligibility is intentionally recomputed from Stage 7's rules.
    eligibility = []
    for reg in regions:
        rid = reg["canonical_region_id"]; st = state_by_region[rid]
        present = [x for x in SPECIES if x in st]
        all_pairs, same, diff = state_pairs(st, present)
        groups = Counter(st[x]["normalized_arrangement_state"] for x in present)
        confident = sum(st[x].get("confidence") == "confident" for x in present)
        eligible = len(present) >= 4 and bool(same) and len(diff) >= 2 and len(groups) >= 2
        eligibility.append({
            "canonical_region_id": rid, "start": reg["canonical_start"], "end": reg["canonical_end"],
            "n_state_taxa": len(present), "state_counts": ";".join(f"{k}:{v}" for k, v in sorted(groups.items())),
            "n_same_pairs": len(same), "n_different_state_pairs": len(diff), "n_state_groups": len(groups),
            "n_confident_states": confident, "n_uncertain_states": len(present) - confident,
            "eligible": str(eligible), "analysis_status": "exploratory_eligible" if eligible else "unsupported_insufficient_same_state_replication",
            "same_pairs": ";".join(f"{a}-{b}" for a, b in same), "different_pairs": ";".join(f"{a}-{b}" for a, b in diff),
        })
    write(RES / "region_eligibility_updated.tsv", eligibility, list(eligibility[0]))
    eligible_ids = [r["canonical_region_id"] for r in eligibility if r["eligible"] == "True"]

    # Parse each source tree once and retain only the six exact structural taxa.
    wanted = list(SPECIES.values())
    tree_distances = {}
    tree_taxa = {}
    parse_errors = 0
    with gzip.open(TREES, "rt") as h:
        for line in h:
            if "\t" not in line or not line.startswith("chr4_"):
                continue
            loc, nwk = line.rstrip("\n").split("\t", 1)
            gene = expected_to_gene.get(loc.removesuffix("_fasta.treefile"))
            if gene is None or gene in tree_distances:
                continue
            try:
                d, n = parse_newick_distances(nwk, wanted)
            except (ValueError, IndexError):
                parse_errors += 1; continue
            tree_distances[gene] = d; tree_taxa[gene] = n
    ordered_genes = sorted((g for g in tree_distances if g in unique), key=lambda g: float(unique[g]["midpoint"]))

    # Source-region masks. Strict excludes every canonical region; broad excludes only focal region.
    region_bounds = {r["canonical_region_id"]: (float(r["canonical_start"]), float(r["canonical_end"])) for r in regions}
    def in_interval(mid, interval):
        return interval[0] <= mid <= interval[1]
    def in_any_region(mid):
        return any(in_interval(mid, x) for x in region_bounds.values())

    bg_rows = []
    baselines = {}
    pair_meta = {}
    residual_rows = []
    adjusted_rows = []
    relative_rows = []
    all_struct_species = list(SPECIES)
    for rid in eligible_ids:
        st = state_by_region[rid]
        all_pairs, same, diff = state_pairs(st)
        pair_meta[rid] = {"all": all_pairs, "same": same, "diff": diff, "state": st}
        focal = region_bounds[rid]
        for baseline_type in ("broad_outside_focal", "strict_outside_all_structural_regions"):
            baseline = {}
            for p in all_pairs:
                vals = []
                for g in ordered_genes:
                    mid = float(unique[g]["midpoint"])
                    if baseline_type == "broad_outside_focal" and in_interval(mid, focal):
                        continue
                    if baseline_type == "strict_outside_all_structural_regions" and in_any_region(mid):
                        continue
                    d = tree_distances[g].get(pair_key(SPECIES[p[0]], SPECIES[p[1]]))
                    if finite(d): vals.append(float(d))
                baseline[p] = median(vals) if vals else math.nan
                bg_rows.append({
                    "canonical_region_id": rid, "baseline_type": baseline_type,
                    "species1": p[0], "species2": p[1],
                    "state1": st[p[0]]["normalized_arrangement_state"], "state2": st[p[1]]["normalized_arrangement_state"],
                    "state_relation": "same" if p in same else "different",
                    "baseline_median": fmt(baseline[p]), "baseline_q25": fmt(np.quantile(vals, .25) if vals else math.nan),
                    "baseline_q75": fmt(np.quantile(vals, .75) if vals else math.nan), "n_baseline_windows": len(vals),
                    "baseline_window_rule": baseline_type,
                })
            baselines[(rid, baseline_type)] = baseline

            # Pair residuals and the balanced/raw window contrasts.
            for gene in ordered_genes:
                r = unique[gene]; mid = float(r["midpoint"]); distances = tree_distances[gene]
                pair_values = {p: distances.get(pair_key(SPECIES[p[0]], SPECIES[p[1]])) for p in all_pairs}
                categories = {"same": [], "different": []}
                rel_categories = {"same": [], "different": []}
                category_groups = {"same": defaultdict(list), "different": defaultdict(list)}
                rel_groups = {"same": defaultdict(list), "different": defaultdict(list)}
                for p in all_pairs:
                    d = pair_values[p]; b = baseline[p]
                    relation = "same" if p in same else "different"
                    if not (finite(d) and finite(b)):
                        continue
                    resid = float(d) - float(b)
                    rel = resid / float(b) if float(b) > 1e-12 else math.nan
                    state_pair = tuple(sorted((st[p[0]]["normalized_arrangement_state"], st[p[1]]["normalized_arrangement_state"])))
                    categories[relation].append(resid); category_groups[relation][state_pair].append(resid)
                    if finite(rel):
                        rel_categories[relation].append(rel); rel_groups[relation][state_pair].append(rel)
                    residual_rows.append({
                        "Gene": gene, "canonical_region_id": rid, "baseline_type": baseline_type,
                        "start": r["start"], "end": r["end"], "midpoint": r["midpoint"],
                        "species1": p[0], "species2": p[1], "state1": st[p[0]]["normalized_arrangement_state"],
                        "state2": st[p[1]]["normalized_arrangement_state"], "state_relation": relation,
                        "distance": fmt(d), "baseline_median": fmt(b), "residual": fmt(resid),
                        "relative_residual": fmt(rel), "inside_focal_region": str(in_interval(mid, focal)),
                    })
                same_mean, diff_mean = fmean(categories["same"]), fmean(categories["different"])
                same_bal, diff_bal = balanced_mean(category_groups["same"]), balanced_mean(category_groups["different"])
                same_rel, diff_rel = fmean(rel_categories["same"]), fmean(rel_categories["different"])
                adjusted_rows.append({
                    "Gene": gene, "canonical_region_id": rid, "baseline_type": baseline_type,
                    "start": r["start"], "end": r["end"], "midpoint": r["midpoint"],
                    "n_same_pairs": len(categories["same"]), "n_different_pairs": len(categories["different"]),
                    "mean_same_residual": fmt(same_mean), "mean_different_residual": fmt(diff_mean),
                    "A_adjusted": fmt(diff_mean - same_mean if finite(same_mean) and finite(diff_mean) else math.nan),
                    "A_balanced": fmt(diff_bal - same_bal if finite(same_bal) and finite(diff_bal) else math.nan),
                    "n_same_state_pair_categories": len(category_groups["same"]), "n_different_state_pair_categories": len(category_groups["different"]),
                    "inside_focal_region": str(in_interval(mid, focal)), "status": "exploratory_pair_specific_residual",
                })
                relative_rows.append({
                    "Gene": gene, "canonical_region_id": rid, "baseline_type": baseline_type,
                    "start": r["start"], "end": r["end"], "midpoint": r["midpoint"],
                    "n_same_pairs": len(rel_categories["same"]), "n_different_pairs": len(rel_categories["different"]),
                    "mean_same_relative_residual": fmt(same_rel), "mean_different_relative_residual": fmt(diff_rel),
                    "A_relative_residual": fmt(diff_rel - same_rel if finite(same_rel) and finite(diff_rel) else math.nan),
                    "inside_focal_region": str(in_interval(mid, focal)), "status": "relative_pair_specific_residual",
                })

    write(RES / "pairwise_background_distances.tsv", bg_rows, list(bg_rows[0]) if bg_rows else ["canonical_region_id"])
    write(RES / "window_pairwise_residuals.tsv", residual_rows, list(residual_rows[0]) if residual_rows else ["Gene"])
    write(RES / "window_adjusted_arrangement_contrast.tsv", adjusted_rows, list(adjusted_rows[0]) if adjusted_rows else ["Gene"])
    write(RES / "window_relative_residual_contrast.tsv", relative_rows, list(relative_rows[0]) if relative_rows else ["Gene"])

    # Inside/outside summaries, with the circularity explicitly recorded.
    io_rows = []
    for rid in eligible_ids:
        reg = region_bounds[rid]
        for bt in ("broad_outside_focal", "strict_outside_all_structural_regions"):
            z = [x for x in adjusted_rows if x["canonical_region_id"] == rid and x["baseline_type"] == bt]
            zr = [x for x in relative_rows if x["canonical_region_id"] == rid and x["baseline_type"] == bt]
            for field, label, rows in [("A_adjusted", "raw_adjusted", z), ("A_balanced", "balanced_adjusted", z), ("A_relative_residual", "relative_adjusted", zr)]:
                vals_in = [float(x[field]) for x in rows if x["inside_focal_region"] == "True" and finite(x[field])]
                vals_out = [float(x[field]) for x in rows if x["inside_focal_region"] == "False" and finite(x[field])]
                io_rows.append({
                    "canonical_region_id": rid, "baseline_type": bt, "contrast": label,
                    "inside_mean": fmt(fmean(vals_in)), "outside_mean": fmt(fmean(vals_out)),
                    "delta_inside_minus_outside": fmt(fmean(vals_in) - fmean(vals_out) if vals_in and vals_out else math.nan),
                    "n_inside_windows": len(vals_in), "n_outside_windows": len(vals_out),
                    "interpretation": "outside residuals are centered by the pair-specific baseline; inside-versus-outside is descriptive and not an independent null",
                })
    write(RES / "inside_outside_adjusted_summary.tsv", io_rows, list(io_rows[0]) if io_rows else ["canonical_region_id"])

    # Leave-one-species-out sensitivity. Baselines and contrasts are recomputed,
    # rather than merely deleting values from a fixed residual table.
    loo = []
    for rid in eligible_ids:
        st = state_by_region[rid]; focal = region_bounds[rid]
        for excluded in SPECIES:
            names = [x for x in SPECIES if x != excluded and x in st]
            allp, same, diff = state_pairs(st, names)
            if not same or not diff:
                loo.append({"canonical_region_id": rid, "excluded_species": excluded, "baseline_type": "broad_outside_focal", "status": "unsupported_after_exclusion"}); continue
            for bt in ("broad_outside_focal", "strict_outside_all_structural_regions"):
                base = {}
                for p in allp:
                    vals = []
                    for g in ordered_genes:
                        mid = float(unique[g]["midpoint"])
                        if bt == "broad_outside_focal" and in_interval(mid, focal): continue
                        if bt == "strict_outside_all_structural_regions" and in_any_region(mid): continue
                        d = tree_distances[g].get(pair_key(SPECIES[p[0]], SPECIES[p[1]]))
                        if finite(d): vals.append(float(d))
                    base[p] = median(vals) if vals else math.nan
                inside, outside = [], []
                for g in ordered_genes:
                    mid = float(unique[g]["midpoint"]); dct = tree_distances[g]
                    sv, dv = [], []
                    for p in same:
                        d, b = dct.get(pair_key(SPECIES[p[0]], SPECIES[p[1]])), base[p]
                        if finite(d) and finite(b): sv.append(float(d)-float(b))
                    for p in diff:
                        d, b = dct.get(pair_key(SPECIES[p[0]], SPECIES[p[1]])), base[p]
                        if finite(d) and finite(b): dv.append(float(d)-float(b))
                    if sv and dv:
                        val = mean(dv) - mean(sv)
                        (inside if in_interval(mid, focal) else outside).append(val)
                loo.append({"canonical_region_id": rid, "excluded_species": excluded, "baseline_type": bt,
                            "inside_mean_A_adjusted": fmt(fmean(inside)), "outside_mean_A_adjusted": fmt(fmean(outside)),
                            "delta_inside_minus_outside": fmt(fmean(inside)-fmean(outside) if inside and outside else math.nan),
                            "n_inside_windows": len(inside), "n_outside_windows": len(outside), "status": "exploratory_leave_one_species_out"})
    write(RES / "leave_one_species_out.tsv", loo, list(loo[0]) if loo else ["canonical_region_id"])

    # Spatial concordance is descriptive; no independent p-value is assigned to
    # adjacent windows. The same region-specific adjusted track is compared to
    # each source-defined QQS design without duplicating the divergence track.
    concord = []
    for rid in eligible_ids:
        for bt in ("broad_outside_focal", "strict_outside_all_structural_regions"):
            arows = {x["Gene"]: x for x in adjusted_rows if x["canonical_region_id"] == rid and x["baseline_type"] == bt}
            for cl in CLADES:
                q = [x for x in support if x["clade"] == cl and x["Gene"] in arows]
                for qfield, qlabel in [("D_q2_minus_q1", "D2"), ("D_q3_minus_q1", "D3")]:
                    pairs = [(float(arows[x["Gene"]]["A_adjusted"]), float(x[qfield])) for x in q if finite(arows[x["Gene"]]["A_adjusted"]) and finite(x[qfield])]
                    corr = np.corrcoef(np.array(pairs).T)[0, 1] if len(pairs) >= 3 and np.std(np.array(pairs)[:,0]) > 0 and np.std(np.array(pairs)[:,1]) > 0 else math.nan
                    inside = [float(x[qfield]) for x in q if arows[x["Gene"]]["inside_focal_region"] == "True" and finite(x[qfield])]
                    outside = [float(x[qfield]) for x in q if arows[x["Gene"]]["inside_focal_region"] == "False" and finite(x[qfield])]
                    concord.append({"canonical_region_id": rid, "baseline_type": bt, "clade": cl, "topology_contrast": qlabel,
                                    "n_windows": len(pairs), "pearson_window_correlation": fmt(corr),
                                    "inside_topology_mean": fmt(fmean(inside)), "outside_topology_mean": fmt(fmean(outside)),
                                    "note": "descriptive spatial concordance; adjacent loci are dependent and no independent significance claim is made"})
    write(RES / "spatial_concordance_summary.tsv", concord, list(concord[0]) if concord else ["canonical_region_id"])

    # Raw Stage 7 reproduction validation for the candidate regions.
    old = { (r["Gene"], r["canonical_region_id"]): r for r in read(STAGE7 / "window_arrangement_distance.tsv") }
    raw_validation = []
    for rid in eligible_ids:
        st = state_by_region[rid]; allp, same, diff = state_pairs(st)
        for g in ordered_genes:
            dct = tree_distances[g]; sv = [dct.get(pair_key(SPECIES[a], SPECIES[b])) for a,b in same]; dv = [dct.get(pair_key(SPECIES[a], SPECIES[b])) for a,b in diff]
            sv = [float(x) for x in sv if finite(x)]; dv = [float(x) for x in dv if finite(x)]
            newa = fmean(dv)-fmean(sv) if sv and dv else math.nan
            olda = old.get((g, rid), {}).get("A_distance", math.nan)
            raw_validation.append({"Gene": g, "canonical_region_id": rid, "stage7_A_distance": fmt(olda), "stage8_recomputed_A_distance": fmt(newa), "absolute_difference": fmt(abs(float(olda)-newa) if finite(olda) and finite(newa) else math.nan), "status": "raw_stage7_reproduction"})
    write(RES / "raw_stage7_reproduction.tsv", raw_validation, list(raw_validation[0]) if raw_validation else ["Gene"])

    write_docs(eligibility, len(tree_distances), len(residual_rows), parse_errors)
    make_figures(support, regions, adjusted_rows, residual_rows)
    inputs = [LOCUS, MANIFEST, STATES, REGIONS, SUPPORT, TREES, ROLES, J48]
    output_files = list(RES.glob("*.tsv")) + list(RES.glob("*.md"))
    manifest = {
        "stage": "Neoaves Stage 8 phylogenetically controlled arrangement divergence",
        "script": str(Path(__file__).relative_to(ROOT)),
        "command": "MPLBACKEND=Agg python empirical/neoaves_chr4/scripts/08_phylogenetic_distance_adjustment.py",
        "inputs_sha256": {str(p.relative_to(ROOT)): sha(p) for p in inputs if p.exists()},
        "n_source_trees_with_six_species": len(tree_distances), "n_residual_rows": len(residual_rows),
        "parse_errors": parse_errors, "eligible_regions": eligible_ids,
        "structural_taxon_count": len(SPECIES), "additional_taxa_added": 0,
        "baseline_types": ["broad_outside_focal", "strict_outside_all_structural_regions"],
        "interpretation": "exploratory species-level multistate synteny contrast; not calibrated time, coalescence time, or recombination-rate inference",
        "output_sha256": {str(p.relative_to(ROOT)): sha(p) for p in output_files if p.exists()},
    }
    (RES / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def write_docs(eligibility, ntrees, nresid, parse_errors):
    cand = [x["canonical_region_id"] for x in eligibility if x["eligible"] == "True"]
    (RES / "validation_report.md").write_text(
        "# Validation report\n\n"
        f"Parsed {ntrees} chromosome-4 local ML trees with at least the six exact structural taxa; {parse_errors} source records had parser errors. The Stage 7 eligibility rules reproduce candidate regions {', '.join(cand) or 'none'}. The Stage 8 residual table contains {nresid} pair-region-window rows.\n\n"
        "Raw Stage 7 A-distance reproduction is written to `raw_stage7_reproduction.tsv`; it is used as a regression audit and does not overwrite Stage 7. Pair-specific baselines are calculated from windows outside the focal region (broad) or outside all canonical regions (strict).\n\n"
        "The outside baseline is estimated from the same chromosome-4 local ML distance series used to form residuals. Therefore outside residuals are centered by construction and inside-versus-outside summaries are descriptive, not independent tests of the baseline model.\n"
    )
    (RES / "methods_text.md").write_text(
        "# Methods\n\n"
        "For each eligible canonical structural region, the six exact Stage 7 structural species were paired using their independently supplied region-specific multistate synteny labels. For each species pair, the background distance was the median local ML patristic distance in windows outside the focal region (broad baseline) or outside every canonical structural region (strict baseline). Pair-specific residuals were `R_ij(w)=d_ij(w)-d0_ij`. The adjusted arrangement contrast was the mean residual for different-state pairs minus the mean residual for same-state pairs. Relative residuals divide each residual by its positive pair-specific baseline. Equal-weight state-pair-category summaries are retained as a sensitivity analysis.\n\n"
        "Distances are in ML substitution units. They are not calibrated divergence times or realized coalescence times. Pairwise observations share species and loci, so no iid pairwise standard errors are reported. Leave-one-species-out recomputation quantifies sensitivity to the limited six-species sample.\n"
    )
    (RES / "results_text.md").write_text(
        "# Results\n\n"
        f"The independent structural-state sample could not be expanded beyond six exact local-tree taxa without unsupported state imputation. The recomputed Stage 7 eligibility set contains {', '.join(cand) or 'no'} eligible exploratory regions. Adjusted contrasts are reported for CR_01, CR_10, and CR_11 where the same-versus-different state comparison is numerically estimable. Because CR_01 has only one same-state pair and CR_10/CR_11 pairs share the same six species, these are exploratory species-level residual contrasts rather than independently replicated evidence for arrangement-associated divergence.\n\n"
        "The three QQS designs remain separate in spatial concordance summaries. A common six-species divergence track is not relabeled as three focal-design-specific A(w) tracks.\n"
    )
    (RES / "report.md").write_text(
        "# Neoaves chromosome 4 Stage 8 report\n\n"
        "## Scientific scope\n\n"
        "Stage 8 asks whether independently supplied structural states explain local ML-tree patristic distances beyond each species pair's chromosome-4 background distance. It does not estimate recombination rates or coalescence times.\n\n"
        "## Findings\n\n"
        f"The structural-state sample remains six species; no additional labels were imputed. The region-specific multistate synteny table supports exploratory comparisons for {', '.join(cand) or 'no regions'}. Pair-specific broad and strict outside baselines, residual contrasts, relative residuals, and leave-one-species-out sensitivity are provided.\n\n"
        "The adjusted statistic is related in purpose to the arrangement-associated contrasts in other datasets, but it is not directly analogous in inferential strength: Neoaves has deep species-level comparisons, multistate states, and limited same-state replication. The topology tracks D2 and D3 remain source-defined and are compared descriptively to the one common six-species divergence analysis.\n\n"
        "## Interpretation\n\n"
        "Any large adjusted signal must be checked against pair composition, leave-one-species-out results, long branches, missing loci, and the broad-versus-strict baseline. A residual contrast that persists numerically does not resolve the lack of independent structural-state replication.\n"
    )


def make_figures(support, regions, adjusted, residuals):
    FIG.mkdir(parents=True, exist_ok=True)
    # Main topology overview, using the frozen Stage 6 source-defined designs.
    fig, axes = plt.subplots(7, 1, figsize=(13, 13), sharex=True, constrained_layout=True)
    colors = {"q1_species_history": "#4C78A8", "q2_alternative_1": "#E45756", "q3_alternative_2": "#72B7B2"}
    for idx, cl in enumerate(CLADES):
        q = sorted([x for x in support if x["clade"] == cl], key=lambda x: float(x["midpoint"]))
        x = np.array([float(x["midpoint"]) / 1e6 for x in q])
        labels = ["q1 species history", "q2 alternative", "q3 alternative"] if idx == 0 else ["", "", ""]
        axes[2*idx].stackplot(x, *[np.array([float(z[k]) for z in q]) for k in colors], labels=labels, colors=list(colors.values()), alpha=.82, linewidth=0)
        axes[2*idx].set_ylim(0, 1); axes[2*idx].set_ylabel(f"{cl}\n support")
        axes[2*idx+1].plot(x, [float(z["D_q2_minus_q1"]) for z in q], color="#E45756", lw=.7, label="D2 = q2−q1" if idx == 0 else None)
        axes[2*idx+1].plot(x, [float(z["D_q3_minus_q1"]) for z in q], color="#72B7B2", lw=.7, label="D3 = q3−q1" if idx == 0 else None)
        axes[2*idx+1].axhline(0, color="0.35", lw=.5); axes[2*idx+1].set_ylabel(f"{cl}\n D")
    for ax in axes:
        for r in regions:
            ax.axvspan(float(r["canonical_start"])/1e6, float(r["canonical_end"])/1e6, color="0.85", alpha=.14)
    axes[0].legend(loc="upper right", ncol=3, fontsize=8); axes[1].legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("chromosome 4 position (Mb)")
    fig.savefig(FIG / "neoaves_chr4_stage8_topology_overview.pdf"); fig.savefig(FIG / "neoaves_chr4_stage8_topology_overview.png", dpi=220); plt.close(fig)

    # Adjusted divergence: one common species-level track, not repeated per QQS design.
    fig, axes = plt.subplots(3, 1, figsize=(13, 8), sharex=True, constrained_layout=True)
    for ax, rid in zip(axes, ["CR_01", "CR_10", "CR_11"]):
        for bt, color, label in [("broad_outside_focal", "#4C78A8", "broad residual"), ("strict_outside_all_structural_regions", "#F58518", "strict residual")]:
            rows = sorted([x for x in adjusted if x["canonical_region_id"] == rid and x["baseline_type"] == bt], key=lambda z: float(z["midpoint"]))
            ax.plot([float(z["midpoint"])/1e6 for z in rows], [float(z["A_adjusted"]) if finite(z["A_adjusted"]) else np.nan for z in rows], lw=.65, color=color, label=label)
        rr = sorted([x for x in adjusted if x["canonical_region_id"] == rid and x["baseline_type"] == "broad_outside_focal"], key=lambda z: float(z["midpoint"]))
        ax.plot([float(z["midpoint"])/1e6 for z in rr], [float(z["A_balanced"]) if finite(z["A_balanced"]) else np.nan for z in rr], lw=.55, ls="--", color="#54A24B", label="broad state-pair balanced")
        reg = next(r for r in regions if r["canonical_region_id"] == rid)
        ax.axvspan(float(reg["canonical_start"])/1e6, float(reg["canonical_end"])/1e6, color="0.85", alpha=.3)
        ax.axhline(0, color="0.35", lw=.5); ax.set_ylabel(f"{rid}\nA residual")
    axes[0].legend(loc="upper right", fontsize=8); axes[-1].set_xlabel("chromosome 4 position (Mb)")
    fig.savefig(FIG / "neoaves_chr4_stage8_adjusted_divergence.pdf"); fig.savefig(FIG / "neoaves_chr4_stage8_adjusted_divergence.png", dpi=220); plt.close(fig)

    # Pair diagnostics expose species-pair influence.
    fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True, constrained_layout=True)
    for ax, rid in zip(axes, ["CR_01", "CR_10", "CR_11"]):
        rows = [x for x in residuals if x["canonical_region_id"] == rid and x["baseline_type"] == "broad_outside_focal"]
        labels = sorted({f"{x['species1']}-{x['species2']}" for x in rows})
        cmap = plt.get_cmap("tab20")
        for j, label in enumerate(labels):
            z = [x for x in rows if f"{x['species1']}-{x['species2']}" == label]
            ax.scatter([float(x["midpoint"])/1e6 for x in z], [float(x["residual"]) for x in z], s=2, alpha=.18, color=cmap(j % 20), label=label)
        reg = next(r for r in regions if r["canonical_region_id"] == rid)
        ax.axvspan(float(reg["canonical_start"])/1e6, float(reg["canonical_end"])/1e6, color="0.85", alpha=.3)
        ax.axhline(0, color="0.35", lw=.5); ax.set_ylabel(f"{rid}\npair residual")
    axes[-1].set_xlabel("chromosome 4 position (Mb)")
    axes[0].legend(loc="upper right", fontsize=6, ncol=3, markerscale=3)
    fig.savefig(FIG / "neoaves_chr4_stage8_pair_diagnostics.pdf"); fig.savefig(FIG / "neoaves_chr4_stage8_pair_diagnostics.png", dpi=220); plt.close(fig)


if __name__ == "__main__":
    main()
