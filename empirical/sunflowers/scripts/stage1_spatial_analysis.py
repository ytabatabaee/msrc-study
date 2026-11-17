#!/usr/bin/env python3
"""Arrangement-blind Stage 1 audit of sunflower local gene trees.

The script scores fixed species quartets in the deposited Newicks and retains
the observed coordinate gaps. Rounded published intervals are shrunk by
500 kb at each end for the conservative interior analysis.
"""
from __future__ import annotations

import csv
import math
import os
import re
from collections import Counter, defaultdict
from io import StringIO
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc-paper-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data/sunflowers/raw"
EMP = ROOT / "empirical/sunflowers"
RES = EMP / "results"
FIG = EMP / "figures"
TREE_FILE = RAW / "10.2kbphysical.trees.txt"
WINDOW_FILE = RAW / "10kb_window_sizes.txt"
BLOCK_FILE = RES / "haploblocks_reported.tsv"
EDGE_UNCERTAINTY = 500_000
FLANK_BP = 5_000_000
ADJACENCY_BP = 20_000
BLOCKS = ("ann05.01", "ann13.01")
Q_NAMES = ("q1_ann_arg", "q2_ann_pet", "q3_arg_pet")
COLORS = {"q1_ann_arg": "#276FBF", "q2_ann_pet": "#D1495B", "q3_arg_pet": "#4C956C"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_blocks() -> dict[str, dict]:
    rows = read_tsv(BLOCK_FILE)
    result = {}
    for row in rows:
        if row["haploblock_id"] not in BLOCKS:
            continue
        parts = []
        for interval in row["reported_regions_mbp"].split(","):
            left, right = map(float, interval.strip().split("-"))
            parts.append((round(left * 1_000_000), round(right * 1_000_000)))
        result[row["haploblock_id"]] = {**row, "parts": parts}
    if set(result) != set(BLOCKS):
        raise ValueError(f"Missing focal block rows: {set(BLOCKS) - set(result)}")
    return result


def tip_group(name: str) -> str | None:
    if name.startswith("ANN"):
        return "ann"
    if name.startswith("ARG"):
        return "arg"
    if name.startswith("PET"):
        return "pet"
    if name.startswith("DEB"):
        return "deb"
    if name.startswith("PAR"):
        return "par"
    return None


def quartet_combinations(tips: list[str], fourth_group: str = "deb") -> list[tuple[str, str, str, str]]:
    groups = {g: sorted(t for t in tips if tip_group(t) == g) for g in ("ann", "arg", "pet", fourth_group)}
    if any(not groups[g] for g in groups):
        raise ValueError(f"A tree lacks one of the fixed species groups: {groups}")
    return [(a, b, p, d) for a in groups["ann"] for b in groups["arg"]
            for p in groups["pet"] for d in groups[fourth_group]]


def score_quartet_distances(d_ab: float, d_ap: float, d_ad: float,
                            d_bp: float, d_bd: float, d_pd: float) -> tuple[int | None, float, float]:
    """Return fixed split index, internal edge, and q2-minus-q1 paired distance."""
    sums = (d_ab + d_pd, d_ap + d_bd, d_ad + d_bp)
    order = sorted(range(3), key=lambda i: sums[i])
    if abs(sums[order[1]] - sums[order[0]]) <= 1e-12:
        return None, 0.0, (sums[1] - sums[0]) / 2.0
    split = order[0]
    internal = max(0.0, (sum(sums) - 3.0 * sums[split]) / 4.0)
    return split, internal, (sums[1] - sums[0]) / 2.0


def mean_pairwise_distance(values: list[float]) -> float:
    return float(np.mean(values)) if values else math.nan


def score_tree(newick: str, fourth_group: str = "deb") -> dict:
    tree = Phylo.read(StringIO(newick), "newick")
    terminals = tree.get_terminals()
    tips = [t.name for t in terminals]
    if len(tips) != 22 or len(set(tips)) != len(tips):
        raise ValueError(f"Expected 22 unique tips, found {len(tips)}")
    branches = [clade.branch_length for clade in tree.find_clades() if clade != tree.root]
    if any(v is None or not math.isfinite(v) or v < 0 for v in branches):
        raise ValueError("Missing, negative, or nonfinite branch length")
    combos = quartet_combinations(tips, fourth_group)
    dist_cache: dict[tuple[str, str], float] = {}

    def dist(x: str, y: str) -> float:
        key = tuple(sorted((x, y)))
        if key not in dist_cache:
            dist_cache[key] = float(tree.distance(x, y))
        return dist_cache[key]

    topo = Counter()
    unresolved = 0
    internal = []
    a_values = []
    pairwise_values = []
    # All combinations of the sampled individuals in four fixed species groups;
    # each tree contributes one within-window average, not 48 independent loci.
    for a, b, p, d in combos:
        d_ab, d_ap, d_ad = dist(a, b), dist(a, p), dist(a, d)
        d_bp, d_bd, d_pd = dist(b, p), dist(b, d), dist(p, d)
        sums = (d_ab + d_pd, d_ap + d_bd, d_ad + d_bp)
        split, internal_edge, a_value = score_quartet_distances(d_ab, d_ap, d_ad, d_bp, d_bd, d_pd)
        if split is None:
            unresolved += 1
            continue
        name = Q_NAMES[split]
        topo[name] += 1
        internal.append(internal_edge)
        # A(w)-like, species-fixed paired-distance contrast. This is explicitly
        # not the arrangement-defined A(w), because tip genotypes are unknown.
        a_values.append(a_value)
        pairwise_values.extend((d_ab, d_ap, d_ad, d_bp, d_bd, d_pd))
    n_resolved = sum(topo.values())
    denom = len(combos)
    q = {k: topo[k] / denom for k in Q_NAMES}
    entropy = -sum(x * math.log(x) for x in q.values() if x > 0) / math.log(3)
    maxima = max(q.values())
    dominant = Q_NAMES[[i for i, k in enumerate(Q_NAMES) if abs(q[k] - maxima) < 1e-12][0]] if maxima > 0 else "unresolved"
    return {
        **q,
        "q_unresolved": unresolved / denom,
        "n_induced_quartets": denom,
        "n_resolved_quartets": n_resolved,
        "D_q2_minus_q1": q["q2_ann_pet"] - q["q1_ann_arg"],
        "discordant_frequency": q["q2_ann_pet"] + q["q3_arg_pet"],
        "topology_entropy_normalized": entropy,
        "dominant_topology": dominant,
        "A_species_pair_distance_q2_minus_q1": float(np.mean(a_values)) if a_values else math.nan,
        "B_mean_pairwise_distance": mean_pairwise_distance(pairwise_values),
        "mean_quartet_internal_branch": float(np.mean(internal)) if internal else math.nan,
        "n_tips": len(tips),
        "n_branch_lengths": len(branches),
    }


def load_trees() -> tuple[list[dict], dict[str, int], dict[str, tuple[int, int, int]]]:
    called = {}
    manifest = {}
    for line in WINDOW_FILE.read_text().splitlines():
        name, n = line.split()
        key = name[:-6] if name.endswith(".fasta") else name
        called[key] = int(n)
        m = re.fullmatch(r"Ha412HOChr(\d+)-(\d+)-(\d+)", key)
        if m:
            manifest[key] = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
    trees = []
    seen = set()
    reference_tips = None
    for line in TREE_FILE.read_text().splitlines():
        key, newick = line.split("\t", 1)
        m = re.fullmatch(r"Ha412HOChr(\d+)-(\d+)-(\d+)", key)
        if not m:
            raise ValueError(f"Bad coordinate key {key}")
        chrom, start, end = map(int, m.groups())
        if key in seen or end <= start:
            raise ValueError(f"Duplicate or invalid interval {key}")
        seen.add(key)
        if key not in called:
            raise ValueError(f"Tree absent from called-base manifest: {key}")
        scores = score_tree(newick)
        par_scores = score_tree(newick, "par")
        this_tips = tuple(sorted(t.name for t in Phylo.read(StringIO(newick), "newick").get_terminals()))
        if reference_tips is None:
            reference_tips = this_tips
        elif this_tips != reference_tips:
            raise ValueError(f"Taxon labels differ at {key}")
        trees.append({"tree_id": key, "chrom": chrom, "start": start, "end": end,
                      "mid": (start + end) / 2, "called_bases": called[key], **scores,
                      **{"par_"+k:v for k,v in par_scores.items() if k in (*Q_NAMES,"D_q2_minus_q1","discordant_frequency","topology_entropy_normalized","dominant_topology","A_species_pair_distance_q2_minus_q1","B_mean_pairwise_distance","mean_quartet_internal_branch","n_induced_quartets","n_resolved_quartets")}})
    trees.sort(key=lambda r: (r["chrom"], r["start"], r["end"]))
    return trees, called, manifest


def tree_qc(trees: list[dict], called: dict[str, int], manifest: dict[str, tuple[int, int, int]]) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    tree_ids = {t["tree_id"] for t in trees}
    missing = []
    for key, (chrom, start, end) in manifest.items():
        if called[key] >= 2_000 and key not in tree_ids:
            missing.append({"tree_id": key, "chromosome": chrom, "start_bp": start, "end_bp": end,
                            "called_bases": called[key], "status": "qualifying manifest window without deposited Newick"})
    overlap_pairs = 0
    overlaps = []
    per_chrom = defaultdict(list)
    for t in trees:
        per_chrom[t["chrom"]].append(t)
    chrom_counts = []
    for chrom, rows in sorted(per_chrom.items()):
        rows.sort(key=lambda r: (r["start"], r["end"]))
        n_overlap = 0
        active = []
        for row in rows:
            active = [prior for prior in active if prior["end"] > row["start"]]
            n_overlap += len(active)
            for prior in active:
                overlaps.append({"tree_id_1": prior["tree_id"], "tree_id_2": row["tree_id"],
                                 "chromosome": chrom, "overlap_start_bp": row["start"],
                                 "overlap_end_bp": min(prior["end"],row["end"]),
                                 "overlap_bp": min(prior["end"],row["end"])-row["start"]})
            active.append(row)
        overlap_pairs += n_overlap
        manifest_chrom = [(s, e) for c, s, e in manifest.values() if c == chrom and called.get("Ha412HOChr%02d-%d-%d" % (c,s,e),0) >= 2_000]
        chrom_counts.append({"scope": "chromosome", "chromosome": chrom, "n_trees": len(rows),
                             "n_qualifying_manifest_windows": len(manifest_chrom),
                             "n_qualifying_manifest_windows_without_tree": sum(1 for m in missing if m["chromosome"]==chrom),
                             "n_overlapping_tree_interval_pairs": n_overlap,
                             "observed_coordinate_span_bp": max(r["end"] for r in rows)-min(r["start"] for r in rows)})
    summary = [{"scope": "all_chromosomes", "chromosome": "all", "n_trees": len(trees),
                "n_unique_tree_ids": len(tree_ids), "n_tip_count_values": len(set(t["n_tips"] for t in trees)),
                "n_branch_count_values": len(set(t["n_branch_lengths"] for t in trees)),
                "n_trees_with_nonnegative_finite_branches": sum(1 for t in trees if t["n_branch_lengths"] >= 22),
                "n_qualifying_manifest_windows": sum(1 for k in manifest if called[k] >= 2_000),
                "n_qualifying_manifest_windows_without_tree": len(missing),
                "n_overlapping_tree_interval_pairs": overlap_pairs,
                "taxon_labels_identical_across_trees": "yes"}]
    return summary, chrom_counts, missing, overlaps


def overlap_bp(start: int, end: int, intervals: list[tuple[int, int]]) -> int:
    clipped = sorted((max(start, a), min(end, b)) for a, b in intervals if start < b and end > a)
    total = 0
    last = None
    for a, b in clipped:
        if last is None or a > last:
            total += b - a
            last = b
        elif b > last:
            total += b - last
            last = b
    return total


def calc_coverage(trees: list[dict], blocks: dict[str, dict]) -> tuple[list[dict], list[dict]]:
    coords, coverage = [], []
    for block_id, block in blocks.items():
        chrom = int(block["chromosome"])
        raw_parts = block["parts"]
        core_parts = [(a + EDGE_UNCERTAINTY, b - EDGE_UNCERTAINTY) for a, b in raw_parts if b - a > 2 * EDGE_UNCERTAINTY]
        uncertainty_parts = [(a - EDGE_UNCERTAINTY, a + EDGE_UNCERTAINTY) for a, _ in raw_parts]
        uncertainty_parts += [(b - EDGE_UNCERTAINTY, b + EDGE_UNCERTAINTY) for _, b in raw_parts]
        span_left, span_right = min(a for a, _ in raw_parts), max(b for _, b in raw_parts)
        flanks = [(max(0, span_left - FLANK_BP), span_left), (span_right, span_right + FLANK_BP)]
        selected = [t for t in trees if t["chrom"] == chrom]
        region_length = sum(b - a for a, b in raw_parts)
        tree_clips = []
        for t in selected:
            tree_clips.extend((max(t["start"], a), min(t["end"], b)) for a,b in raw_parts if t["start"] < b and t["end"] > a)
        tree_clips.sort()
        union_cover = 0
        covered_to = None
        for a,b in tree_clips:
            if covered_to is None or a > covered_to:
                union_cover += b-a; covered_to=b
            elif b > covered_to:
                union_cover += b-covered_to; covered_to=b
        inside = [t for t in selected if overlap_bp(t["start"], t["end"], raw_parts)]
        core = [t for t in selected if any(t["start"] >= a and t["end"] <= b for a, b in core_parts)]
        uncertain = [t for t in selected if overlap_bp(t["start"], t["end"], uncertainty_parts)]
        flank_counts = []
        for a, b in flanks:
            flank_counts.append(sum(1 for t in selected if t["start"] >= a and t["end"] <= b))
        coords.append({"haploblock_id": block_id, "chromosome": chrom,
                       "reported_interval_mbp": block["reported_regions_mbp"],
                       "interval_meaning": "published haploblock region; not fine-mapped inversion breakpoints",
                       "support_codes": block["support_codes"],
                       "genetic_map_suppression_marked": "g" in {x.strip() for x in block["support_codes"].split(",")},
                       "breakpoint_precision": "rounded to nearest Mbp; +/-0.5 Mb at each endpoint used as sensitivity band",
                       "interior_used_for_primary_track": "each reported segment shrunk 0.5 Mb at both ends",
                       "n_parts": len(raw_parts), "reported_length_bp": region_length,
                       "boundary_uncertainty_bp_each_end": EDGE_UNCERTAINTY})
        coverage.append({"haploblock_id": block_id, "reported_tree_overlaps": len(inside),
                         "conservative_interior_trees": len(core), "boundary_uncertain_trees": len(uncertain),
                         "left_5Mb_flank_trees": flank_counts[0], "right_5Mb_flank_trees": flank_counts[1],
                         "informative_tree_span_bp_union": union_cover, "reported_interval_bp": region_length,
                         "fraction_interval_spanned_by_tree_intervals": union_cover / region_length,
                         "mean_called_bases_per_overlapping_tree": float(np.mean([t["called_bases"] for t in inside])) if inside else "NA",
                         "overlapping_tree_called_bases_sum": sum(t["called_bases"] for t in inside),
                         "interior_called_bases_sum": sum(t["called_bases"] for t in core),
                         "tree_density_per_10Mb_reported_interval": len(inside) / (region_length / 10_000_000),
                         "flank_length_each_bp": FLANK_BP})
    return coords, coverage


def in_block(tree: dict, block: dict) -> bool:
    return tree["chrom"] == int(block["chromosome"]) and overlap_bp(tree["start"], tree["end"], block["parts"]) > 0


def in_core(tree: dict, block: dict) -> bool:
    return tree["chrom"] == int(block["chromosome"]) and any(tree["start"] >= a + EDGE_UNCERTAINTY and tree["end"] <= b - EDGE_UNCERTAINTY for a, b in block["parts"])


def in_flank(tree: dict, block: dict) -> bool:
    if tree["chrom"] != int(block["chromosome"]):
        return False
    a = min(x[0] for x in block["parts"]); b = max(x[1] for x in block["parts"])
    return (tree["start"] >= max(0, a - FLANK_BP) and tree["end"] <= a) or (tree["start"] >= b and tree["end"] <= b + FLANK_BP)


def adjacent(a: dict, b: dict) -> bool:
    return 0 <= b["start"] - a["end"] <= ADJACENCY_BP


def detect_runs(rows: list[dict], block_id: str, region: str) -> list[dict]:
    out = []
    for q in Q_NAMES:
        current = []
        for row in rows:
            unique_max = sum(1 for key in Q_NAMES if abs(row[key] - max(row[k] for k in Q_NAMES)) <= 1e-12) == 1
            if row["dominant_topology"] == q and row[q] == max(row[k] for k in Q_NAMES) and unique_max:
                if current and adjacent(current[-1], row):
                    current.append(row)
                else:
                    if current:
                        out.append(run_row(block_id, region, q, current))
                    current = [row]
            else:
                if current:
                    out.append(run_row(block_id, region, q, current)); current = []
        if current:
            out.append(run_row(block_id, region, q, current))
    for q in Q_NAMES:
        n = sum(1 for r in out if r["haploblock_id"] == block_id and r["region"] == region and r["topology"] == q)
        if not n:
            out.append({"haploblock_id": block_id, "region": region, "topology": q, "n_adjacent_windows": 0,
                        "start_bp": "NA", "end_bp": "NA", "span_bp": 0,
                        "adjacency_rule": f"inter-tree gap <= {ADJACENCY_BP} bp", "n_independent_loci": 0,
                        "biological_independence_accounted_for": "no; none observed under adjacency rule"})
    return out


def run_row(block_id: str, region: str, q: str, rows: list[dict]) -> dict:
    return {"haploblock_id": block_id, "region": region, "topology": q, "n_adjacent_windows": len(rows),
            "start_bp": rows[0]["start"], "end_bp": rows[-1]["end"],
            "span_bp": rows[-1]["end"] - rows[0]["start"], "adjacency_rule": f"inter-tree gap <= {ADJACENCY_BP} bp",
            "n_independent_loci": len(rows), "biological_independence_accounted_for": "no; windows are linked"}


def autocorrelation(rows: list[dict], block_id: str) -> list[dict]:
    bins = [(0, 50_000), (50_000, 200_000), (200_000, 1_000_000), (1_000_000, 5_000_000), (5_000_000, 50_000_000)]
    result = []
    usable = [r for r in rows if math.isfinite(r["D_q2_minus_q1"])]
    for lo, hi in bins:
        x, y = [], []
        for i, a in enumerate(usable):
            for b in usable[i+1:]:
                d = b["mid"] - a["mid"]
                if lo <= d < hi:
                    x.append(a["D_q2_minus_q1"]); y.append(b["D_q2_minus_q1"])
        corr = float(np.corrcoef(x, y)[0, 1]) if len(x) >= 3 and np.std(x) and np.std(y) else "NA"
        result.append({"haploblock_id": block_id, "distance_min_bp": lo, "distance_max_bp_exclusive": hi,
                       "n_tree_pairs": len(x), "pearson_correlation_D": corr,
                       "interpretation": "descriptive pairwise autocorrelation; overlapping pairs are not independent"})
    return result


def mean_stats(rows: list[dict], block_id: str, region: str) -> dict:
    out = {"haploblock_id": block_id, "region": region, "n_windows": len(rows)}
    for key in (*Q_NAMES, "D_q2_minus_q1", "discordant_frequency", "topology_entropy_normalized",
                "A_species_pair_distance_q2_minus_q1", "B_mean_pairwise_distance", "mean_quartet_internal_branch"):
        vals = [r[key] for r in rows if isinstance(r[key], (float, int)) and math.isfinite(r[key])]
        out[f"mean_{key}"] = float(np.mean(vals)) if vals else "NA"
        out[f"median_{key}"] = float(np.median(vals)) if vals else "NA"
    return out


def build_region_outputs(trees: list[dict], blocks: dict[str, dict]):
    tracks, summaries, runs, autocorr = [], [], [], []
    regions_by_block = {}
    for block_id, block in blocks.items():
        local = [t for t in trees if t["chrom"] == int(block["chromosome"])]
        groups = {"interior": [t for t in local if in_core(t, block)],
                  "left_flank": [t for t in local if t["start"] >= max(0, min(a for a, _ in block["parts"]) - FLANK_BP) and t["end"] <= min(a for a, _ in block["parts"])],
                  "right_flank": [t for t in local if t["start"] >= max(b for _, b in block["parts"]) and t["end"] <= max(b for _, b in block["parts"]) + FLANK_BP],
                  "full_reported_interval": [t for t in local if in_block(t, block)]}
        regions_by_block[block_id] = groups
        for region, rows in groups.items():
            summaries.append(mean_stats(rows, block_id, region))
            for row in rows:
                tracks.append({"tree_id": row["tree_id"], "chromosome": row["chrom"], "start_bp": row["start"], "end_bp": row["end"],
                               "midpoint_bp": row["mid"], "called_bases": row["called_bases"], "haploblock_id": block_id, "region": region,
                               **{k: row[k] for k in (*Q_NAMES, "q_unresolved", "D_q2_minus_q1", "discordant_frequency", "topology_entropy_normalized", "dominant_topology", "A_species_pair_distance_q2_minus_q1", "B_mean_pairwise_distance", "mean_quartet_internal_branch", "n_induced_quartets", "n_resolved_quartets")}})
            runs.extend(detect_runs(rows, block_id, region))
        autocorr.extend(autocorrelation(groups["interior"], block_id))
    return tracks, summaries, runs, autocorr, regions_by_block


def boundary_sensitivity(trees: list[dict], blocks: dict[str, dict]) -> list[dict]:
    result=[]
    for block_id,block in blocks.items():
        local=[t for t in trees if t["chrom"]==int(block["chromosome"])]
        for uncertainty in (0,500_000,1_000_000):
            selected=[t for t in local if any(t["start"] >= a+uncertainty and t["end"] <= b-uncertainty for a,b in block["parts"])]
            row={"haploblock_id":block_id,"endpoint_uncertainty_bp":uncertainty,"n_interior_trees":len(selected),
                 "mean_called_bases":float(np.mean([t["called_bases"] for t in selected])) if selected else "NA"}
            for q in Q_NAMES:
                row["mean_"+q]=float(np.mean([t[q] for t in selected])) if selected else "NA"
            row["mean_D_q2_minus_q1"]=float(np.mean([t["D_q2_minus_q1"] for t in selected])) if selected else "NA"
            row["mean_B_subs_per_site"]=float(np.mean([t["B_mean_pairwise_distance"] for t in selected])) if selected else "NA"
            result.append(row)
    return result


def bin_sensitivity(trees: list[dict], blocks: dict[str, dict]) -> list[dict]:
    result=[]
    for block_id,block in blocks.items():
        core=[t for t in trees if in_core(t,block)]
        for size in (1_000_000,5_000_000,10_000_000):
            bins=defaultdict(list)
            for t in core: bins[int(t["mid"]//size)].append(t)
            d_bin=[float(np.mean([t["D_q2_minus_q1"] for t in rows])) for rows in bins.values()]
            b_bin=[float(np.mean([t["B_mean_pairwise_distance"] for t in rows])) for rows in bins.values()]
            result.append({"haploblock_id":block_id,"bin_size_bp":size,"n_nonempty_bins":len(bins),
                           "min_trees_per_nonempty_bin":min((len(x) for x in bins.values()),default=0),
                           "max_trees_per_nonempty_bin":max((len(x) for x in bins.values()),default=0),
                           "equal_bin_mean_D_q2_minus_q1":float(np.mean(d_bin)) if d_bin else "NA",
                           "window_weighted_mean_D_q2_minus_q1":float(np.mean([t["D_q2_minus_q1"] for t in core])) if core else "NA",
                           "equal_bin_mean_B_subs_per_site":float(np.mean(b_bin)) if b_bin else "NA",
                           "empty_bins_omitted": "yes; no interpolation"})
    return result


def transition_summary(trees: list[dict], blocks: dict[str, dict]) -> list[dict]:
    result=[]
    for block_id,block in blocks.items():
        rows=sorted([t for t in trees if in_core(t,block)],key=lambda r:r["start"])
        adjacent_pairs=[]; components=[]; current=[]
        for a,b in zip(rows,rows[1:]):
            if adjacent(a,b):
                adjacent_pairs.append((a,b))
                if not current: current=[a,b]
                else: current.append(b)
            elif current:
                components.append(current);current=[]
        if current: components.append(current)
        transitions=0
        changes=Counter()
        for a,b in adjacent_pairs:
            da,db=a["dominant_topology"],b["dominant_topology"]
            if da in Q_NAMES and db in Q_NAMES and da!=db:
                transitions+=1;changes[da+"->"+db]+=1
        span=sum(c[-1]["end"]-c[0]["start"] for c in components)
        result.append({"haploblock_id":block_id,"n_interior_windows":len(rows),"n_conservative_adjacent_pairs":len(adjacent_pairs),
                       "n_strict_dominant_topology_transitions":transitions,"adjacent_component_span_bp":span,
                       "transitions_per_Mb_of_observed_adjacent_span":transitions/(span/1e6) if span else "NA",
                       "transition_direction_counts": ";".join("%s:%d"%(k,v) for k,v in sorted(changes.items())) or "none",
                       "adjacency_rule":f"inter-tree gap <= {ADJACENCY_BP} bp"})
    return result


def global_sensitivity(trees: list[dict], blocks: dict[str, dict]):
    excluded = {}
    for name in ("all_verified_trees", "exclude_ann05_01", "exclude_ann13_01", "exclude_both"):
        remove = set()
        if name in ("exclude_ann05_01", "exclude_both"):
            remove.add("ann05.01")
        if name in ("exclude_ann13_01", "exclude_both"):
            remove.add("ann13.01")
        kept = [t for t in trees if not any(in_block(t, blocks[b]) for b in remove)]
        row = {"analysis_set": name, "n_windows": len(kept), "quartet_unit": "window-averaged 48 induced individual quartets; windows linked",
               "species_tree_inferred": "no; fixed four-species quartet concordance sensitivity only",
               "coalescent_unit_branch_lengths_estimable": "no; substitution-unit local gene trees and reticulation do not identify MSC CU branch lengths"}
        for q in Q_NAMES:
            row[q] = float(np.mean([t[q] for t in kept])) if kept else "NA"
        row["discordant_frequency"] = float(np.mean([t["discordant_frequency"] for t in kept])) if kept else "NA"
        row["mean_D_q2_minus_q1"] = float(np.mean([t["D_q2_minus_q1"] for t in kept])) if kept else "NA"
        excluded[name] = row
    return list(excluded.values())


def quartet_selection_sensitivity(trees: list[dict], blocks: dict[str, dict]) -> list[dict]:
    result=[]
    groups={"H. debilis":("",),"H. paradoxus":("par_",)}
    sets={"all_verified_trees":trees}
    for block_id,block in blocks.items():
        sets[block_id+"_interior"]=[t for t in trees if in_core(t,block)]
    for set_name,rows in sets.items():
        for fourth_species,prefixes in groups.items():
            prefix=prefixes[0]
            row={"analysis_set":set_name,"fourth_species":fourth_species,"n_windows":len(rows),
                 "quartet_definition":"q1=ann+arg; q2=ann+pet; q3=arg+pet",
                 "individual_combinations_per_tree":48,"window_is_analysis_unit":"yes; individual quartets averaged within tree"}
            for q in Q_NAMES:
                row[q]=float(np.mean([t[prefix+q] for t in rows])) if rows else "NA"
            for src,dest in (("D_q2_minus_q1","mean_D_q2_minus_q1"),("B_mean_pairwise_distance","mean_B_subs_per_site"),("mean_quartet_internal_branch","mean_internal_branch_subs_per_site")):
                row[dest]=float(np.mean([t[prefix+src] for t in rows])) if rows else "NA"
            result.append(row)
    return result


def matched_control_status(blocks: dict[str, dict], trees: list[dict]) -> list[dict]:
    # Honest feasibility check. A chromosome-matched region must be outside all
    # same-species g-marked blocks, physically equal length and within 20% of
    # focal tree count/density. Given these giant intervals and sparse windows,
    # do not manufacture pseudo-replicates if no candidates exist.
    out = []
    all_rows = read_tsv(BLOCK_FILE)
    exclusions = defaultdict(list)
    for row in all_rows:
        codes = {x.strip() for x in row["support_codes"].split(",")}
        if "g" in codes:
            for part in row["reported_regions_mbp"].split(","):
                a, b = map(float, part.strip().split("-"))
                exclusions[int(row["chromosome"])].append((round(a*1_000_000), round(b*1_000_000)))
    for block_id, block in blocks.items():
        chrom = int(block["chromosome"])
        length = sum(b-a for a,b in block["parts"])
        focal = [t for t in trees if in_block(t, block)]
        chromosome = [t for t in trees if t["chrom"] == chrom]
        eligible = []
        # Non-overlapping candidate windows of focal physical length.
        step = 1_000_000
        maxpos = max((t["end"] for t in chromosome), default=0)
        for start in range(0, maxpos-length+1, step):
            end = start + length
            if any(start < b and end > a for a,b in exclusions[chrom]):
                continue
            cand = [t for t in chromosome if start <= t["start"] and t["end"] <= end]
            if not cand or len(focal)==0:
                continue
            if abs(len(cand)-len(focal))/len(focal) <= .20:
                eligible.append((start, end, cand))
        out.append({"haploblock_id": block_id, "status": "estimable" if len(eligible) >= 20 else "not_estimable",
                    "n_eligible_same_chromosome_controls": len(eligible), "n_controls_required_for_comparison": 20,
                    "physical_length_matched": "yes" if length else "no", "tree_count_match_tolerance": "within 20%",
                    "tree_density_matched": "not separately optimized", "exclusion_mask": "all g-marked rounded haploblocks on the chromosome",
                    "reason": "too few eligible same-chromosome intervals of equal physical length and similar observed-tree count" if len(eligible)<20 else "candidate placements available; spatial dependence still requires block-aware null"})
    return out


def fig_save(fig, stem: str) -> None:
    fig.savefig(FIG / f"{stem}.png", dpi=240, bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    svg=FIG / f"{stem}.svg"
    fig.savefig(svg, bbox_inches="tight")
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines())+"\n")
    plt.close(fig)


def plot_tracks(tracks: list[dict], blocks: dict[str, dict]) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    # Figure 1: intervals and actually observed trees, not interpolated coverage.
    fig, axes = plt.subplots(2, 1, figsize=(12, 5.2), constrained_layout=True)
    for ax, block_id in zip(axes, BLOCKS):
        b = blocks[block_id]
        lo = min(a for a,_ in b["parts"])/1e6; hi=max(z for _,z in b["parts"])/1e6
        ax.axvspan(lo, hi, color="#f0c36e", alpha=.20, label="published rounded haploblock interval")
        ax.axvspan(lo, lo+EDGE_UNCERTAINTY/1e6, color="#d98e04", alpha=.32, label="±0.5-Mb endpoint uncertainty" if block_id==BLOCKS[0] else None)
        ax.axvspan(hi-EDGE_UNCERTAINTY/1e6, hi, color="#d98e04", alpha=.32)
        ax.vlines([r["midpoint_bp"]/1e6 for r in tracks if r["haploblock_id"]==block_id and r["region"]=="full_reported_interval"], .1, .9, color="#234f70", lw=.75)
        ax.set_xlim(lo-FLANK_BP/1e6, hi+FLANK_BP/1e6); ax.set_ylim(0,1); ax.set_yticks([])
        ax.set_title(f"{block_id}, chromosome {b['chromosome']}: observed tree positions (rug)", loc="left", fontsize=10)
        ax.set_xlabel("HA412-HOv2 physical coordinate (Mbp)"); ax.grid(axis="x", alpha=.18)
    fig.suptitle("Sunflower haploblock context and sparse tree coverage", weight="bold")
    fig_save(fig, "sunflower_stage1_chromosome_context")

    # Figure 2: no lines across genomic gaps.
    fig, axes=plt.subplots(2,2,figsize=(13,7),constrained_layout=True)
    for col,block_id in enumerate(BLOCKS):
        b=blocks[block_id]; subset=[r for r in tracks if r["haploblock_id"]==block_id and r["region"] in ("interior","left_flank","right_flank")]
        for row,region in enumerate(("interior","flanks")):
            ax=axes[row,col]
            vals=[r for r in subset if r["region"]==region or (region=="flanks" and r["region"].endswith("flank"))]
            for q in Q_NAMES:
                rr=sorted((r for r in vals if r[q] is not None),key=lambda r:r["midpoint_bp"])
                x=np.array([r["midpoint_bp"]/1e6 for r in rr]); y=np.array([r[q] for r in rr])
                ax.scatter(x,y,s=11,color=COLORS[q],label=q.replace("q", "q="),alpha=.85)
                for i in range(len(rr)-1):
                    if rr[i+1]["start_bp"]-rr[i]["end_bp"]<=ADJACENCY_BP:
                        ax.plot(x[i:i+2],y[i:i+2],color=COLORS[q],lw=.65,alpha=.7)
            ax.set_title(f"{block_id}: {region}",loc="left",fontsize=10); ax.set_ylim(-.03,1.03); ax.set_ylabel("induced quartet fraction")
            ax.grid(alpha=.18); ax.set_xlabel("Mbp; gaps retained")
    handles,labels=axes[0,0].get_legend_handles_labels()
    legend=axes[0,0].get_legend()
    if legend is not None:
        legend.remove()
    fig.legend(handles,labels,loc="lower center",bbox_to_anchor=(.5,-.025),frameon=False,fontsize=8,ncol=3)
    fig.suptitle("Fixed species quartet topologies across focal regions",weight="bold")
    fig_save(fig,"sunflower_stage1_quartet_topology_tracks")

    # Figure 3: branch summaries in substitution/site units.
    fig,axes=plt.subplots(4,2,figsize=(15,10),sharey="row",constrained_layout=True)
    metrics=[("A_species_pair_distance_q2_minus_q1","A-like fixed species-pair distance contrast"),
             ("B_mean_pairwise_distance","B: mean quartet-tip pairwise distance"),
             ("mean_quartet_internal_branch","Mean quartet internal branch"),
             ("D_q2_minus_q1","D(w)=q2-q1")]
    for row,(metric,label) in enumerate(metrics):
        for col,block_id in enumerate(BLOCKS):
            ax=axes[row,col]
            rr=sorted((r for r in tracks if r["haploblock_id"]==block_id and r["region"] in ("interior","left_flank","right_flank")),key=lambda r:r["midpoint_bp"])
            x=np.array([r["midpoint_bp"]/1e6 for r in rr]); y=np.array([r[metric] for r in rr],dtype=float)
            ax.scatter(x,y,s=10,color=COLORS["q2_ann_pet"],alpha=.8)
            for i in range(len(rr)-1):
                if rr[i+1]["start_bp"]-rr[i]["end_bp"]<=ADJACENCY_BP:
                    ax.plot(x[i:i+2],y[i:i+2],color=COLORS["q2_ann_pet"],lw=.6,alpha=.55)
            ax.set_title(block_id+" interior + 5-Mb flanks",loc="left",fontsize=9)
            if col==0:
                short={"A_species_pair_distance_q2_minus_q1":"A-like pair contrast\n(subs/site)","B_mean_pairwise_distance":"B mean pair distance\n(subs/site)","mean_quartet_internal_branch":"Quartet internal edge\n(subs/site)","D_q2_minus_q1":"D = q2 - q1"}[metric]
                ax.set_ylabel(short)
            else:
                ax.tick_params(labelleft=False)
            ax.set_xlabel("Mbp");ax.grid(alpha=.18)
    fig.suptitle("Arrangement-blind topology and branch-length tracks",weight="bold")
    fig_save(fig,"sunflower_stage1_branch_length_tracks")

    # Figure 4: descriptive interior/flank contrasts; no matched null claimed.
    summaries=read_tsv(RES/"region_comparison.tsv")
    fig,axes=plt.subplots(1,2,figsize=(10,4.5),constrained_layout=True)
    metrics=[("mean_D_q2_minus_q1","D = q2 - q1"),("mean_B_mean_pairwise_distance","B (subs/site)")]
    for ax,(col,label) in zip(axes,metrics):
        for i,block_id in enumerate(BLOCKS):
            group=[r for r in summaries if r["haploblock_id"]==block_id and r["region"] in ("interior","left_flank","right_flank")]
            vals={r["region"]:float(r[col]) for r in group if r[col]!="NA"}
            for j,region in enumerate(("left_flank","interior","right_flank")):
                if region in vals: ax.scatter(i*4+j,vals[region],s=38,color={"left_flank":"#777777","interior":"#D1495B","right_flank":"#999999"}[region])
            if "interior" in vals and "left_flank" in vals and "right_flank" in vals:
                ax.plot([i*4, i*4+1],[vals["left_flank"],vals["interior"]],color="#777",lw=.8)
                ax.plot([i*4+1, i*4+2],[vals["interior"],vals["right_flank"]],color="#777",lw=.8)
        ax.set_xticks([.5,4.5],BLOCKS);ax.set_ylabel(label);ax.grid(axis="y",alpha=.2)
    from matplotlib.lines import Line2D
    axes[0].legend(handles=[Line2D([0],[0],marker="o",color="w",markerfacecolor="#777777",label="left flank"),
                            Line2D([0],[0],marker="o",color="w",markerfacecolor="#D1495B",label="interior"),
                            Line2D([0],[0],marker="o",color="w",markerfacecolor="#999999",label="right flank")],
                   frameon=False,loc="lower left",fontsize=8)
    fig.suptitle("Descriptive interior and immediate-flank contrasts; matched controls unavailable",weight="bold")
    fig_save(fig,"sunflower_stage1_region_effect_sizes")

    # Figure 5: fixed-quartet sensitivity only, not species-tree inference.
    sens=read_tsv(RES/"phylogenomic_sensitivity.tsv")
    fig,ax=plt.subplots(figsize=(9,4.5),constrained_layout=True)
    xs=np.arange(len(sens));bottom=np.zeros(len(sens))
    for q in Q_NAMES:
        vals=np.array([float(r[q]) for r in sens]); ax.bar(xs,vals,bottom=bottom,label={"q1_ann_arg":"q1: ann + arg","q2_ann_pet":"q2: ann + pet","q3_arg_pet":"q3: arg + pet"}[q],color=COLORS[q]);bottom+=vals
    labels=["All", "− ann05.01", "− ann13.01", "− both"]
    ax.set_xticks(xs,labels);ax.set_ylim(0,1);ax.set_ylabel("window-averaged quartet fraction")
    ax.legend(frameon=False,ncol=3,loc="upper center",bbox_to_anchor=(.5,-.12));ax.set_title("Fixed four-species quartet sensitivity to focal-region exclusion",weight="bold")
    fig_save(fig,"sunflower_stage1_phylogenomic_sensitivity")

    # Figure 6: prespecified species-quartet choice sensitivity.
    selection=read_tsv(RES/"quartet_selection_sensitivity.tsv")
    fig,axes=plt.subplots(1,3,figsize=(13,4.5),constrained_layout=True)
    set_names=("all_verified_trees","ann05.01_interior","ann13.01_interior")
    for ax,set_name in zip(axes,set_names):
        rows=[r for r in selection if r["analysis_set"]==set_name]
        x=np.arange(2);bottom=np.zeros(2)
        for q in Q_NAMES:
            vals=np.array([float(r[q]) for r in rows]);ax.bar(x,vals,bottom=bottom,label={"q1_ann_arg":"q1 ann+arg","q2_ann_pet":"q2 ann+pet","q3_arg_pet":"q3 arg+pet"}[q],color=COLORS[q]);bottom+=vals
        ax.set_xticks(x,["DEB","PAR"]);ax.set_ylim(0,1);ax.set_title(set_name.replace("_"," "),fontsize=9);ax.grid(axis="y",alpha=.18)
    axes[0].set_ylabel("window-averaged quartet fraction")
    axes[-1].legend(frameon=False,loc="upper center",bbox_to_anchor=(1.18,1.0),fontsize=8)
    fig.suptitle("Fixed quartet result sensitivity to the fourth species",weight="bold")
    fig_save(fig,"sunflower_stage1_quartet_selection_sensitivity")


def write_report(blocks, coverage, summaries, sens, controls, tracks, all_trees):
    def fmt(v): return f"{float(v):.4f}" if isinstance(v,(float,int)) else str(v)
    c={r["haploblock_id"]:r for r in coverage}; s={(r["haploblock_id"],r["region"]):r for r in summaries}
    cctrl={r["haploblock_id"]:r for r in controls}
    lines=["# Sunflower Stage 1 spatial genealogical audit", "",
           "## Decision", "",
           "**NO-GO for claiming a defensible recombination-suppression-associated spatial genealogical effect or a sunflower MSRC success from these inputs.** The deposited Newicks support an arrangement-blind descriptive pilot, but the published regions are rounded to whole Mbp, observed windows are sparse across the very large blocks, flank/control replication is weak or unavailable, and linked windows are not independent replicates. Results below are descriptive; no significance claims are made.", "",
           "## Data, definitions, and fixed quartets", "",
           f"The analysis re-parsed {len(all_trees):,} coordinate-indexed Newicks from the local published tree table and joined called-base counts from the local window manifest. Each tree had 22 unique tips, finite nonnegative numeric branch lengths, and both fixed four-species sets. Two quartet sets were frozen: ANN (H. annuus), ARG (H. argophyllus), PET (H. petiolaris), plus either DEB (H. debilis) or PAR (H. paradoxus). All 48 observed individual combinations are scored separately within each species set and averaged to one window-level observation; they are not counted as independent loci. For both sets, q1 is the independently supported H. annuus + H. argophyllus sister pair; q2 and q3 are fixed alternatives ann-pet and arg-pet. The PAR set includes a hybrid-origin species; species histories in this radiation are reticulate and introgression is widespread, so q1 is a reference contrast, not an unambiguous bifurcating history [Owens et al. 2023](https://doi.org/10.1093/molbev/msad013).", "",
           "Todesco et al. define the reported intervals as haploblock regions; their table notes that some regions comprise rearranged pieces. Support codes `g` mean recombination suppression in genetic maps, `G` means orientation differences between genetic maps, and `H`/`h` refer to Hi-C comparisons. These are not interchangeable: the table endpoints are not fine-mapped inversion breakpoints. The focal ann05.01 carries A,g (assembly differences plus map suppression), and ann13.01 h,g (reference-relative Hi-C plus map suppression) [Todesco et al. 2020, Extended Data Table 1](https://doi.org/10.1038/s41586-020-2467-6). Published integer-Mbp endpoints were treated as ±0.5-Mb uncertainty; the primary interior excludes 0.5 Mb at every edge. Arrangement genotypes for the exact tree tips remain unavailable, so interpretation is arrangement-blind.", "",
           "## Coverage and coordinate audit", ""]
    for block in BLOCKS:
        row=c[block]; lines.append(f"- **{block}:** {row['reported_tree_overlaps']} tree windows overlap the rounded interval; {row['conservative_interior_trees']} lie wholly inside the conservative interior; {row['boundary_uncertain_trees']} touch endpoint uncertainty bands; left/right 5-Mb flanks contain {row['left_5Mb_flank_trees']}/{row['right_5Mb_flank_trees']} windows. The union of observed tree intervals spans {100*float(row['fraction_interval_spanned_by_tree_intervals']):.3f}% of the reported physical interval.")
    qc=read_tsv(RES/"tree_validation_summary.tsv")[0]
    missing=read_tsv(RES/"missing_tree_manifest_windows.tsv")
    lines += ["", f"Validation found {qc['n_unique_tree_ids']} unique tree IDs, {qc['n_tip_count_values']} distinct tip counts and {qc['n_branch_count_values']} distinct branch counts. There were {qc['n_overlapping_tree_interval_pairs']} overlapping tree-interval pairs across chromosomes. The manifest contains {qc['n_qualifying_manifest_windows']} windows with ≥2,000 called bases; {len(missing)} of these lack a deposited Newick and are listed explicitly in `missing_tree_manifest_windows.tsv`.",
              "Window spacing is highly uneven. Trees are discrete observations separated by unobserved genomic sequence; no interpolation is used. Called-base counts are retained per tree in the track table. The tree source's 3,242 records are used; no missing trees are inferred from the paper's 3,267 count.", "",
              "## Spatial topology and branch lengths", ""]
    for block in BLOCKS:
        interior=s[(block,"interior")]; left=s[(block,"left_flank")]; right=s[(block,"right_flank")]
        lines.append(f"- **{block}:** interior n={interior['n_windows']}, flanks n={left['n_windows']} left and {right['n_windows']} right. Mean q1/q2/q3 in the interior = {fmt(interior['mean_q1_ann_arg'])}/{fmt(interior['mean_q2_ann_pet'])}/{fmt(interior['mean_q3_arg_pet'])}; mean D=q2−q1 = {fmt(interior['mean_D_q2_minus_q1'])}. Mean branch summaries: A-like fixed species-pair distance contrast {fmt(interior['mean_A_species_pair_distance_q2_minus_q1'])} subs/site, B mean pairwise distance {fmt(interior['mean_B_mean_pairwise_distance'])} subs/site, and mean quartet internal branch {fmt(interior['mean_mean_quartet_internal_branch'])} subs/site.")
    lines += ["", "Here `A_species_pair_distance_q2_minus_q1` is an explicitly arrangement-blind contrast of fixed species-pair distances; it is not the arrangement-defined A(w) used in MSRC analyses. `B_mean_pairwise_distance` is the mean patristic distance across the six pairs in each fixed four-tip combination, averaged over combinations. Quartet internal branches and all distances are substitutions/site from local ML gene trees, not time or coalescent units.", "",
              "Endpoint sensitivity (0/0.5/1-Mb shrinkage) gives ann05.01 n=57/57/57 and mean D={d05}; ann13.01 n=97/94/92 and mean D={d13}. Equal-weighted nonempty 1/5/10-Mb bin means range {b05} in ann05.01 and {b13} in ann13.01; empty bins are omitted, not filled. These are sensitivity summaries, not independent replications.".format(
                  d05="/".join("%.3f"%float(r["mean_D_q2_minus_q1"]) for r in read_tsv(RES/"boundary_sensitivity.tsv") if r["haploblock_id"]=="ann05.01"),
                  d13="/".join("%.3f"%float(r["mean_D_q2_minus_q1"]) for r in read_tsv(RES/"boundary_sensitivity.tsv") if r["haploblock_id"]=="ann13.01"),
                  b05=" to ".join("%.3f"%x for x in (min(float(r["equal_bin_mean_D_q2_minus_q1"]) for r in read_tsv(RES/"bin_size_sensitivity.tsv") if r["haploblock_id"]=="ann05.01"),max(float(r["equal_bin_mean_D_q2_minus_q1"]) for r in read_tsv(RES/"bin_size_sensitivity.tsv") if r["haploblock_id"]=="ann05.01"))),
                  b13=" to ".join("%.3f"%x for x in (min(float(r["equal_bin_mean_D_q2_minus_q1"]) for r in read_tsv(RES/"bin_size_sensitivity.tsv") if r["haploblock_id"]=="ann13.01"),max(float(r["equal_bin_mean_D_q2_minus_q1"]) for r in read_tsv(RES/"bin_size_sensitivity.tsv") if r["haploblock_id"]=="ann13.01")))), "",
              "Conservative adjacent-window transitions per observed adjacent Mb are listed separately. A run absent under the 20-kb adjacency threshold does not indicate a lack of change across intervening unobserved sequence.", "",
              "Run detection treats two tree windows as adjacent only if the inter-tree gap is ≤20 kb. This conservative threshold makes almost all widely spaced tree observations isolated; an absence of adjacent runs is a lack of sampled adjacency, not evidence of genealogical continuity. Spatial autocorrelation is reported only with pair counts and is descriptive because pairs overlap and loci are linked.", "",
              "## Controls and genome-wide quartet sensitivity", ""]
    for row in controls:
        lines.append(f"- **{row['haploblock_id']} matched controls:** {row['status']}; {row['n_eligible_same_chromosome_controls']} eligible placements (20 required), matched on same chromosome and reported physical length with tree count within 20%, while excluding g-marked rounded intervals. No matched-control effect or null distribution is reported.")
    lines += ["", "The fixed q1 fraction is 0.6495 with all 3,242 trees, 0.6518 after excluding ann05.01, 0.6446 after excluding ann13.01, and 0.6469 after excluding both. These small descriptive changes do not suggest a large impact on this selected quartet summary, but no uncertainty or independent matched-removal baseline is available. The exclusion analysis recomputes a fixed four-species quartet concordance summary; it is not a full ASTRAL species-tree inference, cannot estimate coalescent-unit branches, and counts linked windows only as input records. Because suppression-associated arrangement genotypes are unavailable and reticulation violates a simple MSC interpretation, a change in these descriptive frequencies is not called phylogenomic bias. Matched-control removals could not be constructed under the predeclared physical-length and same-chromosome constraints.", "",
              "The ann13.01 q1 enrichment is similar across the two fourth-species choices: 0.806 with DEB versus 0.779 with PAR, compared with their respective genome-wide 0.650 and 0.550 baselines. Ann05.01 shows q1=0.520 (DEB set) and 0.488 (PAR set), with discordant alternatives carrying more support. The direction of ann13's q1 enrichment is therefore insensitive to this prespecified fourth-species choice, while ann05's specific discordant topology depends on which quartet set is used; neither result has matched controls.", "",
              "## Answers to the Stage 1 questions", "",
              "1. **Sustained topology shift?** No sustained alternative-topology replacement is established. Ann13.01 has q1=0.806 in its interior versus 0.650 genome-wide; ann05.01 is more discordant (q1=0.520, q2=0.295). Conservative adjacency yields only 6 and 7 adjacent interior tree pairs, respectively, and matched controls are unavailable, so these regional compositions cannot yet be attributed to recombination suppression.",
              "2. **Branch lengths inside versus outside?** The report provides descriptive interior/flank distances in substitutions/site, but does not claim a robust regional difference from sparse flanks.",
              "3. **Exceptional versus matched genomic controls?** Not estimable: eligible same-chromosome, same-length, similar-tree-count controls are insufficient.",
              "4. **Independent suppression evidence?** Both focal rows are marked `g` for genetic-map recombination suppression. ann05.01 additionally has `A`; ann13.01 has `h`. The interval coordinates describe published haploblocks, not exact inversion breakpoints.",
              "5. **Consistency across quartets?** Ann13 q1 enrichment is directionally consistent with either DEB or PAR as the fourth species; ann05 discordance is less stable in which alternative topology is favored. Both sets include all 48 sampled individual combinations averaged within each tree, not 96 independent replicates.",
              "6. **Species-tree topology or CU lengths after exclusion?** Only fixed quartet concordance sensitivity is estimable here. No species-tree topology was re-inferred and CU branch lengths are unavailable.",
              "7. **Arrangement-specific interpretation?** No. The exact 22 tree tips cannot be assigned to arrangements with available local metadata.",
              "8. **Suppression, introgression, linked selection?** These data do not distinguish mechanisms. Introgression is independently documented in this annual sunflower radiation; suppression and linked selection remain plausible co-occurring processes.",
              "9. **Strong enough plant dataset for MSRC?** No-go on current evidence. Reconsider only after exact boundaries, arrangement crosswalks, and adequately sampled matched controls become available.", "",
              "## Outputs and next step", "",
              "Coordinate and coverage audits, missing manifest windows, validation summaries, per-window quartet/distance tracks, region summaries, run lengths, transition counts, boundary/bin-size and quartet-selection sensitivity, distance-bin autocorrelation, control feasibility, and fixed-quartet exclusion sensitivity are stored in the adjacent TSVs. Six figures are provided as PNG/PDF/SVG and show only observed points; connecting segments appear only across gaps ≤20 kb. No raw inputs were downloaded or added to Git. The limiting next step is sample-to-arrangement crosswalk recovery plus finer independent interval definitions, not expansion to more haploblocks.", ""]
    (RES/"stage1_analysis_report.md").write_text("\n".join(lines))


def arrangement_audit() -> None:
    tip_names = sorted(Phylo.read(StringIO(TREE_FILE.read_text().splitlines()[0].split("\t", 1)[1]), "newick").get_terminals(), key=lambda c: c.name)
    tip_names = [tip.name for tip in tip_names]
    (RES/"arrangement_genotype_audit.md").write_text("""# Sunflower arrangement-genotype audit

The deposited local trees contain these 22 tip labels: %s. The local coordinate/called-base manifest has no individual sample metadata or haploblock genotype columns. The source study's sample list and the independent haploblock paper's population genotypes do not provide a verified crosswalk for these exact 22 tree tips at ann05.01 or ann13.01.

Therefore the Stage 1 quartet definitions aggregate tips by species prefixes only: ANN, ARG, PET, and either DEB or PAR. Both fixed species sets (ann/arg/pet/deb and ann/arg/pet/paradoxus) are analyzed as a quartet-selection sensitivity. All combinations among observed individuals in each group are scored within each tree and averaged per window. PET tips are pooled at H. petiolaris level because the six tip labels are not independently assigned to H. p. petiolaris versus H. p. fallax in the local manifest. No tip is labeled with a block arrangement, no arrangement-specific q2 is asserted, and no topology-derived arrangement label is constructed.

No public genotype matrix, diagnostic-marker assignment table, or sample-to-arrangement crosswalk was present among local inputs. No raw reads or large genotype reconstruction were attempted. If the exact 2023 sequence sample accessions can be mapped to the population-level haploblock marker data, this audit can be revised; until then all regional interpretations are arrangement-blind.
""" % ", ".join("`%s`" % name for name in tip_names))


def main() -> None:
    RES.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    blocks=parse_blocks(); trees,called,manifest=load_trees()
    qc,chrom_qc,missing,overlaps=tree_qc(trees,called,manifest)
    write_tsv(RES/"tree_validation_summary.tsv",qc)
    write_tsv(RES/"tree_validation_by_chromosome.tsv",chrom_qc,["scope","chromosome","n_trees","n_qualifying_manifest_windows","n_qualifying_manifest_windows_without_tree","n_overlapping_tree_interval_pairs","observed_coordinate_span_bp"])
    write_tsv(RES/"missing_tree_manifest_windows.tsv",missing,["tree_id","chromosome","start_bp","end_bp","called_bases","status"])
    write_tsv(RES/"overlapping_tree_intervals.tsv",overlaps,["tree_id_1","tree_id_2","chromosome","overlap_start_bp","overlap_end_bp","overlap_bp"])
    coords,coverage=calc_coverage(trees,blocks)
    write_tsv(RES/"haploblock_coordinate_audit.tsv",coords)
    write_tsv(RES/"haploblock_tree_coverage.tsv",coverage)
    tracks,summaries,runs,autocorr,regions=build_region_outputs(trees,blocks)
    write_tsv(RES/"quartet_topology_tracks.tsv",tracks)
    write_tsv(RES/"branch_length_tracks.tsv",tracks)
    write_tsv(RES/"region_comparison.tsv",summaries)
    write_tsv(RES/"spatial_run_lengths.tsv",runs,["haploblock_id","region","topology","n_adjacent_windows","start_bp","end_bp","span_bp","adjacency_rule","n_independent_loci","biological_independence_accounted_for"])
    write_tsv(RES/"spatial_autocorrelation.tsv",autocorr)
    boundary=boundary_sensitivity(trees,blocks)
    write_tsv(RES/"boundary_sensitivity.tsv",boundary)
    bin_rows=bin_sensitivity(trees,blocks)
    write_tsv(RES/"bin_size_sensitivity.tsv",bin_rows)
    transitions=transition_summary(trees,blocks)
    write_tsv(RES/"spatial_transition_counts.tsv",transitions)
    controls=matched_control_status(blocks,trees)
    write_tsv(RES/"matched_control_comparison.tsv",controls)
    sens=global_sensitivity(trees,blocks)
    write_tsv(RES/"phylogenomic_sensitivity.tsv",sens)
    selection=quartet_selection_sensitivity(trees,blocks)
    write_tsv(RES/"quartet_selection_sensitivity.tsv",selection)
    arrangement_audit()
    write_report(blocks,coverage,summaries,sens,controls,tracks,trees)
    plot_tracks(tracks,blocks)
    print(f"trees={len(trees)} focal_trees={sum(1 for t in trees if any(in_block(t,b) for b in blocks.values()))}")
    for row in coverage:
        print(f"{row['haploblock_id']}: overlap={row['reported_tree_overlaps']} core={row['conservative_interior_trees']} coverage={row['fraction_interval_spanned_by_tree_intervals']:.4%} flanks={row['left_5Mb_flank_trees']},{row['right_5Mb_flank_trees']}")


if __name__ == "__main__":
    main()
