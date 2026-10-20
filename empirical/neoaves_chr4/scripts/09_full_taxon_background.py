#!/usr/bin/env python3
"""Neoaves chromosome-4 full-taxon background-distance audit.

This stage extends the frozen Stage-6 fixed mapped-taxon distance audit.  It
does not alter the Stage-6 files and does not infer arrangement labels.  The
full-taxon quantities are explicitly branch-length summaries in ML
substitution units.
"""
from __future__ import annotations

import csv
import gc
import gzip
import hashlib
import json
import math
import statistics
import subprocess
import html
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical" / "neoaves_chr4"
OUT = EMP / "results" / "stage9_full_taxon_background"
FIG = EMP / "figures"
LOCUS = ROOT / "data/neoaves_chr4/processed/locus_table_chr4.tsv"
REGIONS = ROOT / "data/neoaves_chr4/processed/canonical_structural_regions.tsv"
MANIFEST = EMP / "results/stage4d_locus_manifest.tsv"
J48 = EMP / "data/stage4d/jarvis48_taxa.tsv"
TREES = EMP / "external/stiller2024/63430.named.gene.trees.gz"
STAGE6 = EMP / "results/stage6_genomic_tracks/window_distance_tracks.tsv"
CLADES = ("Columbea", "N61", "N62")


def read_rows(path: Path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write_rows(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def f(x, default=math.nan):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


class Tree:
    """Small Newick reader retaining all branch lengths and leaf paths."""

    def __init__(self, text: str):
        self.text = text
        self.parent = []
        self.edge = []
        self.label = []
        self.children = []
        self.i = 0
        self.root = self._parse(-1)
        self.depth = [0.0] * len(self.parent)
        for k in range(1, len(self.parent)):
            p = self.parent[k]
            self.depth[k] = self.depth[p] + self.edge[k] if math.isfinite(self.edge[k]) else math.nan
        self.leaf_ids = [k for k, ch in enumerate(self.children) if not ch and self.label[k]]
        self.labels = [self.label[k] for k in self.leaf_ids]
        self.leaf_map = {}
        self.duplicate_labels = 0
        for k in self.leaf_ids:
            if self.label[k] in self.leaf_map:
                self.duplicate_labels += 1
            else:
                self.leaf_map[self.label[k]] = k

    def _node(self, parent):
        k = len(self.parent)
        self.parent.append(parent)
        self.edge.append(0.0)
        self.label.append("")
        self.children.append([])
        return k

    def _name(self):
        j = self.i
        while self.i < len(self.text) and self.text[self.i] not in ",():;":
            self.i += 1
        return self.text[j:self.i].strip()

    def _length(self, k):
        if self.i < len(self.text) and self.text[self.i] == ":":
            self.i += 1
            j = self.i
            while self.i < len(self.text) and self.text[self.i] not in ",();":
                self.i += 1
            try:
                self.edge[k] = float(self.text[j:self.i])
            except ValueError:
                self.edge[k] = math.nan

    def _parse(self, parent):
        k = self._node(parent)
        if self.i < len(self.text) and self.text[self.i] == "(":
            self.i += 1
            self.children[k].append(self._parse(k))
            while self.i < len(self.text) and self.text[self.i] == ",":
                self.i += 1
                self.children[k].append(self._parse(k))
            if self.i >= len(self.text) or self.text[self.i] != ")":
                raise ValueError("unbalanced Newick")
            self.i += 1
            if self.i < len(self.text) and self.text[self.i] not in ":,();":
                self.label[k] = self._name()
        else:
            self.label[k] = self._name()
        self._length(k)
        return k

    def paths_for(self, taxa):
        result = {}
        for taxon in taxa:
            k = self.leaf_map.get(taxon)
            if k is None:
                continue
            path = {}
            while k >= 0:
                path[k] = self.depth[k]
                k = self.parent[k]
            result[taxon] = path
        return result

    def pair_distances(self, taxa):
        present = [x for x in taxa if x in self.leaf_map]
        paths = self.paths_for(present)
        vals = {}
        for i, a in enumerate(present):
            for b in present[i + 1:]:
                common = set(paths[a]).intersection(paths[b])
                if not common:
                    continue
                lca = max(common, key=lambda k: self.depth[k])
                v = self.depth[self.leaf_map[a]] + self.depth[self.leaf_map[b]] - 2 * self.depth[lca]
                if math.isfinite(v):
                    vals[(a, b)] = v
        return vals

    def all_pair_mean(self):
        """Mean pair distance without materializing O(n^2) pairs."""
        return self.subset_pair_mean(self.labels)

    def subset_pair_mean(self, taxa):
        """Mean distance among a selected taxon set via branch contributions."""
        selected = {self.leaf_map[x] for x in taxa if x in self.leaf_map}
        n = len(selected)
        if n < 2:
            return math.nan
        counts = [0] * len(self.parent)
        for leaf in selected:
            k = leaf
            while k >= 0:
                counts[k] += 1
                k = self.parent[k]
        total = 0.0
        for k in range(1, len(self.parent)):
            if math.isfinite(self.edge[k]):
                total += self.edge[k] * counts[k] * (n - counts[k])
        return total / (n * (n - 1) / 2)

    def metrics(self):
        finite = [x for x in self.edge if math.isfinite(x)]
        terminal = [self.edge[k] for k in self.leaf_ids if math.isfinite(self.edge[k])]
        internal = [self.edge[k] for k, ch in enumerate(self.children) if ch and k != self.root and math.isfinite(self.edge[k])]
        return {
            "n_tree_taxa": len(self.leaf_ids),
            "n_distinct_taxa": len(self.leaf_map),
            "n_duplicate_labels": self.duplicate_labels,
            "total_tree_branch_length": sum(finite) if finite else math.nan,
            "mean_terminal_branch_length": statistics.mean(terminal) if terminal else math.nan,
            "median_terminal_branch_length": statistics.median(terminal) if terminal else math.nan,
            "mean_internal_branch_length": statistics.mean(internal) if internal else math.nan,
            "n_nonfinite_branches": len(self.edge) - len(finite),
            "n_negative_branches": sum(x < 0 for x in finite),
            "n_zero_branches": sum(x == 0 for x in finite),
            "B_observed": self.all_pair_mean(),
        }


def source_records(manifest_rows, genes):
    by_name = {}
    for r in manifest_rows:
        if r["chromosome"] == "chr4" and r["original_gene_id"] in genes:
            name = (
                f"chr4_{int(f(r['start']))}_{int(f(r['end']))}_1k_start"
                f"{r['locus_start_within_window']}_fasta.treefile"
            )
            by_name[name] = r
    seen = set()
    with gzip.open(TREES, "rt") as h:
        for line in h:
            if "\t" not in line:
                continue
            name, nwk = line.rstrip("\n").split("\t", 1)
            if not name.startswith("chr4_"):
                continue
            meta = by_name.get(Path(name).name)
            if not meta or meta["original_gene_id"] in seen:
                continue
            seen.add(meta["original_gene_id"])
            yield meta, nwk


def region_membership(midpoint, regions):
    return [r["canonical_region_id"] for r in regions if f(r["canonical_start"]) <= midpoint <= f(r["canonical_end"])]


def quantiles(vals):
    vals = sorted(x for x in vals if math.isfinite(x))
    if not vals:
        return {"n": 0, "mean": math.nan, "median": math.nan, "sd": math.nan, "q05": math.nan, "q95": math.nan}
    return {"n": len(vals), "mean": statistics.mean(vals), "median": statistics.median(vals), "sd": statistics.stdev(vals) if len(vals) > 1 else 0.0, "q05": vals[max(0, int(.05 * (len(vals) - 1)))], "q95": vals[min(len(vals) - 1, int(.95 * (len(vals) - 1)))]}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    locus_rows = read_rows(LOCUS)
    genes = {r["Gene"] for r in locus_rows}
    locus_by_gene = {}
    for r in locus_rows:
        locus_by_gene.setdefault(r["Gene"], r)
    regions = read_rows(REGIONS)
    manifest_rows = read_rows(MANIFEST)
    exact_fixed = [r["stiller2024_taxon"] for r in read_rows(J48) if r["mapping_status"] == "exact"]

    # Pass 1: establish the observed full-tree taxon universe and presence.
    union = set()
    presence = Counter()
    duplicate_by_gene = {}
    n_records = 0
    for idx, (meta, text) in enumerate(source_records(manifest_rows, genes), 1):
        tree = Tree(text)
        n_records += 1
        labels = set(tree.labels)
        union.update(labels)
        for x in labels:
            presence[x] += 1
        duplicate_by_gene[meta["original_gene_id"]] = tree.duplicate_labels
        if idx % 1000 == 0:
            print(f"taxon-universe pass: {idx}/{len(genes)}", flush=True)
    expected = sorted(union)
    n_loci = n_records
    broad = sorted(x for x in expected if presence[x] / n_loci >= 0.95)
    # Fixed Stage-6 complete-case comparison uses the exact mapped labels that
    # are actually present in the local-tree universe (normally 47 of 48).
    fixed47 = sorted(x for x in exact_fixed if x in union)
    absent_fixed = sorted(x for x in exact_fixed if x not in union)

    # Pass 2: calculate candidate tracks and strict outside-region pair baselines.
    rows = {}
    baseline_sum = defaultdict(float)
    baseline_n = Counter()
    for idx, (meta, text) in enumerate(source_records(manifest_rows, genes), 1):
        gene = meta["original_gene_id"]
        tree = Tree(text)
        loc = locus_by_gene[gene]
        midpoint = f(loc["midpoint"])
        in_regions = region_membership(midpoint, regions)
        m = tree.metrics()
        d47 = tree.pair_distances(fixed47)
        broad_present = [x for x in broad if x in tree.leaf_map]
        m["B47"] = statistics.mean(d47.values()) if len(d47) == len(fixed47) * (len(fixed47) - 1) // 2 else math.nan
        m["B47_observed"] = statistics.mean(d47.values()) if d47 else math.nan
        m["n_fixed47_present"] = len({x for x in fixed47 if x in tree.leaf_map})
        m["n_broad_present"] = len({x for x in broad if x in tree.leaf_map})
        m["B_fixed"] = tree.subset_pair_mean(broad) if len(broad_present) == len(broad) else math.nan
        m["B_fixed_observed"] = tree.subset_pair_mean(broad_present)
        m["fixed47_complete"] = len(d47) == len(fixed47) * (len(fixed47) - 1) // 2
        m["fixed_broad_complete"] = len(broad_present) == len(broad)
        m["n_observed_pairs"] = len(tree.leaf_map) * (len(tree.leaf_map) - 1) // 2
        m["n_fixed47_pairs"] = len(d47)
        m["n_fixed_broad_pairs"] = len(broad_present) * (len(broad_present) - 1) // 2
        # Pair-standardized residuals use the explicit Stage-6 fixed taxon
        # universe; this avoids a quadratic broad-set materialization while
        # retaining a clearly defined, independently comparable pair universe.
        for pair, dist in d47.items():
            if not in_regions:
                baseline_sum[pair] += dist
                baseline_n[pair] += 1
        rows[gene] = {**loc, **m, "region_ids": ";".join(in_regions)}
        if idx % 1000 == 0:
            print(f"distance pass: {idx}/{n_loci}", flush=True)
            gc.collect()

    pair_mu = {pair: baseline_sum[pair] / baseline_n[pair] for pair in baseline_sum if baseline_n[pair] >= 20}

    # Pass 3: pair-standardized residuals, summarized per locus to avoid an
    # enormous pair-by-window table while retaining counts and the weighting.
    residual_rows = []
    for idx, (meta, text) in enumerate(source_records(manifest_rows, genes), 1):
        gene = meta["original_gene_id"]
        tree = Tree(text)
        dist = tree.pair_distances(fixed47)
        residual = [v - pair_mu[p] for p, v in dist.items() if p in pair_mu]
        rel = [(v - pair_mu[p]) / pair_mu[p] for p, v in dist.items() if p in pair_mu and pair_mu[p] > 1e-12]
        rows[gene]["B_residual"] = statistics.mean(residual) if residual else math.nan
        rows[gene]["B_residual_relative"] = statistics.mean(rel) if rel else math.nan
        rows[gene]["n_residual_pairs"] = len(residual)
        residual_rows.append({"Gene": gene, "start": rows[gene]["start"], "end": rows[gene]["end"], "midpoint": rows[gene]["midpoint"], "n_pairs": len(residual), "mean_residual": rows[gene]["B_residual"], "mean_relative_residual": rows[gene]["B_residual_relative"], "baseline_windows_min": min((baseline_n[p] for p in dist if p in pair_mu), default=0)})
        if idx % 1000 == 0:
            print(f"residual pass: {idx}/{n_loci}", flush=True)
            gc.collect()

    ordered = sorted(rows.values(), key=lambda x: f(x["midpoint"]))
    estimator_fields = ["Gene", "start", "end", "midpoint", "n_tree_taxa", "n_distinct_taxa", "n_duplicate_labels", "n_observed_pairs", "n_fixed47_present", "n_fixed47_pairs", "fixed47_complete", "B47", "B47_observed", "n_broad_present", "n_fixed_broad_pairs", "fixed_broad_complete", "B_fixed", "B_fixed_observed", "B_observed", "B_residual", "B_residual_relative", "n_residual_pairs", "total_tree_branch_length", "mean_terminal_branch_length", "median_terminal_branch_length", "mean_internal_branch_length", "n_nonfinite_branches", "n_negative_branches", "n_zero_branches", "region_ids"]
    write_rows(OUT / "background_distance_estimators.tsv", ordered, estimator_fields)
    write_rows(OUT / "pair_standardized_residuals.tsv", residual_rows, list(residual_rows[0]))

    cov_rows = []
    for r in ordered:
        cov_rows.append({"Gene": r["Gene"], "start": r["start"], "end": r["end"], "midpoint": r["midpoint"], "n_tree_taxa": r["n_tree_taxa"], "n_expected_taxa": len(expected), "n_missing_expected": len(expected) - int(r["n_distinct_taxa"]), "n_extra_taxa": 0, "n_duplicate_labels": r["n_duplicate_labels"], "n_fixed47_present": r["n_fixed47_present"], "fixed47_complete": r["fixed47_complete"], "n_broad_present": r["n_broad_present"], "fixed_broad_complete": r["fixed_broad_complete"], "n_observed_pairs": r["n_observed_pairs"], "in_canonical_region": bool(r["region_ids"])})
    write_rows(OUT / "per_window_taxon_coverage.tsv", cov_rows, list(cov_rows[0]))
    presence_rows = [{"taxon": x, "presence_count": presence[x], "presence_fraction": presence[x] / n_loci, "in_stage6_exact_set": x in exact_fixed, "in_fixed_broad_set": x in broad} for x in expected]
    write_rows(OUT / "species_presence_summary.tsv", presence_rows, list(presence_rows[0]))

    def summary_row(name, vals):
        q = quantiles(vals)
        return {"metric": name, **q}
    summary = [summary_row("n_tree_taxa", [f(r["n_tree_taxa"]) for r in ordered]), summary_row("B47_complete", [f(r["B47"]) for r in ordered]), summary_row("B_observed", [f(r["B_observed"]) for r in ordered]), summary_row("B_fixed_broad_complete", [f(r["B_fixed"]) for r in ordered]), summary_row("B_residual", [f(r["B_residual"]) for r in ordered])]
    for label, subset in (("coverage_inside_canonical", [r for r in ordered if r["region_ids"]]), ("coverage_outside_canonical", [r for r in ordered if not r["region_ids"]])):
        summary.append(summary_row(label, [f(r["n_tree_taxa"]) for r in subset]))
    write_rows(OUT / "taxon_coverage_summary.tsv", summary, list(summary[0]))

    # Stage-6 regression: compare B47 with the frozen complete-case track.
    old = {r["Gene"]: r for r in read_rows(STAGE6)}
    diffs = [abs(f(r["B47"]) - f(old[r["Gene"]].get("B_fixed_mapped"))) for r in ordered if r["Gene"] in old and math.isfinite(f(r["B47"])) and math.isfinite(f(old[r["Gene"]].get("B_fixed_mapped")))]
    comparisons = []
    names = ["B47", "B_observed", "B_fixed", "B_residual"]
    for a in names:
        for b in names:
            pairs = [(f(r[a]), f(r[b])) for r in ordered if math.isfinite(f(r[a])) and math.isfinite(f(r[b]))]
            if pairs:
                xa, xb = zip(*pairs)
                ma, mb = statistics.mean(xa), statistics.mean(xb)
                cov = sum((x - ma) * (y - mb) for x, y in pairs)
                sa = math.sqrt(sum((x - ma) ** 2 for x in xa)); sb = math.sqrt(sum((y - mb) ** 2 for y in xb))
                comparisons.append({"estimator_a": a, "estimator_b": b, "n_windows": len(pairs), "pearson": cov / (sa * sb) if sa and sb else math.nan, "mean_abs_difference": statistics.mean(abs(x-y) for x,y in pairs), "mean_a": ma, "mean_b": mb})
    write_rows(OUT / "estimator_comparison.tsv", comparisons, list(comparisons[0]))

    sens = []
    for region in regions:
        s, e = f(region["canonical_start"]), f(region["canonical_end"])
        for est in names:
            inside = [f(r[est]) for r in ordered if s <= f(r["midpoint"]) <= e and math.isfinite(f(r[est]))]
            outside = [f(r[est]) for r in ordered if not (s <= f(r["midpoint"]) <= e) and math.isfinite(f(r[est]))]
            sens.append({"region": region["canonical_region_id"], "estimator": est, "inside_n": len(inside), "outside_n": len(outside), "inside_mean": statistics.mean(inside) if inside else math.nan, "outside_mean": statistics.mean(outside) if outside else math.nan, "inside_minus_outside": (statistics.mean(inside)-statistics.mean(outside)) if inside and outside else math.nan})
    write_rows(OUT / "structural_region_sensitivity.tsv", sens, list(sens[0]))

    # Fixed broad selection is based only on taxon presence, with a predeclared
    # 95% locus coverage threshold. It is never selected by signal strength.
    write_rows(OUT / "fixed_taxon_selection.tsv", [{"selection": "Stage6_exact_observed", "threshold": "exact mapping and observed in universe", "n_taxa": len(fixed47), "n_complete_windows": sum(bool(r["fixed47_complete"]) for r in ordered), "excluded_taxa": ";".join(absent_fixed)}, {"selection": "fixed_broad_95pct", "threshold": "presence in >=95% of chr4 local trees", "n_taxa": len(broad), "n_complete_windows": sum(bool(r["fixed_broad_complete"]) for r in ordered), "excluded_taxa": ""}], ["selection", "threshold", "n_taxa", "n_complete_windows", "excluded_taxa"])

    # Frozen support values, retained as separate focal designs.
    support = []
    for r in locus_rows:
        q = [f(r[k]) for k in ("q1", "q2", "q3")]
        support.append({"Gene": r["Gene"], "clade": r["clade"], "start": r["start"], "end": r["end"], "midpoint": r["midpoint"], "q1": q[0], "q2": q[1], "q3": q[2], "D2": q[1]-q[0], "D3": q[2]-q[0]})
    preferred = "B_observed"
    # Prefer the broad fixed set only if it retains at least half the loci;
    # this rule is based on coverage, not structural signal.
    if sum(bool(r["fixed_broad_complete"]) for r in ordered) >= 0.5 * len(ordered):
        preferred = "B_fixed"

    make_figure(support, ordered, regions, preferred)

    audit_lines = [
        "# Full-taxon data and estimator audit",
        "",
        f"The chromosome-4 local-tree universe contains {len(expected)} distinct observed taxon labels across {n_loci} unique loci. This is the empirical expected universe used here; no absent 363-species labels were invented.",
        f"The source file contains up to {len(expected)} observed taxa per local tree. The complete expected-universe count is {sum(int(r['n_distinct_taxa']) == len(expected) for r in ordered)} loci.",
        f"Stage-6 exact mapping has {len(exact_fixed)} exact entries; {len(fixed47)} occur in the local-tree universe and {len(absent_fixed)} are absent. Stage-6 complete-case regression used the observed fixed set.",
        f"The fixed broad set is defined before inspecting signal as taxa present in >=95% of chr4 local trees: {len(broad)} taxa.",
        "B_observed is the mean patristic distance over all distinct observed local-tree taxa. B_fixed requires complete representation of the fixed broad set. B_residual is a pair-standardized residual centered by pair-specific strict-outside-region medians; it is not an absolute distance.",
        "All distances are ML substitution units. They are not calibrated divergence times or coalescence times. Missing taxa are retained as missing observations and are never filled by interpolation or label guessing.",
        f"Stage-6 B regression maximum absolute difference over matched complete rows: {max(diffs) if diffs else math.nan:.6g}.",
    ]
    (OUT / "report.md").write_text("\n".join(audit_lines) + "\n")
    (OUT / "methods_text.md").write_text("""# Methods

We reused the named chromosome-4 ML local trees and the frozen QQS locus table. The expected full-taxon universe was defined as the union of exact leaf labels observed across the unique chromosome-4 loci. All-observed mean pairwise distance was calculated directly from branch-length contributions, while fixed-set estimates used patristic paths for the selected taxa. A fixed broad set was prespecified by >=95% locus presence. Pair-standardized residuals used pair-specific medians estimated from windows outside all canonical structural regions, requiring at least 20 observations per pair. QQS support tracks remain separate for Columbea, N61, and N62. No arrangement labels were inferred and no Stage-6–8 output was overwritten.
""")
    (OUT / "results_text.md").write_text(f"""# Results

The full local-tree label universe contains {len(expected)} taxa across {n_loci} loci. Local taxon composition varies, so B_observed and fixed-set estimators answer different questions. The Stage-6 fixed-set values reproduce with maximum absolute difference {max(diffs) if diffs else math.nan:.6g}. The preferred background track for the new figure is {preferred}, chosen by the predeclared coverage rule rather than structural-region signal.

The 47-taxon arrangement-blind distance restriction is retained as a regression/reference analysis. The six-species Stage-8 arrangement-adjusted analysis remains separate and is not merged into this full-taxon background track.
""")
    (OUT / "validation_report.md").write_text(f"""# Validation report

- Unique locus records analyzed: {n_loci}
- Observed taxon universe: {len(expected)}
- Full-universe complete local trees: {sum(int(r['n_distinct_taxa']) == len(expected) for r in ordered)}
- Stage-6 fixed observed taxon count: {len(fixed47)}
- Stage-6 regression rows: {len(diffs)}
- Stage-6 regression maximum absolute difference: {max(diffs) if diffs else math.nan:.6g}
- Fixed broad set (>=95% presence): {len(broad)} taxa
- Preferred new background estimator: {preferred}

This validates numerical reproduction and coverage accounting; it does not establish a causal rearrangement effect.
""")
    (OUT / "figure_caption.md").write_text("""# Figure caption

Chromosome-4 structural regions, source-defined focal quartet support, topology contrasts, and a full-taxon background divergence track. Columbea, N61, and N62 remain separate quartet designs. The background track is a mean ML patristic distance in substitution units, with incomplete taxon coverage retained as missing observations. Structural intervals are canonical regions, not nucleotide-resolved causal breakpoints. No arrangement-informed A(w) is shown because the full taxon set lacks independent binary arrangement labels.
""")
    manifest = {"stage": "Neoaves chr4 Stage 9 full-taxon background divergence", "script": str(Path(__file__).relative_to(ROOT)), "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in [LOCUS, REGIONS, MANIFEST, J48, TREES, STAGE6]}, "n_loci": n_loci, "n_expected_observed_taxa": len(expected), "n_stage6_exact_mapped": len(exact_fixed), "n_stage6_exact_observed": len(fixed47), "n_fixed_broad_95pct": len(broad), "preferred_estimator": preferred, "stage6_regression_max_abs_difference": max(diffs) if diffs else None, "random_seed": None}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def plot_with_gaps(ax, x, y, color, lw=0.55):
    xx, yy = [], []
    for i, (a, b) in enumerate(zip(x, y)):
        if i and math.isfinite(x[i-1]) and a - x[i-1] > 0.25:
            if xx:
                ax.plot(xx, yy, color=color, lw=lw)
            xx, yy = [], []
        if math.isfinite(b):
            xx.append(a); yy.append(b)
    if xx:
        ax.plot(xx, yy, color=color, lw=lw)


def make_figure(support, ordered, regions, preferred):
    """Write a lightweight vector figure without depending on a GUI backend."""
    width, left, right = 1600, 130, 35
    panel_h, gap, top = 120, 18, 55
    n_panels = 8
    height = top + n_panels * panel_h + (n_panels - 1) * gap + 60
    xmax = max(f(r["end"]) for r in ordered) / 1e6
    colors = {"q1": "#2673B8", "q2": "#D65F5F", "q3": "#4A9B68", "d2": "#D65F5F", "d3": "#4A9B68", "b": "#6A3D9A"}
    usable = width - left - right

    def esc(v):
        return html.escape(str(v), quote=True)

    def xmap(x):
        return left + usable * x / xmax

    def ymap(v, lo, hi, y0):
        return y0 + panel_h - panel_h * (v - lo) / (hi - lo)

    def path_for(xs, ys, lo, hi, y0):
        segments, current = [], []
        for i, (x, y) in enumerate(zip(xs, ys)):
            if not math.isfinite(x) or not math.isfinite(y) or (i and xs[i] - xs[i-1] > 0.25):
                if current: segments.append(" ".join(current)); current = []
                continue
            current.append(("M" if not current else "L") + f" {xmap(x):.2f},{ymap(y, lo, hi, y0):.2f}")
        if current: segments.append(" ".join(current))
        return segments

    def region_marks(y0):
        out = []
        for rr in regions:
            s, e = f(rr["canonical_start"]) / 1e6, f(rr["canonical_end"]) / 1e6
            out.append(f'<rect x="{xmap(s):.2f}" y="{y0:.2f}" width="{max(1,xmap(e)-xmap(s)):.2f}" height="{panel_h}" fill="#bdbdbd" opacity="0.16"/>')
            out.append(f'<line x1="{xmap(s):.2f}" x2="{xmap(s):.2f}" y1="{y0:.2f}" y2="{y0+panel_h:.2f}" stroke="#777" stroke-dasharray="4,4" stroke-width="0.7"/>')
            out.append(f'<line x1="{xmap(e):.2f}" x2="{xmap(e):.2f}" y1="{y0:.2f}" y2="{y0+panel_h:.2f}" stroke="#777" stroke-dasharray="4,4" stroke-width="0.7"/>')
        return "".join(out)

    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>', '<style>text{font-family:Arial,sans-serif;fill:#222} .axis{stroke:#333;stroke-width:1} .grid{stroke:#ddd;stroke-width:.6}</style>']
    svg.append(f'<text x="{left}" y="28" font-size="18" font-weight="bold">Neoaves chromosome 4: full-taxon background divergence and focal quartet designs</text>')
    y0 = top
    # structural-region annotation
    svg.append(region_marks(y0)); svg.append(f'<text x="{left-10}" y="{y0+panel_h/2+4}" text-anchor="end" font-size="12">regions</text>')
    for rr in regions:
        s, e = f(rr["canonical_start"]) / 1e6, f(rr["canonical_end"]) / 1e6
        if e-s > 1:
            svg.append(f'<text x="{(xmap(s)+xmap(e))/2:.2f}" y="{y0+panel_h/2+4:.2f}" text-anchor="middle" font-size="9">{esc(rr["canonical_region_id"])}</text>')
    y0 += panel_h + gap
    for cl in CLADES:
        z = sorted([r for r in support if r["clade"] == cl], key=lambda r: f(r["midpoint"]))
        xs = [f(r["midpoint"]) / 1e6 for r in z]
        svg.append(region_marks(y0))
        for key, col in (("q1", colors["q1"]), ("q2", colors["q2"]), ("q3", colors["q3"])):
            for path in path_for(xs, [r[key] for r in z], 0, 1, y0):
                svg.append(f'<path d="{path}" fill="none" stroke="{col}" stroke-width="0.7"/>')
        svg += [f'<line x1="{left}" x2="{width-right}" y1="{y0+panel_h}" y2="{y0+panel_h}" class="axis"/>', f'<text x="{left-10}" y="{y0+panel_h/2+4}" text-anchor="end" font-size="11">{esc(cl)} q support</text>']
        if cl == "Columbea":
            svg += [f'<text x="{width-right-250}" y="{y0+18}" font-size="11" fill="{colors["q1"]}">q1 species-history</text>', f'<text x="{width-right-130}" y="{y0+18}" font-size="11" fill="{colors["q2"]}">q2 / q3 alternatives</text>']
        y0 += panel_h + gap
        svg.append(region_marks(y0))
        svg.append(f'<line x1="{left}" x2="{width-right}" y1="{ymap(0,-1,1,y0):.2f}" y2="{ymap(0,-1,1,y0):.2f}" class="grid"/>')
        for key, col in (("D2", colors["d2"]), ("D3", colors["d3"])):
            for path in path_for(xs, [r[key] for r in z], -1, 1, y0):
                svg.append(f'<path d="{path}" fill="none" stroke="{col}" stroke-width="0.7"/>')
        svg += [f'<text x="{left-10}" y="{y0+panel_h/2+4}" text-anchor="end" font-size="11">{esc(cl)} D2 / D3</text>', f'<text x="{width-right-140}" y="{y0+18}" font-size="10" fill="#555">D2=q2−q1; D3=q3−q1</text>']
        y0 += panel_h + gap
    svg.append(region_marks(y0))
    bvals = [f(r[preferred]) for r in ordered if math.isfinite(f(r[preferred]))]
    blo, bhi = (min(bvals), max(bvals)) if bvals else (0, 1)
    if bhi <= blo: bhi = blo + 1
    xs = [f(r["midpoint"]) / 1e6 for r in ordered]
    for path in path_for(xs, [f(r[preferred]) for r in ordered], blo, bhi, y0):
        svg.append(f'<path d="{path}" fill="none" stroke="{colors["b"]}" stroke-width="1.0"/>')
    svg += [f'<text x="12" y="{y0+panel_h/2-2}" text-anchor="start" font-size="10">{esc(preferred)}</text>', f'<text x="12" y="{y0+panel_h/2+12}" text-anchor="start" font-size="9">substitution units</text>', f'<text x="{width-right-330}" y="{y0+18}" font-size="10" fill="#555">ML patristic distance; incomplete taxa retained</text>']
    for tick in range(0, int(xmax)+1, 10):
        x = xmap(tick)
        svg += [f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{y0+panel_h}" y2="{y0+panel_h+5}" class="axis"/>', f'<text x="{x:.2f}" y="{y0+panel_h+22}" text-anchor="middle" font-size="10">{tick}</text>']
    svg.append(f'<text x="{left+usable/2}" y="{height-10}" text-anchor="middle" font-size="12">GalGal6 chromosome 4 position (Mb)</text>')
    svg.append('</svg>')
    svg_path = OUT / "neoaves_chr4_stage9_full_taxon_multitrack.svg"
    svg_path.write_text("".join(svg))
    subprocess.run(["rsvg-convert", "-f", "pdf", "-o", str(FIG / "neoaves_chr4_stage9_full_taxon_multitrack.pdf"), str(svg_path)], check=True)
    subprocess.run(["rsvg-convert", "-f", "png", "-d", "300", "-p", "300", "-o", str(FIG / "neoaves_chr4_stage9_full_taxon_multitrack.png"), str(svg_path)], check=True)


if __name__ == "__main__":
    main()
