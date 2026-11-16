#!/usr/bin/env python3
"""Calculate site-frequency Patterson's D from deposited haploid calls."""
from __future__ import annotations

import csv
import gzip
import json
import math
import zipfile
from array import array
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from validate_outgroup_alignment import alignment_columns, read_fasta_member

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / "empirical/neurospora/analysis/patterson_d"
RESULTS = BASE / "results"
GENO = ROOT / "data/neurospora/raw/neurospora/geno/cras_sit_hist_tet_L10.min20MAC1.geno.gz"
WINDOWS = ROOT / "data/neurospora/raw/neurospora/trees/cras_sit_hist_tet_L10.min20MAC1.phyml_GTR_bionj_w50sites.data.tsv.gz"
TREE_WEIGHTS = ROOT / "empirical/neurospora/results/gene_tree_inventory.tsv"
BRANCHES = ROOT / "empirical/neurospora/results/branch_length_inventory.tsv"
LINKAGE_GROUPS = ["LGI", "LGII", "LGIII", "LGIV", "LGV", "LGVI", "LGVII"]
BASES = set(b"ACGT")
PRIMARY_P1 = ["CJ02", "CJ04", "CJ06", "CJ08", "CJ10", "CJ14", "CJ16", "CJ20", "CJ22", "CJ24"]
PRIMARY_P2 = ["CJ01", "CJ03", "CJ05", "CJ09", "CJ13", "CJ15", "CJ19", "CJ23"]
P3 = ["8817@999"]
OUTGROUP = "Ndis@999"
PAIR_P1 = ["CJ04", "CJ16"]
PAIR_P2 = ["CJ03", "CJ15"]
BIN_WIDTHS = [50_000, 100_000, 200_000, 250_000, 500_000]
BLOCK_WIDTHS = [250_000, 500_000, 1_000_000]


class Acc:
    def __init__(self):
        self.v = Counter()

    def add(self, result, sample_calls=None):
        self.v["source_sites"] += 1
        if sample_calls:
            for k, val in sample_calls.items():
                self.v[k] += val
        self.v[result["status"]] += 1
        if result["status"] == "polarizable":
            self.v["polarizable_sites"] += 1
            self.v["abba"] += result["abba"]
            self.v["baba"] += result["baba"]
            self.v["transversion_abba"] += result["tv_abba"]
            self.v["transversion_baba"] += result["tv_baba"]
            if result["abba"] + result["baba"] > 0:
                self.v["informative_sites"] += 1
            if result["tv_abba"] + result["tv_baba"] > 0:
                self.v["transversion_informative_sites"] += 1

    def add_arrays(self, stats, start=0, end=None, mask=None, coverage=True):
        sl = slice(start, end)
        n = len(stats["polarizable"][sl])
        if mask is None:
            mask = np.ones(n, dtype=bool)
        v = self.v
        v["source_sites"] += int(mask.sum())
        for name in ["outgroup_gap", "outgroup_ambiguous", "insufficient_calls", "monomorphic", "multiallelic", "outgroup_other"]:
            v[name] += int((stats[name][sl] & mask).sum())
        valid = stats["polarizable"][sl] & mask
        v["polarizable_sites"] += int(valid.sum())
        abba, baba = stats["abba"][sl], stats["baba"][sl]
        tv_abba, tv_baba = stats["tv_abba"][sl], stats["tv_baba"][sl]
        v["abba"] += float(abba[mask].sum())
        v["baba"] += float(baba[mask].sum())
        v["transversion_abba"] += float(tv_abba[mask].sum())
        v["transversion_baba"] += float(tv_baba[mask].sum())
        v["informative_sites"] += int(((abba + baba > 0) & valid).sum())
        v["transversion_informative_sites"] += int(((tv_abba + tv_baba > 0) & valid).sum())
        if coverage:
            for field in ["P1_called_alleles", "P1_possible_alleles", "P2_called_alleles", "P2_possible_alleles", "P3_called_alleles", "P3_possible_alleles"]:
                v[field] += int(stats[field][sl][mask].sum())

    def d(self, abba_key="abba", baba_key="baba"):
        a, b = self.v[abba_key], self.v[baba_key]
        return (a - b) / (a + b) if a + b else float("nan")

    def row(self, **keys):
        v = self.v
        return {
            **keys,
            "source_sites": v["source_sites"],
            "polarizable_sites": v["polarizable_sites"],
            "informative_sites": v["informative_sites"],
            "ABBA": v["abba"],
            "BABA": v["baba"],
            "D": self.d(),
            "transversion_informative_sites": v["transversion_informative_sites"],
            "transversion_ABBA": v["transversion_abba"],
            "transversion_BABA": v["transversion_baba"],
            "transversion_D": self.d("transversion_abba", "transversion_baba"),
            "outgroup_gap_excluded": v["outgroup_gap"],
            "outgroup_ambiguous_excluded": v["outgroup_ambiguous"],
            "insufficient_ingroup_calls_excluded": v["insufficient_calls"],
            "monomorphic_ingroup_excluded": v["monomorphic"],
            "multiallelic_ingroup_excluded": v["multiallelic"],
            "outgroup_not_ingroup_allele_excluded": v["outgroup_other"],
        }


def site_contribution(p1, p2, p3, outgroup, min_fraction=0.5):
    """Return pooled frequency contributions for one haploid nucleotide site."""
    if outgroup not in BASES:
        return {"status": "outgroup_gap" if outgroup == ord("-") else "outgroup_ambiguous", "abba": 0.0, "baba": 0.0, "tv_abba": 0.0, "tv_baba": 0.0}
    groups = [p1, p2, p3]
    called = []
    for values in groups:
        vals = [x for x in values if x in BASES]
        minimum = max(1, math.ceil(len(values) * min_fraction))
        if len(vals) < minimum:
            return {"status": "insufficient_calls", "abba": 0.0, "baba": 0.0, "tv_abba": 0.0, "tv_baba": 0.0}
        called.append(vals)
    alleles = set(x for values in called for x in values)
    if len(alleles) < 2:
        return {"status": "monomorphic", "abba": 0.0, "baba": 0.0, "tv_abba": 0.0, "tv_baba": 0.0}
    if len(alleles) > 2:
        return {"status": "multiallelic", "abba": 0.0, "baba": 0.0, "tv_abba": 0.0, "tv_baba": 0.0}
    if outgroup not in alleles:
        return {"status": "outgroup_other", "abba": 0.0, "baba": 0.0, "tv_abba": 0.0, "tv_baba": 0.0}
    derived = next(x for x in alleles if x != outgroup)
    freq = [values.count(derived) / len(values) for values in called]
    abba = (1 - freq[0]) * freq[1] * freq[2]
    baba = freq[0] * (1 - freq[1]) * freq[2]
    transition = {frozenset((ord("A"), ord("G"))), frozenset((ord("C"), ord("T")))}
    is_tv = frozenset((outgroup, derived)) not in transition
    return {
        "status": "polarizable", "abba": abba, "baba": baba,
        "tv_abba": abba if is_tv else 0.0,
        "tv_baba": baba if is_tv else 0.0,
        "ancestral": outgroup, "derived": derived,
    }


def vector_site_contributions(matrix, outgroup, groups, min_fraction=0.5):
    """Vectorized estimator; matches site_contribution's haploid frequency rules."""
    base_codes = np.array([ord("A"), ord("C"), ord("G"), ord("T")], dtype=np.uint8)
    count_arrays = []
    call_counts = []
    for indices in groups:
        g = matrix[:, indices]
        counts = np.stack([(g == b).sum(axis=1) for b in base_codes], axis=1).astype(np.int16)
        count_arrays.append(counts)
        call_counts.append(counts.sum(axis=1))
    total = count_arrays[0] + count_arrays[1] + count_arrays[2]
    present = total > 0
    allele_n = present.sum(axis=1)
    oi = np.full(len(outgroup), -1, dtype=np.int8)
    for i, b in enumerate(base_codes):
        oi[outgroup == b] = i
    outgap = outgroup == ord("-")
    outamb = (oi < 0) & ~outgap
    good_calls = np.ones(len(outgroup), dtype=bool)
    for ix, n in zip(groups, call_counts):
        minimum = max(1, math.ceil(len(ix) * min_fraction))
        good_calls &= n >= minimum
    safe = np.maximum(oi, 0)
    out_present = present[np.arange(len(outgroup)), safe]
    remaining = ~(outgap | outamb)
    insuff = remaining & ~good_calls
    remaining &= good_calls
    mono = remaining & (allele_n < 2)
    remaining &= allele_n >= 2
    multi = remaining & (allele_n > 2)
    remaining &= allele_n == 2
    outother = remaining & ~out_present
    polar = remaining & out_present

    def derived_frequency(counts, n):
        ancestral = np.take_along_axis(counts, safe[:, None], axis=1)[:, 0]
        return np.divide(n - ancestral, n, out=np.zeros(len(n), dtype=float), where=n > 0)

    p1, p2, p3 = [derived_frequency(c, n) for c, n in zip(count_arrays, call_counts)]
    abba = np.where(polar, (1 - p1) * p2 * p3, 0.0)
    baba = np.where(polar, p1 * (1 - p2) * p3, 0.0)
    derived = np.zeros(len(outgroup), dtype=np.uint8)
    for i, b in enumerate(base_codes):
        is_nonanc_present = present[:, i] & (oi != i)
        derived[is_nonanc_present] = b
    is_transition = ((outgroup == ord("A")) & (derived == ord("G"))) | ((outgroup == ord("G")) & (derived == ord("A"))) | ((outgroup == ord("C")) & (derived == ord("T"))) | ((outgroup == ord("T")) & (derived == ord("C")))
    tv_abba = np.where(polar & ~is_transition, abba, 0.0)
    tv_baba = np.where(polar & ~is_transition, baba, 0.0)
    return {
        "outgroup_gap": outgap, "outgroup_ambiguous": outamb,
        "insufficient_calls": insuff, "monomorphic": mono,
        "multiallelic": multi, "outgroup_other": outother,
        "polarizable": polar, "abba": abba, "baba": baba,
        "tv_abba": tv_abba, "tv_baba": tv_baba,
        "P1_called_alleles": call_counts[0], "P1_possible_alleles": np.full(len(outgroup), len(groups[0])),
        "P2_called_alleles": call_counts[1], "P2_possible_alleles": np.full(len(outgroup), len(groups[1])),
        "P3_called_alleles": call_counts[2], "P3_possible_alleles": np.full(len(outgroup), len(groups[2])),
    }


def read_genotypes():
    pos: dict[str, array] = defaultdict(lambda: array("I"))
    calls: dict[str, bytearray] = defaultdict(bytearray)
    with gzip.open(GENO, "rt") as f:
        samples = f.readline().strip().split()[2:]
        for line in f:
            fields = line.split()
            chrom = fields[0]
            pos[chrom].append(int(fields[1]))
            calls[chrom].extend("".join(fields[2:]).encode("ascii"))
    return samples, pos, calls


def load_windows():
    out = {}
    with gzip.open(WINDOWS, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out.setdefault(r["scaffold"], []).append((int(r["start"]), int(r["end"]), int(r["mid"]), int(r["sites"])))
    return out


def load_genealogy_bins():
    data = defaultdict(lambda: defaultdict(list))
    with TREE_WEIGHTS.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["linkage_group"] != "LGI":
                continue
            b = (int(float(r["mid"])) - 1) // 100_000
            for field, target in [("q1_species", "q1"), ("q2_other", "q2"), ("q3_hisp_A", "q3")]:
                try: data[b][target].append(float(r[field]))
                except (ValueError, KeyError): pass
    with BRANCHES.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["linkage_group"] != "LGI":
                continue
            b = (int(float(r["mid"])) - 1) // 100_000
            for field, target in [("pair1_A_hisp_minus_a_hisp", "A_pair1"), ("pair1_B_mean_six_pairwise", "B_pair1"), ("pair1_internal_branch", "internal_pair1")]:
                try:
                    x = float(r[field])
                    if math.isfinite(x): data[b][target].append(x)
                except (ValueError, KeyError): pass
    return data


def mean(values):
    values = [x for x in values if math.isfinite(x)]
    return sum(values) / len(values) if values else float("nan")


def safe_div(a, b):
    return a / b if b else float("nan")


def write_dict_tsv(path, rows):
    if not rows:
        return
    with path.open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main():
    if "frozen before D inspection" not in (BASE / "config/frozen_quartets.tsv").read_text():
        raise RuntimeError("Quartet configuration is not frozen; stop before calculating D")
    samples, positions, calls = read_genotypes()
    ix = {s: i for i, s in enumerate(samples)}
    groups = {
        "primary": ([ix[s] for s in PRIMARY_P1], [ix[s] for s in PRIMARY_P2], [ix[s] for s in P3]),
        "matched_pairs": ([ix[s] for s in PAIR_P1], [ix[s] for s in PAIR_P2], [ix[s] for s in P3]),
    }
    windows = load_windows()
    all_acc = defaultdict(Acc)
    window_acc = defaultdict(Acc)
    bin_acc = defaultdict(Acc)
    block_acc = defaultdict(Acc)
    site_mapping_check = Counter()
    with zipfile.ZipFile(ROOT / "data/neurospora/raw/outgroup_alignment.zip") as zf:
        for lg in LINKAGE_GROUPS:
            rec = read_fasta_member(zf, f"outgroup_alignment/{lg}.outgroups.fasta")
            ref, ndis = rec["2509@6"], rec["Ndis@999"]
            coords = np.frombuffer(positions[lg], dtype=np.uint32)
            matrix = np.frombuffer(calls[lg], dtype=np.uint8).reshape((-1, len(samples)))
            cols = np.frombuffer(alignment_columns(ref, positions[lg]), dtype=np.uint32)
            site_mapping_check["reference_checked"] += len(coords)
            site_mapping_check["reference_match"] += int((np.frombuffer(ref, dtype=np.uint8)[cols] == matrix[:, ix["2509@6"]]).sum())
            out = np.frombuffer(ndis, dtype=np.uint8)[cols]
            primary = vector_site_contributions(matrix, out, groups["primary"], 0.5)
            pair = vector_site_contributions(matrix, out, groups["matched_pairs"], 1.0)
            whole = np.ones(len(coords), dtype=bool)
            all_acc[("primary", "all_genome")].add_arrays(primary, mask=whole)
            all_acc[("matched_pairs", "all_genome")].add_arrays(pair, mask=whole, coverage=False)
            all_acc[("primary", lg)].add_arrays(primary, mask=whole)
            all_acc[("matched_pairs", lg)].add_arrays(pair, mask=whole, coverage=False)
            if lg == "LGI":
                inside = (coords >= 1_500_000) & (coords <= 8_000_000)
                outside = ~inside
                left = (coords >= 1_250_000) & (coords < 1_500_000)
                right = (coords > 8_000_000) & (coords <= 8_250_000)
                sens_inside = (coords >= 1_200_000) & (coords <= 8_200_000)
                for region, mask in [("LGI_inside", inside), ("LGI_outside", outside),
                                     ("LGI_left_flank_250kb", left), ("LGI_right_flank_250kb", right),
                                     ("inside_1.2_8.2Mb", sens_inside), ("outside_1.2_8.2Mb", ~sens_inside)]:
                    all_acc[("primary", region)].add_arrays(primary, mask=mask)
                    if region in ["LGI_inside", "LGI_outside", "LGI_left_flank_250kb", "LGI_right_flank_250kb"]:
                        all_acc[("matched_pairs", region)].add_arrays(pair, mask=mask, coverage=False)

            # Original variable-width 50-SNP windows.
            for start, end, mid, nsites in windows[lg]:
                lo = int(np.searchsorted(coords, start, side="left"))
                hi = int(np.searchsorted(coords, end, side="right"))
                key = (lg, start, end, mid)
                window_acc[("primary", key)].add_arrays(primary, lo, hi)
                window_acc[("matched_pairs", key)].add_arrays(pair, lo, hi, coverage=False)

            # Fixed-width physical bins, with pooled contributions and no averaging of D.
            for width in BIN_WIDTHS:
                starts = ((coords - 1) // width) * width + 1
                unique_starts = np.unique(starts)
                for start in unique_starts:
                    lo = int(np.searchsorted(starts, start, side="left"))
                    hi = int(np.searchsorted(starts, start, side="right"))
                    end = int(start + width - 1)
                    bin_acc[("primary", width, lg, int(start), end)].add_arrays(primary, lo, hi)
                    bin_acc[("matched_pairs", width, lg, int(start), end)].add_arrays(pair, lo, hi, coverage=False)

            if lg == "LGI":
                for bw in BLOCK_WIDTHS:
                    inside = (coords >= 1_500_000) & (coords <= 8_000_000)
                    left = coords < 1_500_000
                    right = coords > 8_000_000
                    block_starts = np.empty(len(coords), dtype=np.uint32)
                    block_starts[inside] = 1_500_000 + ((coords[inside] - 1_500_000) // bw) * bw
                    block_starts[left] = ((coords[left] - 1) // bw) * bw + 1
                    block_starts[right] = 8_000_001 + ((coords[right] - 8_000_001) // bw) * bw
                    for region, region_mask in [("inside", inside), ("outside", ~inside)]:
                        starts_region = np.unique(block_starts[region_mask])
                        for start in starts_region:
                            mask = region_mask & (block_starts == start)
                            block_acc[("primary", bw, region, int(start))].add_arrays(primary, mask=mask)

    if site_mapping_check["reference_checked"] != site_mapping_check["reference_match"]:
        raise RuntimeError(f"reference coordinate validation failed: {site_mapping_check}")

    # Per-source-window estimates, retaining windows with zero informative sites.
    win_rows = []
    for lg in LINKAGE_GROUPS:
        for start, end, mid, nsites in windows[lg]:
            key = (lg, start, end, mid)
            acc = window_acc[("primary", key)]
            win_rows.append(acc.row(linkage_group=lg, start=start, end=end, midpoint=mid, window_sites_assigned=nsites, inside_LGI_1p5_8Mb=int(lg == "LGI" and 1_500_000 <= mid <= 8_000_000)))
    write_dict_tsv(RESULTS / "window_patterson_d.tsv", win_rows)

    multi_rows = []
    for (design, width, lg, start, end), acc in sorted(bin_acc.items()):
        multi_rows.append(acc.row(design=design, linkage_group=lg, bin_width_bp=width, start=start, end=end, midpoint=(start + end) / 2))
    write_dict_tsv(RESULTS / "multiresolution_patterson_d.tsv", multi_rows)

    summary_rows = []
    keys = ["all_genome", *LINKAGE_GROUPS, "LGI_inside", "LGI_outside", "LGI_left_flank_250kb", "LGI_right_flank_250kb", "inside_1.2_8.2Mb", "outside_1.2_8.2Mb"]
    for design in ["primary", "matched_pairs"]:
        for key in keys:
            if (design, key) in all_acc:
                summary_rows.append(all_acc[(design, key)].row(design=design, region=key))
    for key in ["LGI_inside", "LGI_outside"]:
        acc = all_acc[("primary", key)]
        row = acc.row(design="primary_transversions_only", region=key)
        row["ABBA"], row["BABA"], row["D"] = acc.v["transversion_abba"], acc.v["transversion_baba"], acc.d("transversion_abba", "transversion_baba")
        row["transversion_ABBA"], row["transversion_BABA"], row["transversion_D"] = row["ABBA"], row["BABA"], row["D"]
        row["informative_sites"] = acc.v["transversion_informative_sites"]
        row["polarizable_sites"] = acc.v["transversion_informative_sites"]
        summary_rows.append(row)
    write_dict_tsv(RESULTS / "inside_outside_summary.tsv", summary_rows)

    # Block-jackknife uncertainty calculated separately within inside/outside strata.
    uncertainty = []
    for bw in BLOCK_WIDTHS:
        region_values = {}
        for region in ["inside", "outside"]:
            blockset = [(k, a) for k, a in block_acc.items() if k[1] == bw and k[2] == region]
            total_abba = sum(a.v["abba"] for _, a in blockset)
            total_baba = sum(a.v["baba"] for _, a in blockset)
            d_full = safe_div(total_abba - total_baba, total_abba + total_baba)
            loo = []
            for k, a in blockset:
                ab, ba = total_abba - a.v["abba"], total_baba - a.v["baba"]
                val = safe_div(ab - ba, ab + ba)
                loo.append(val)
                uncertainty.append({"block_width_bp": bw, "region": region, "block_start": k[3], "block_end": k[3] + bw - 1, "block_ABBA": a.v["abba"], "block_BABA": a.v["baba"], "block_informative_sites": a.v["informative_sites"], "leave_one_block_out_D": val, "estimate_D": d_full, "jackknife_se": "NA", "ci95_low": "NA", "ci95_high": "NA", "n_blocks": len(blockset)})
            if len(loo) >= 3:
                avg = mean(loo)
                se = math.sqrt((len(loo) - 1) / len(loo) * sum((x - avg) ** 2 for x in loo))
            else:
                se = float("nan")
            region_values[region] = (d_full, se, len(blockset))
            uncertainty.append({"block_width_bp": bw, "region": region, "block_start": "SUMMARY", "block_end": "SUMMARY", "block_ABBA": total_abba, "block_BABA": total_baba, "block_informative_sites": sum(a.v["informative_sites"] for _, a in blockset), "leave_one_block_out_D": "NA", "estimate_D": d_full, "jackknife_se": se, "ci95_low": d_full - 1.96 * se if math.isfinite(se) else "NA", "ci95_high": d_full + 1.96 * se if math.isfinite(se) else "NA", "n_blocks": len(blockset)})
        din, sein, nin = region_values["inside"]
        dout, seout, nout = region_values["outside"]
        diff = din - dout
        sed = math.sqrt(sein * sein + seout * seout) if math.isfinite(sein) and math.isfinite(seout) else float("nan")
        uncertainty.append({"block_width_bp": bw, "region": "inside_minus_outside", "block_start": "SUMMARY", "block_end": "SUMMARY", "block_ABBA": "NA", "block_BABA": "NA", "block_informative_sites": "NA", "leave_one_block_out_D": "NA", "estimate_D": diff, "jackknife_se": sed, "ci95_low": diff - 1.96 * sed if math.isfinite(sed) else "NA", "ci95_high": diff + 1.96 * sed if math.isfinite(sed) else "NA", "n_blocks": min(nin, nout)})
    write_dict_tsv(RESULTS / "block_uncertainty.tsv", uncertainty)

    # 100-kb LGI circular-shift null for the pooled inside/outside contrast.
    lg_bins = sorted((k, a) for k, a in bin_acc.items() if k[0] == "primary" and k[1] == 100_000 and k[2] == "LGI")
    max_end = max(k[4] for k, _ in lg_bins)
    nbins = (max_end + 99_999) // 100_000
    bin_lookup = {((k[3] - 1) // 100_000): a for k, a in lg_bins}
    vals = [(bin_lookup.get(i, Acc()).v["abba"], bin_lookup.get(i, Acc()).v["baba"]) for i in range(nbins)]
    observed = all_acc[("primary", "LGI_inside")].d() - all_acc[("primary", "LGI_outside")].d()
    perm_rows = []
    contrasts = []
    for shift in range(nbins):
        inside_a = inside_b = outside_a = outside_b = 0.0
        for dest in range(nbins):
            src = (dest - shift) % nbins
            a, b = vals[src]
            pos_mid = dest * 100_000 + 50_000
            if 1_500_000 <= pos_mid <= 8_000_000:
                inside_a += a; inside_b += b
            else:
                outside_a += a; outside_b += b
        din = safe_div(inside_a - inside_b, inside_a + inside_b)
        dout = safe_div(outside_a - outside_b, outside_a + outside_b)
        diff = din - dout
        contrasts.append(diff)
        perm_rows.append({"null_type": "LGI_100kb_circular_shift_inside_outside", "shift_bins": shift, "shift_bp": shift * 100_000, "D_inside": din, "D_outside": dout, "inside_minus_outside": diff, "observed_contrast": observed})
    p_spatial = sum(abs(x) + 1e-12 >= abs(observed) for x in contrasts) / len(contrasts)
    for r in perm_rows:
        r["two_sided_exact_circular_shift_p"] = p_spatial
    write_dict_tsv(RESULTS / "spatial_null_validation.tsv", perm_rows)

    # Correlate fixed physical 100-kb D bins with distinct published tree metrics.
    gene = load_genealogy_bins()
    all_bins = []
    for b in range(nbins):
        a = bin_lookup.get(b)
        d = a.d() if a else float("nan")
        g = gene.get(b, {})
        all_bins.append({
            "bin": b, "start": b * 100_000, "end": (b + 1) * 100_000,
            "D": d,
            "q1": mean(g.get("q1", [])), "q2": mean(g.get("q2", [])), "q3": mean(g.get("q3", [])),
            "A_pair1": mean(g.get("A_pair1", [])), "B_pair1": mean(g.get("B_pair1", [])), "internal_pair1": mean(g.get("internal_pair1", [])),
        })
    corr_rows = []
    for region, predicate in [("LGI_all", lambda r: True), ("LGI_inside", lambda r: 1_500_000 <= (r["start"] + r["end"]) / 2 <= 8_000_000), ("LGI_outside", lambda r: not 1_500_000 <= (r["start"] + r["end"]) / 2 <= 8_000_000)]:
        subset = [r for r in all_bins if predicate(r)]
        for metric in ["q1", "q2", "q3", "A_pair1", "B_pair1", "internal_pair1"]:
            pairs = [(r["D"], r[metric]) for r in subset if math.isfinite(r["D"]) and math.isfinite(r[metric])]
            pearson, spearman = correlations(pairs)
            corr_rows.append({"region": region, "metric": metric, "n_100kb_bins": len(pairs), "pearson_r": pearson, "spearman_rho": spearman, "interpretation": "descriptive; separate estimands, linked bins; no exact Ndis-matched tree topology"})
    write_dict_tsv(RESULTS / "spatial_correlations.tsv", corr_rows)

    # Preserve the matched coordinate table for visualization and auditing.
    matched_rows = []
    for r in all_bins:
        matched_rows.append({"linkage_group": "LGI", **r, "inside_LGI_1p5_8Mb": int(1_500_000 <= (r["start"] + r["end"]) / 2 <= 8_000_000), "genealogy_definition": "published TWISST q weights and matched CJ03/CJ04 branch metrics; Ndis is absent from Newick trees"})
    write_dict_tsv(RESULTS / "matched_nucleotide_genealogy.tsv", matched_rows)

    manifest = json.loads((RESULTS / "run_manifest.json").read_text())
    manifest.update({
        "status": "complete_primary_nucleotide_D_analysis",
        "patterson_d_calculated": True,
        "quartet_frozen_before_D": True,
        "primary_P1": PRIMARY_P1,
        "primary_P2": PRIMARY_P2,
        "P3": P3,
        "outgroup": OUTGROUP,
        "frequency_estimator": "haploid calls; group allele frequencies; minimum 50% called per multi-sample ingroup; P3 singleton required; outgroup must be unambiguous and match one of exactly two ingroup alleles",
        "mapping_reference_matches": site_mapping_check["reference_match"],
        "mapping_reference_checked": site_mapping_check["reference_checked"],
        "inside_outside_D_difference": observed,
        "circular_shift_p_100kb": p_spatial,
        "exact_Ndis_D_topology_available": False,
        "figures_pending": True,
    })
    (RESULTS / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Finished D; primary summary:")
    for k in ["all_genome", "LGI_inside", "LGI_outside", "LGI_left_flank_250kb", "LGI_right_flank_250kb"]:
        a = all_acc[("primary", k)]
        print(k, "D", a.d(), "ABBA", a.v["abba"], "BABA", a.v["baba"], "informative", a.v["informative_sites"])
    print("inside-outside", observed, "circular-shift p", p_spatial)


def ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + 1 + j) / 2
        for k in order[i:j]: out[k] = rank
        i = j
    return out


def pearson(x, y):
    if len(x) < 3: return float("nan")
    mx, my = mean(x), mean(y)
    xx = sum((v - mx) ** 2 for v in x); yy = sum((v - my) ** 2 for v in y)
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(xx * yy) if xx and yy else float("nan")


def correlations(pairs):
    if not pairs: return float("nan"), float("nan")
    x, y = zip(*pairs)
    return pearson(x, y), pearson(ranks(x), ranks(y))


if __name__ == "__main__":
    main()
