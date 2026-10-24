#!/usr/bin/env python3
"""Independent nucleotide-site Patterson's D analysis for Anopheles 2La.

The default command is design-only. ``--pilot`` performs only the eight
preselected Stage-1R windows. ``--run-full`` is explicit and processes the
frozen 988-window Stage-2R grid in 1-Mb regional chunks.
"""
from __future__ import annotations

import argparse, csv, hashlib, importlib.metadata, json, math, os, platform, sys, time
from collections import defaultdict
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[3]
BASE = REPO / "empirical" / "anopheles_2la"
DATA = REPO / "data" / "anopheles_2la"
PROC = DATA / "processed"; META = DATA / "metadata"
OUT = BASE / "results" / "patterson_d"; FIG = BASE / "figures"
SAMPLE_MANIFEST = META / "sample_manifest.tsv"
PRIMARY = PROC / "stage1r_primary_homozygote_samples.tsv"
PILOT = PROC / "stage1r_pilot_windows.tsv"
GRID = PROC / "stage2r_window_grid.tsv"
CHROM, INV_START, INV_END, WIN = "2L", 20_524_058, 42_165_532, 50_000
SAMPLE_SET, RELEASE, BUCKET = "fontaine-2015-rebuild", "3.10", "gs://vo_agam_release_master_us_central1"
LITERATURE = {
    "fontaine2015": "https://doi.org/10.1126/science.1258524",
    "thawornwattana2018": "https://doi.org/10.1093/molbev/msy137",
}

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""): h.update(block)
    return h.hexdigest()

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f: return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n"); w.writeheader()
        for row in rows: w.writerow({k: row.get(k, "") for k in fields})

def group_ids() -> tuple[dict[str, list[str]], list[dict[str, str]]]:
    all_rows, primary = read_tsv(SAMPLE_MANIFEST), read_tsv(PRIMARY)
    groups = {
        "G_s": [r["sample_id"] for r in primary if r["species"] == "gambiae" and r["arrangement"] == "2L+a/2L+a"],
        "G_i": [r["sample_id"] for r in primary if r["species"] == "gambiae" and r["arrangement"] == "2La/2La"],
        "C_s": [r["sample_id"] for r in primary if r["species"] == "coluzzii" and r["arrangement"] == "2L+a/2L+a"],
        "C_i": [r["sample_id"] for r in primary if r["species"] == "coluzzii" and r["arrangement"] == "2La/2La"],
        "O_merus": [r["sample_id"] for r in all_rows if r["species"] == "merus"],
        "O_melas": [r["sample_id"] for r in all_rows if r["species"] == "melas"],
    }
    return groups, all_rows

def design_rows() -> list[dict[str, str]]:
    groups, _ = group_ids(); rows = []
    candidates = [("A", "G_s", "G_i", "C_s"), ("B", "G_s", "G_i", "C_i"), ("C", "C_s", "C_i", "G_s"), ("D", "C_s", "C_i", "G_i")]
    for test, p1, p2, p3 in candidates:
        for outgroup, cls in (("O_merus", "primary"), ("O_melas", "sensitivity")):
            classification = cls if test in {"A", "B"} else ("sensitivity_small_P2" if cls == "primary" else "sensitivity_small_P2_and_outgroup")
            eligibility = "eligible_as_population_test" if classification == "primary" else "eligible_sensitivity_small_group" if classification == "sensitivity_small_P2" else "eligible_sensitivity_small_group_and_outgroup"
            rows.append({"test_id": f"{test}_{outgroup.removeprefix('O_')}", "P1": p1, "P2": p2, "P3": p3, "outgroup": outgroup.removeprefix("O_"), "P1_sample_ids": ";".join(groups[p1]), "P2_sample_ids": ";".join(groups[p2]), "P3_sample_ids": ";".join(groups[p3]), "outgroup_sample_ids": ";".join(groups[outgroup]), "group_sizes": f"{len(groups[p1])},{len(groups[p2])},{len(groups[p3])},{len(groups[outgroup])}", "arrangement_states": "P1/P2/P3 are standard or inverted as encoded by sample-level Ag3 karyotypes; outgroup state not used for polarization", "assumed_rooted_topology": "((P1,P2),P3),O", "literature_justification": f"gambiae and coluzzii are sister species; merus is earliest branching in the complex; {LITERATURE['thawornwattana2018']}; {LITERATURE['fontaine2015']}", "potential_confounding": "P1/P2 are arrangement-defined within-species samples and are geographically structured; ancestral inversion polymorphism and gene flow can violate a simple tree", "classification": classification, "eligibility": eligibility})
    return rows

def design_audit() -> None:
    rows = design_rows(); write_tsv(OUT / "frozen_quartet_designs.tsv", rows, list(rows[0]))
    groups, all_rows = group_ids()
    text = f"""# Patterson's D design audit

The four candidate tests use the ordered population design `((P1,P2),P3),O`.
The primary designs use A. merus as outgroup; A. melas is retained as a
sensitivity outgroup. Groups are frozen from the existing sample-level Ag3.10
karyotype calls and Stage-1R primary homozygote table. Heterokaryotypes are
not used as ingroup P1/P2/P3 samples. Outgroup arrangement is not used to
polarize sites; only reliable outgroup nucleotide calls determine the
ancestral allele at each site.

The species-tree baseline is supported by Fontaine et al. (2015) and
Thawornwattana et al. (2018): A. gambiae and A. coluzzii are sister taxa and
A. merus is the earliest branching member of the complex. The design therefore
tests population-level asymmetry conditional on that baseline; it does not
assert that arrangement-defined samples are a species-tree quartet.

The main confounding is geographic structure: standard and inverted samples
are not geographically exchangeable within either species in this frozen
sample set. Ancestral arrangement polymorphism, incomplete lineage sorting,
and introgression can also produce asymmetry. Designs A and B are primary
exploratory tests because the gambiae groups have 8 and 12 samples and use
coluzzii as P3. Designs C and D are sensitivity tests; C_i has only 3 samples.

A. christyi is a more external outgroup in the published whole-genome
framework, but no A. christyi sample is present in the authenticated
`fontaine-2015-rebuild` sample set. It is therefore not downloaded or used in
this analysis. The existing sample metadata contain 9 merus and 4 melas
samples; their genotypes, rather than arrangement labels, determine per-site
eligibility. A fixed-outgroup estimator excludes sites where the selected
outgroup is missing or polymorphic.

The complete frozen table is `frozen_quartet_designs.tsv`. A nonzero result is
not evidence by itself for introgression or MSRC causation.
"""
    (OUT / "design_audit.md").write_text(text)

def windows(pilot: bool) -> list[dict[str, str]]:
    rows = read_tsv(PILOT if pilot else GRID)
    if pilot and len(rows) != 8: raise ValueError("Stage-1R pilot window set must contain 8 windows")
    return rows

def design_index(sample_ids: list[str], metadata: object) -> dict[str, int]:
    ids = [str(x) for x in metadata["sample_id"].tolist()]
    missing = sorted(set(sample_ids) - set(ids))
    if missing: raise ValueError(f"sample IDs absent from Ag3 sample metadata: {missing}")
    return {sid: ids.index(sid) for sid in sample_ids}

def calculate_row(gt: np.ndarray, positions: np.ndarray, window: dict[str, str], design: dict[str, str], ids: dict[str, int], min_mac: int) -> dict[str, object]:
    from patterson_d_core import calculate_window
    mask = (positions >= int(window["start"])) & (positions <= int(window["end"]))
    selected = [x for key in ("P1_sample_ids", "P2_sample_ids", "P3_sample_ids", "outgroup_sample_ids") for x in design[key].split(";") if x]
    take = [ids[x] for x in selected]
    arr = gt[mask][:, take, :]
    offsets = {}; cursor = 0
    for key in ("P1_sample_ids", "P2_sample_ids", "P3_sample_ids", "outgroup_sample_ids"):
        n = len([x for x in design[key].split(";") if x]); offsets[key] = list(range(cursor, cursor + n)); cursor += n
    result = calculate_window(arr, {"P1": offsets["P1_sample_ids"], "P2": offsets["P2_sample_ids"], "P3": offsets["P3_sample_ids"]}, offsets["outgroup_sample_ids"], min_mac=min_mac)
    row = {"window_id": window["window_id"], "chrom": window["chrom"], "start": window["start"], "end": window["end"], "region_class": window.get("region_class", window.get("position_class", "")), "test_id": design["test_id"], "P1": design["P1"], "P2": design["P2"], "P3": design["P3"], "outgroup": design["outgroup"], "filter_label": "MAC2" if min_mac else "primary", "abba": result.abba, "baba": result.baba, "abba_plus_baba": result.denominator, "D_patterson": result.d, "n_sites_examined": result.n_examined, "n_sites_pass": result.n_pass, "n_sites_informative": result.n_informative, "n_outgroup_fixed": result.outgroup_fixed_ancestral, "mean_callable_P1": result.mean_callable_chromosomes["P1"], "mean_callable_P2": result.mean_callable_chromosomes["P2"], "mean_callable_P3": result.mean_callable_chromosomes["P3"], "mean_callable_O": result.mean_callable_chromosomes["O"], "exclusion_counts": json.dumps(result.exclusions, sort_keys=True), "eligibility": "valid" if result.d is not None else "undefined_denominator"}
    return row

def load_stage2r_comparison(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    sources = {"C": BASE / "results" / "stage2r_crossed_distance_signal.tsv", "M": BASE / "results" / "stage2r_quartet_margin_signal.tsv", "D_topology": BASE / "results" / "stage2r_quartet_topology_signal.tsv"}
    maps = {key: {r["window_id"]: r for r in read_tsv(path)} for key, path in sources.items()}
    out = []
    for row in rows:
        joined = dict(row)
        joined["C"] = maps["C"].get(row["window_id"], {}).get("C", "")
        joined["M"] = maps["M"].get(row["window_id"], {}).get("mean_M", "")
        joined["D_topology"] = maps["D_topology"].get(row["window_id"], {}).get("D", "")
        out.append(joined)
    return out

def summarize(rows: list[dict[str, object]]) -> None:
    from patterson_d_core import jackknife
    valid = [r for r in rows if r["filter_label"] == "primary"]
    summary = []; jkrows = []
    for test, region in sorted({(r["test_id"], r["region_class"]) for r in valid}):
        values = [r for r in valid if r["test_id"] == test and r["region_class"] == region]
        if region == "boundary": continue
        abba, baba = sum(float(r["abba"]) for r in values), sum(float(r["baba"]) for r in values)
        d = (abba - baba) / (abba + baba) if abba + baba else None
        ds = [float(r["D_patterson"]) for r in values if r["D_patterson"] is not None and r["D_patterson"] != ""]
        summary.append({"test_id": test, "filter_label": "primary", "region": region, "n_windows": len(values), "pooled_ABBA": abba, "pooled_BABA": baba, "pooled_denominator": abba + baba, "pooled_D": d, "window_mean_D": float(np.mean(ds)) if ds else "", "window_median_D": float(np.median(ds)) if ds else "", "n_informative_windows": len(ds)})
        blocks = defaultdict(lambda: [0.0, 0.0, 0])
        for r in values:
            block = str(((int(r["start"]) - 1) // 1_000_000) * 1_000_000 + 1); blocks[block][0] += float(r["abba"]); blocks[block][1] += float(r["baba"]); blocks[block][2] += int(float(r["n_sites_informative"]))
        jk = jackknife([(a, b, n) for a, b, n in blocks.values()]); jkrows.append({"test_id": test, "filter_label": "primary", "region": region, "block_size_bp": 1_000_000, **jk})
    write_tsv(OUT / "inside_outside_summary.tsv", summary, list(summary[0]) if summary else ["test_id"])
    write_tsv(OUT / "jackknife_results.tsv", jkrows, list(jkrows[0]) if jkrows else ["test_id"])

def run_access(args: argparse.Namespace, pilot: bool) -> list[dict[str, object]]:
    import malariagen_data
    designs = {r["test_id"]: r for r in design_rows()}; selected_designs = list(designs.values())
    all_ids = sorted({x for d in designs.values() for key in ("P1_sample_ids", "P2_sample_ids", "P3_sample_ids", "outgroup_sample_ids") for x in d[key].split(";") if x})
    ag3 = malariagen_data.Ag3(url=args.url or BUCKET, check_location=False, show_progress=False)
    metadata = ag3.sample_metadata(sample_sets=SAMPLE_SET); indices = design_index(all_ids, metadata)
    windows_list = windows(pilot); results = []; output = OUT / ("pilot_results.tsv" if pilot else "window_results.tsv")
    fields = ["window_id","chrom","start","end","region_class","test_id","P1","P2","P3","outgroup","filter_label","abba","baba","abba_plus_baba","D_patterson","n_sites_examined","n_sites_pass","n_sites_informative","n_outgroup_fixed","mean_callable_P1","mean_callable_P2","mean_callable_P3","mean_callable_O","exclusion_counts","eligibility"]
    completed = set()
    if not pilot and output.exists(): completed = {(r["window_id"], r["test_id"], r["filter_label"]) for r in read_tsv(output)}
    for chunk_start in range(0, len(windows_list), max(1, args.chunk_bp // WIN)):
        chunk = windows_list[chunk_start:chunk_start + max(1, args.chunk_bp // WIN)]; region = f"{CHROM}:{chunk[0]['start']}-{chunk[-1]['end']}"; print(f"query {region} ({len(chunk)} windows)", flush=True)
        ds = ag3.snp_calls(region=region, sample_sets=SAMPLE_SET, sample_indices=[indices[x] for x in all_ids], site_mask=None, inline_array=True, chunks="native")
        positions = np.asarray(ds["variant_position"].values); gt_all = np.asarray(ds["call_genotype"].values); local_ids = {sid: i for i, sid in enumerate(all_ids)}
        chunk_rows = []
        for window in chunk:
            for design in selected_designs:
                for min_mac in (0, 2):
                    label = "MAC2" if min_mac else "primary"
                    if (window["window_id"], design["test_id"], label) in completed: continue
                    chunk_rows.append(calculate_row(gt_all, positions, window, design, local_ids, min_mac))
        results.extend(chunk_rows)
        if not pilot:
            existing = read_tsv(output) if output.exists() else []
            write_tsv(output, existing + chunk_rows, fields); (OUT / "checkpoint.json").write_text(json.dumps({"last_chunk_start": chunk_start, "completed_rows": len(existing) + len(chunk_rows), "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, indent=2) + "\n")
    if pilot: write_tsv(output, results, fields)
    if not pilot:
        rows = read_tsv(output); summarize(rows); write_tsv(OUT / "comparison_with_stage2r.tsv", load_stage2r_comparison(rows), list(load_stage2r_comparison(rows)[0]))
        make_figures(rows)
    return results

def make_figures(rows: list[dict[str, object]]) -> None:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    primary = [r for r in rows if r["filter_label"] == "primary" and r["test_id"] == "A_merus"]
    comparison = load_stage2r_comparison(primary); by = {k: {r["window_id"]: r for r in read_tsv(BASE / "results" / f"stage2r_{k}_signal.tsv")} for k in ("crossed_distance", "quartet_margin", "quartet_topology")}
    fig, axes = plt.subplots(6, 1, figsize=(11, 10), sharex=True, constrained_layout=True); tracks = [("C(w)", "C"), ("M(w)", "M"), ("D_topology(w)", "D_topology"), ("Patterson D", "D_patterson"), ("ABBA+BABA", "abba_plus_baba")]
    x = np.array([(int(r["start"]) + int(r["end"])) / 2 for r in comparison])
    axes[0].axvspan(INV_START, INV_END, color="#b39ddb", alpha=.25); axes[0].set_ylabel("2La")
    axes[0].set_yticks([])
    for ax, (label, key) in zip(axes[1:], tracks):
        vals = [float(r[key]) if r.get(key, "") not in ("", "nan") else np.nan for r in comparison]; ax.plot(x, vals, ".-", lw=.6, ms=2); ax.axhline(0, color="0.5", lw=.6); ax.set_ylabel(label); ax.axvspan(INV_START, INV_END, color="#b39ddb", alpha=.25)
    axes[-1].set_xlabel("AgamP4 2L coordinate"); FIG.mkdir(parents=True, exist_ok=True); fig.savefig(FIG / "anopheles_2la_patterson_d.pdf"); fig.savefig(FIG / "anopheles_2la_patterson_d.svg"); fig.savefig(FIG / "anopheles_2la_patterson_d.png", dpi=300); plt.close(fig)
    fig, ax = plt.subplots(figsize=(11, 4), constrained_layout=True)
    for design in sorted({r["test_id"] for r in rows}):
        subset = [r for r in rows if r["test_id"] == design and r["filter_label"] == "primary"]
        subset.sort(key=lambda r: int(r["start"])); xx = [(int(r["start"]) + int(r["end"])) / 2 for r in subset]; yy = [float(r["D_patterson"]) if r.get("D_patterson", "") not in ("", "nan") else np.nan for r in subset]
        ax.plot(xx, yy, ".-", lw=.7, ms=2, label=design)
    ax.axvspan(INV_START, INV_END, color="#b39ddb", alpha=.25); ax.axhline(0, color="0.5", lw=.6); ax.set_xlabel("AgamP4 2L coordinate"); ax.set_ylabel("Patterson's D"); ax.legend(ncol=4, fontsize=7); fig.savefig(FIG / "anopheles_2la_patterson_d_sensitivity.pdf"); fig.savefig(FIG / "anopheles_2la_patterson_d_sensitivity.svg"); fig.savefig(FIG / "anopheles_2la_patterson_d_sensitivity.png", dpi=300); plt.close(fig)

def provenance(args: argparse.Namespace) -> None:
    OUT.mkdir(parents=True, exist_ok=True); design_audit()
    payload = {"analysis": "Patterson D ABBA-BABA", "release": RELEASE, "sample_set": SAMPLE_SET, "bucket": BUCKET, "chromosome": CHROM, "inversion_interval_1based_inclusive": [INV_START, INV_END], "window_size": WIN, "stage2r_grid_sha256": sha256(GRID), "sample_manifest_sha256": sha256(SAMPLE_MANIFEST), "primary_sample_table_sha256": sha256(PRIMARY), "malariagen_data_version": importlib.metadata.version("malariagen-data") if importlib.util.find_spec("malariagen_data") else "not-imported", "host": platform.node(), "python": sys.version, "command": " ".join(sys.argv), "status": "design_frozen"}
    (OUT / "provenance.json").write_text(json.dumps(payload, indent=2) + "\n"); (OUT / "run_manifest.json").write_text(json.dumps({"status": "prepared_not_run", "pilot": "not_run", "full": "not_run", "designs": [r["test_id"] for r in design_rows()]}, indent=2) + "\n")

def main():
    p = argparse.ArgumentParser(); p.add_argument("--url", default=BUCKET); p.add_argument("--design-only", action="store_true"); p.add_argument("--pilot", action="store_true"); p.add_argument("--run-full", action="store_true"); p.add_argument("--chunk-bp", type=int, default=1_000_000); args = p.parse_args(); provenance(args)
    if args.pilot == args.run_full:
        write_tsv(OUT / "pilot_results.tsv", [{"status": "PENDING_AUTHENTICATED_ACCESS", "note": "No SNP-access package is installed in the current interpreter; no empirical pilot was run."}], ["status", "note"])
        for name in ("inside_outside_summary.tsv", "jackknife_results.tsv", "comparison_with_stage2r.tsv"):
            write_tsv(OUT / name, [{"status": "PENDING_PILOT_OR_FULL_RUN"}], ["status"])
        (OUT / "pilot_report.md").write_text("# Patterson's D pilot\n\nThe eight-window empirical pilot is pending because the current configured Python environment does not provide `malariagen_data`. No numerical pilot values are reported. The existing repository provenance records prior successful Ag3.10 metadata access, but the authenticated SNP-access environment must be restored before querying genotype arrays.\n")
        (OUT / "methods_text.md").write_text("# Patterson's D methods\n\nPatterson's D is calculated directly from unphased diploid nucleotide genotypes using population allele frequencies and a fixed-outgroup ancestral allele. ABBA and BABA are summed across sites before taking the ratio. Sites with multiallelic calls, unsuitable missingness, insufficient group callability, missing or polymorphic outgroups, and sensitivity MAC failures are retained in counts and excluded from the estimator. The primary filter uses no MAC threshold; MAC >= 2 is a sensitivity analysis matching Stage 2R.\n")
        (OUT / "results_text.md").write_text("# Patterson's D results\n\nNo empirical numerical result is available. The pilot is blocked by the missing local `malariagen_data` package and must be run in the existing authenticated Ag3.10 environment.\n")
        (OUT / "figure_caption.md").write_text("# Figure caption\n\nAligned 2L tracks will compare Stage 2R arrangement-associated contrasts with nucleotide-site Patterson's D after validated empirical output exists. Boundary-overlap windows remain visible and are excluded from the primary inside/outside comparison.\n")
        return 0
    rows = run_access(args, args.pilot); (OUT / "pilot_report.md").write_text("# Patterson's D pilot\n\nPilot results were generated from the eight coordinate-frozen Stage-1R windows. Interpretive summaries are withheld until authentication, outgroup callability, and informative-site counts have been reviewed.\n" if args.pilot else "# Patterson's D full scan\n\nThe full scan completed using the frozen Stage-2R grid. See `window_results.tsv`, `inside_outside_summary.tsv`, and `jackknife_results.tsv`.\n")
    (OUT / "methods_text.md").write_text("# Patterson's D methods\n\nPatterson's D was calculated directly from unphased diploid nucleotide genotypes using population allele frequencies and a fixed-outgroup ancestral allele. ABBA and BABA were summed across sites before taking the ratio. Sites with multiallelic calls, unsuitable missingness, insufficient group callability, missing or polymorphic outgroups, and sensitivity MAC failures were retained in counts and excluded from the estimator. The primary filter uses no MAC threshold; MAC >= 2 is a sensitivity analysis matching Stage 2R.\n")
    (OUT / "results_text.md").write_text("# Patterson's D results\n\nNumerical interpretation is deferred to the validated pilot or full output tables. Patterson's D and Stage 2R `D_topology` are separate estimands.\n")
    (OUT / "figure_caption.md").write_text("# Figure caption\n\nAligned 2L tracks compare Stage 2R arrangement-associated contrasts with nucleotide-site Patterson's D. The Patterson's D panel uses fixed candidate design A with A. merus as outgroup; alternative valid designs are reported separately. Boundary-overlap windows remain visible and are excluded from the primary inside/outside comparison.\n")
    return 0

if __name__ == "__main__": raise SystemExit(main())
