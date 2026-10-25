#!/usr/bin/env python3
"""Independent nucleotide-site Patterson's D analysis for Anopheles 2La.

The default command is design-only. ``--pilot`` performs only the eight
preselected Stage-1R windows. ``--run-full`` is explicit and processes the
frozen 988-window Stage-2R grid in 1-Mb regional chunks.
"""
from __future__ import annotations

import argparse, csv, hashlib, importlib.metadata, json, math, os, platform, resource, sys, time
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

def query_regions(window_rows: list[dict[str, str]], pilot: bool, chunk_bp: int) -> list[list[dict[str, str]]]:
    """Return exact query groups; pilot windows are always queried singly."""
    if pilot:
        return [[w] for w in window_rows]
    groups = []
    current = []
    for window in window_rows:
        if not current or (int(window["end"]) - int(current[0]["start"]) + 1) <= chunk_bp:
            current.append(window)
        else:
            groups.append(current); current = [window]
    if current: groups.append(current)
    return groups

def design_index(sample_ids: list[str], metadata: object) -> dict[str, int]:
    ids = [str(x) for x in metadata["sample_id"].tolist()]
    if len(ids) != len(set(ids)): raise ValueError("Ag3 sample metadata contains duplicate sample IDs")
    missing = sorted(set(sample_ids) - set(ids))
    if missing: raise ValueError(f"sample IDs absent from Ag3 sample metadata: {missing}")
    return {sid: ids.index(sid) for sid in sample_ids}

def dataset_sample_ids(ds: object) -> list[str]:
    """Read the returned sample coordinate rather than assuming API order."""
    if "sample_id" not in ds:
        raise ValueError("Ag3 SNP dataset did not return a sample_id coordinate")
    values = [str(x) for x in np.asarray(ds["sample_id"].values).tolist()]
    if len(values) != len(set(values)): raise ValueError("Ag3 SNP dataset returned duplicate sample IDs")
    return values

def calculate_row(gt: np.ndarray, positions: np.ndarray, window: dict[str, str], design: dict[str, str], ids: dict[str, int], min_mac: int) -> dict[str, object]:
    from patterson_d_core import calculate_window
    mask = (positions >= int(window["start"])) & (positions <= int(window["end"]))
    selected = [x for key in ("P1_sample_ids", "P2_sample_ids", "P3_sample_ids", "outgroup_sample_ids") for x in design[key].split(";") if x]
    if len(selected) != len(set(selected)): raise ValueError(f"duplicate sample ID in design {design['test_id']}")
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
    from patterson_d_core import jackknife, jackknife_difference
    valid = [r for r in rows if r["filter_label"] == "primary"]
    summary = []; jkrows = []
    for test, region in sorted({(r["test_id"], r["region_class"]) for r in valid}):
        values = [r for r in valid if r["test_id"] == test and r["region_class"] == region]
        if region == "boundary": continue
        abba, baba = sum(float(r["abba"]) for r in values), sum(float(r["baba"]) for r in values)
        d = (abba - baba) / (abba + baba) if abba + baba else None
        ds = [float(r["D_patterson"]) for r in values if r["D_patterson"] is not None and r["D_patterson"] != ""]
        summary.append({"test_id": test, "filter_label": "primary", "region": region, "n_windows": len(values), "pooled_ABBA": abba, "pooled_BABA": baba, "pooled_denominator": abba + baba, "pooled_D": d, "window_mean_D": float(np.mean(ds)) if ds else "", "window_median_D": float(np.median(ds)) if ds else "", "n_informative_windows": len(ds)})
        for block_size in (500_000, 1_000_000, 2_000_000):
            blocks = defaultdict(lambda: [0.0, 0.0, 0])
            for r in values:
                block = str(((int(r["start"]) - 1) // block_size) * block_size + 1); blocks[block][0] += float(r["abba"]); blocks[block][1] += float(r["baba"]); blocks[block][2] += int(float(r["n_sites_informative"]))
            jk = jackknife([(a, b, n) for a, b, n in blocks.values()]); jkrows.append({"test_id": test, "filter_label": "primary", "region": region, "block_size_bp": block_size, **jk})
    for test in sorted({r["test_id"] for r in valid}):
        inside_rows = [r for r in valid if r["test_id"] == test and r["region_class"] == "inside"]
        outside_rows = [r for r in valid if r["test_id"] == test and r["region_class"] == "outside"]
        for block_size in (500_000, 1_000_000, 2_000_000):
            blocks = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0, 0, 0])
            for r in inside_rows:
                block = str(((int(r["start"]) - 1) // block_size) * block_size + 1); blocks[block][0] += float(r["abba"]); blocks[block][1] += float(r["baba"]); blocks[block][4] += int(float(r["n_sites_informative"]))
            for r in outside_rows:
                block = str(((int(r["start"]) - 1) // block_size) * block_size + 1); blocks[block][2] += float(r["abba"]); blocks[block][3] += float(r["baba"]); blocks[block][5] += int(float(r["n_sites_informative"]))
            diff = jackknife_difference([tuple(v) for v in blocks.values()]); jkrows.append({"test_id": test, "filter_label": "primary", "region": "inside_minus_outside", "block_size_bp": block_size, **diff})
    write_tsv(OUT / "inside_outside_summary.tsv", summary, list(summary[0]) if summary else ["test_id"])
    write_tsv(OUT / "jackknife_results.tsv", jkrows, list(jkrows[0]) if jkrows else ["test_id"])

def run_access(args: argparse.Namespace, pilot: bool) -> list[dict[str, object]]:
    import malariagen_data
    designs = {r["test_id"]: r for r in design_rows()}; selected_designs = list(designs.values())
    all_ids = sorted({x for d in designs.values() for key in ("P1_sample_ids", "P2_sample_ids", "P3_sample_ids", "outgroup_sample_ids") for x in d[key].split(";") if x})
    ag3 = malariagen_data.Ag3(url=args.url or BUCKET, check_location=False, show_progress=False)
    metadata = ag3.sample_metadata(sample_sets=SAMPLE_SET); indices = design_index(all_ids, metadata)
    prov_path = OUT / "provenance.json"; prov = json.loads(prov_path.read_text()); prov.update({"actual_metadata_sample_count": len(metadata), "requested_sample_ids": all_ids, "requested_sample_indices": [indices[x] for x in all_ids], "sample_order_verified": True, "pilot_window_ids": [w["window_id"] for w in windows(pilot)], "pilot_query_mode": "one independent 50-kb query per window" if pilot else "contiguous bounded chunks"}); prov_path.write_text(json.dumps(prov, indent=2) + "\n")
    windows_list = windows(pilot); results = []; output = OUT / ("pilot_results.tsv" if pilot else "window_results.tsv")
    fields = ["window_id","chrom","start","end","region_class","test_id","P1","P2","P3","outgroup","filter_label","abba","baba","abba_plus_baba","D_patterson","n_sites_examined","n_sites_pass","n_sites_informative","n_outgroup_fixed","mean_callable_P1","mean_callable_P2","mean_callable_P3","mean_callable_O","exclusion_counts","eligibility","query_seconds","peak_rss_mb"]
    existing_rows = [r for r in read_tsv(output) if r.get("window_id") and r.get("test_id") and r.get("filter_label")] if output.exists() else []
    completed = {(r["window_id"], r["test_id"], r["filter_label"]) for r in existing_rows}
    for chunk in query_regions(windows_list, pilot, args.chunk_bp):
        region = f"{CHROM}:{chunk[0]['start']}-{chunk[-1]['end']}"
        started = time.perf_counter(); print(f"query {region} ({len(chunk)} window{'s' if len(chunk) != 1 else ''})", flush=True)
        ds = ag3.snp_calls(region=region, sample_sets=SAMPLE_SET, sample_indices=[indices[x] for x in all_ids], site_mask=None, inline_array=True, chunks="native")
        returned_ids = dataset_sample_ids(ds)
        if returned_ids != all_ids:
            raise ValueError("Ag3 SNP dataset sample order does not match requested sample order")
        positions = np.asarray(ds["variant_position"].values); gt_all = np.asarray(ds["call_genotype"].values); local_ids = {sid: i for i, sid in enumerate(all_ids)}
        query_seconds = time.perf_counter() - started; peak_rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 if sys.platform == "darwin" else 1024)
        chunk_rows = []
        for window in chunk:
            for design in selected_designs:
                for min_mac in (0, 2):
                    label = "MAC2" if min_mac else "primary"
                    key = (window["window_id"], design["test_id"], label)
                    if key in completed: continue
                    row = calculate_row(gt_all, positions, window, design, local_ids, min_mac)
                    row["query_seconds"], row["peak_rss_mb"] = f"{query_seconds:.6f}", f"{peak_rss_mb:.3f}"
                    chunk_rows.append(row); completed.add(key)
        results.extend(chunk_rows)
        existing = [r for r in read_tsv(output) if r.get("window_id") and r.get("test_id") and r.get("filter_label")] if output.exists() else []
        write_tsv(output, existing + chunk_rows, fields)
        update_status("pilot_running" if pilot else "full_running", completed_windows=sorted({r["window_id"] for r in existing + chunk_rows}))
        (OUT / "checkpoint.json").write_text(json.dumps({"last_region": region, "completed_rows": len(existing) + len(chunk_rows), "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, indent=2) + "\n")
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

def write_pilot_report(rows: list[dict[str, object]]) -> None:
    primary = [r for r in rows if r.get("filter_label") == "primary"]
    lines = ["# Patterson's D pilot", "", "The eight frozen Stage-1R windows were queried independently as eight 50-kb Ag3 regions (400,000 bp total). Each design has a primary no-MAC result and a MAC>=2 sensitivity result. This is a feasibility and estimator-validation pilot; it is not an inside/outside inferential analysis.", "", "## Access and mapping", "", f"- Windows: {len({r['window_id'] for r in rows})}; designs: {len({r['test_id'] for r in rows})}; rows: {len(rows)}.", "- Sample order was checked against the returned Ag3 `sample_id` coordinate for every query.", "- Outgroup polarization required complete diploid calls fixed for one allele; missing or polymorphic outgroups were excluded.", "- The primary designs A_merus and B_merus are exploratory population site-pattern contrasts under the frozen baseline, not automatic introgression tests. C/D and melas are sensitivity contrasts with smaller P2 or outgroup groups.", "", "## Primary results", "", "| test | window | D | ABBA | BABA | informative sites |", "|---|---|---:|---:|---:|---:|"]
    for r in primary:
        lines.append(f"| {r['test_id']} | {r['window_id']} | {float(r['D_patterson']):.6g} | {float(r['abba']):.6g} | {float(r['baba']):.6g} | {r['n_sites_informative']} |")
    lines += ["", "## Runtime and exclusions", "", f"- Maximum recorded resident memory: {max(float(r['peak_rss_mb']) for r in rows) / 1024:.1f} MiB.", f"- Sum of independent regional query times: {sum(float(r['query_seconds']) for r in rows[::16]):.1f} s.", "- Exclusion counts are recorded as JSON in `pilot_results.tsv`; they are design/filter diagnostics and are not evidence for introgression.", "", "Patterson's D and Stage 2R `D_topology` remain separate estimands. No nonzero value is interpreted as proof of introgression or MSRC causation."]
    (OUT / "pilot_report.md").write_text("\n".join(lines) + "\n")

def update_status(status: str, **extra: object) -> None:
    path = OUT / "run_manifest.json"; existing = json.loads(path.read_text()) if path.exists() else {"history": []}
    existing["status"] = status
    if status.startswith("pilot_"): existing["pilot"] = status.removeprefix("pilot_")
    if status.startswith("full_"): existing["full"] = status.removeprefix("full_")
    existing.update(extra); existing.setdefault("history", []).append({"status": status, "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "command": " ".join(sys.argv), "python": sys.executable})
    path.write_text(json.dumps(existing, indent=2, sort_keys=True) + "\n")

def provenance(args: argparse.Namespace) -> None:
    OUT.mkdir(parents=True, exist_ok=True); design_audit()
    payload = {"analysis": "Patterson D ABBA-BABA", "release": RELEASE, "sample_set": SAMPLE_SET, "bucket": BUCKET, "chromosome": CHROM, "inversion_interval_1based_inclusive": [INV_START, INV_END], "window_size": WIN, "stage2r_grid_sha256": sha256(GRID), "stage1r_pilot_windows_sha256": sha256(PILOT), "sample_manifest_sha256": sha256(SAMPLE_MANIFEST), "primary_sample_table_sha256": sha256(PRIMARY), "malariagen_data_version": importlib.metadata.version("malariagen-data") if importlib.util.find_spec("malariagen_data") else "not-imported", "host": platform.node(), "python": sys.version, "python_executable": sys.executable, "command": " ".join(sys.argv), "filter_settings": {"max_site_missing": 0.25, "min_callable_chromosomes": 2, "primary_mac": 0, "sensitivity_mac": 2, "outgroup_requires_complete_fixed_diploids": True}}
    if not (OUT / "provenance.json").exists(): (OUT / "provenance.json").write_text(json.dumps(payload, indent=2) + "\n")
    else:
        current = json.loads((OUT / "provenance.json").read_text()); payload["status"] = current.get("status", payload["status"]); current.update(payload); (OUT / "provenance.json").write_text(json.dumps(current, indent=2) + "\n")
    if not (OUT / "run_manifest.json").exists(): update_status("design_prepared", designs=[r["test_id"] for r in design_rows()], sample_ids=sorted({x for r in design_rows() for key in ("P1_sample_ids", "P2_sample_ids", "P3_sample_ids", "outgroup_sample_ids") for x in r[key].split(";") if x}))

def main():
    p = argparse.ArgumentParser(); p.add_argument("--url", default=BUCKET); p.add_argument("--design-only", action="store_true"); p.add_argument("--pilot", action="store_true"); p.add_argument("--run-full", action="store_true"); p.add_argument("--chunk-bp", type=int, default=1_000_000); args = p.parse_args(); provenance(args)
    if not args.pilot and not args.run_full:
        if json.loads((OUT / "run_manifest.json").read_text()).get("status") not in {"pilot_complete", "full_complete", "pilot_running", "full_running"}: update_status("design_prepared")
        return 0
    update_status("pilot_running" if args.pilot else "full_running")
    try:
        rows = run_access(args, args.pilot)
    except Exception as exc:
        update_status("pilot_failed" if args.pilot else "full_failed", error=repr(exc))
        raise
    update_status("pilot_complete" if args.pilot else "full_complete", completed_windows=sorted({r["window_id"] for r in rows}))
    prov = json.loads((OUT / "provenance.json").read_text()); prov.update({"completed_window_count": len({r["window_id"] for r in rows}), "output_sha256": sha256(OUT / ("pilot_results.tsv" if args.pilot else "window_results.tsv")), "status": "pilot_complete" if args.pilot else "full_complete"}); (OUT / "provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    if args.pilot: write_pilot_report(rows)
    else: (OUT / "pilot_report.md").write_text("# Patterson's D full scan\n\nThe full scan completed using the frozen Stage-2R grid. See `window_results.tsv`, `inside_outside_summary.tsv`, and `jackknife_results.tsv`.\n")
    (OUT / "methods_text.md").write_text("# Patterson's D methods\n\nPatterson's D was calculated directly from unphased diploid nucleotide genotypes using population allele frequencies and a fixed-outgroup ancestral allele. ABBA and BABA were summed across sites before taking the ratio. Sites with multiallelic calls, unsuitable missingness, insufficient group callability, missing or polymorphic outgroups, and sensitivity MAC failures were retained in counts and excluded from the estimator. The primary filter uses no MAC threshold; MAC >= 2 is a sensitivity analysis matching Stage 2R.\n")
    (OUT / "results_text.md").write_text("# Patterson's D results\n\nNumerical interpretation is deferred to the validated pilot or full output tables. Patterson's D and Stage 2R `D_topology` are separate estimands.\n")
    (OUT / "figure_caption.md").write_text("# Figure caption\n\nAligned 2L tracks compare Stage 2R arrangement-associated contrasts with nucleotide-site Patterson's D. The Patterson's D panel uses fixed candidate design A with A. merus as outgroup; alternative valid designs are reported separately. Boundary-overlap windows remain visible and are excluded from the primary inside/outside comparison.\n")
    return 0

if __name__ == "__main__": raise SystemExit(main())
