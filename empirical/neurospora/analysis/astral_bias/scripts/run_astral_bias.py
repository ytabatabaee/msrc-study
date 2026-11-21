#!/usr/bin/env python3
"""Prepare biologically mapped Neurospora trees and run ASTRAL-IV/CASTLES-II."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import subprocess
import time
from copy import deepcopy
from pathlib import Path

from Bio import Phylo

ROOT = Path(__file__).resolve().parents[5]
HERE = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/neurospora/raw/neurospora"
COORD = DATA / "trees/cras_sit_hist_tet_L10.min20MAC1.phyml_GTR_bionj_w50sites.data.tsv.gz"
TREES = DATA / "trees/cras_sit_hist_tet_L10.min20MAC1.phyml_GTR_bionj_w50sites.trees.gz"
MATING = ROOT / "empirical/neurospora/analysis/patterson_d/results/cj_mating_type_resolution.tsv"
BOUNDARY = (1_500_000, 8_000_000)
SUBSETS = ("all", "region_excluded", "region_inside", "LGII_VII", "LGI_outside")
SUBSET_LABELS = {
    "all": "all valid gene trees",
    "region_excluded": "all valid trees excluding LGI 1.5–8.0 Mb",
    "region_inside": "LGI 1.5–8.0 Mb only",
    "LGII_VII": "LGII–LGVII only",
    "LGI_outside": "LGI outside 1.5–8.0 Mb",
}
ASTRAL_VERSION = "v1.25.4.8"
ARCHIVE = ROOT / "data/neurospora/raw/neurospora.tar.gz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def mapping_designs(mating_path: Path, tip_names: set[str]) -> dict[str, dict[str, str]]:
    """Return species and independently supported lineage/background maps."""
    with mating_path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    resolved = {r["tree_tip"]: r["mating_type"] for r in rows
                if r["L10_membership"] == "yes" and r["mating_type"] in {"A", "a"}}
    mapped_l10 = {r["tree_tip"] for r in rows if r["L10_membership"] == "yes"}
    if resolved.keys() - tip_names or mapped_l10 - tip_names:
        raise ValueError("mating-type audit and Newick tip labels do not match")
    species = {}
    for tip in tip_names:
        if tip.startswith("CJ") or tip == "2509@6":
            species[tip] = "N_tetrasperma"
        elif tip == "Ncra@999":
            species[tip] = "N_crassa"
        elif tip == "8817@999":
            species[tip] = "N_hispaniola"
        elif tip == "1135@999":
            species[tip] = "N_sitophila"
        else:
            raise ValueError(f"unmapped tree tip: {tip}")

    shared = {
        "2509@6": "N_tetrasperma_L6",
        "Ncra@999": "N_crassa",
        "8817@999": "N_hispaniola",
        "1135@999": "N_sitophila",
    }
    background = dict(shared)
    for tip, mating_type in resolved.items():
        background[tip] = "L10_mat_A" if mating_type == "A" else "L10_mat_a"
    pair1 = dict(shared, CJ03="L10_mat_A", CJ04="L10_mat_a")
    pair2 = dict(shared, CJ15="L10_mat_A", CJ16="L10_mat_a")
    return {"species_level": species, "background_lineage": background,
            "background_pair1": pair1, "background_pair2": pair2}


def is_inside(row: dict[str, str]) -> bool:
    return row["scaffold"] == "LGI" and BOUNDARY[0] <= int(row["mid"]) <= BOUNDARY[1]


def subset_membership(row: dict[str, str]) -> dict[str, bool]:
    inside = is_inside(row)
    lgi = row["scaffold"] == "LGI"
    return {"all": True, "region_excluded": not inside, "region_inside": inside,
            "LGII_VII": not lgi, "LGI_outside": lgi and not inside}


def write_mapping(path: Path, mapping: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        for tip, species in sorted(mapping.items()):
            out.write(f"{tip}\t{species}\n")


def prepare_inputs(outdir: Path) -> dict:
    """Validate source alignment and write exact pruned ASTRAL inputs."""
    if not all(p.exists() for p in (COORD, TREES, MATING, ARCHIVE)):
        raise FileNotFoundError("required Neurospora archive or assignment audit file is missing")
    with gzip.open(COORD, "rt") as handle:
        coords = list(csv.DictReader(handle, delimiter="\t"))
    with gzip.open(TREES, "rt") as handle:
        tree_lines = [line.rstrip("\r\n") for line in handle]
    if len(coords) != 45_303 or len(tree_lines) != len(coords):
        raise ValueError(f"coordinate/tree row mismatch: {len(coords)} / {len(tree_lines)}")
    first_nwk = next(t for t in tree_lines if t not in {"NA", ""})
    first_tree = Phylo.read(io.StringIO(first_nwk), "newick")
    tip_names = {tip.name for tip in first_tree.get_terminals()}
    if len(tip_names) != 26:
        raise ValueError(f"expected 26 tree tips, got {len(tip_names)}")
    designs = mapping_designs(MATING, tip_names)
    for name, mapping in designs.items():
        if not set(mapping) <= tip_names:
            raise ValueError(f"mapping {name} contains absent tips")
        write_mapping(outdir / "inputs" / name / "taxon_mapping.tsv", mapping)

    input_counts = {design: {subset: 0 for subset in SUBSETS} for design in designs}
    missing = 0
    observed_tipsets: set[tuple[str, ...]] = set()
    paths = {d: {s: outdir / "inputs" / d / f"{s}.trees" for s in SUBSETS}
             for d in designs}
    for dpaths in paths.values():
        for path in dpaths.values():
            path.parent.mkdir(parents=True, exist_ok=True)
    writers = {p: p.open("w") for d in paths.values() for p in d.values()}
    try:
        for i, (row, nwk) in enumerate(zip(coords, tree_lines)):
            if nwk in {"NA", ""}:
                missing += 1
                continue
            tree = Phylo.read(io.StringIO(nwk), "newick")
            current_tips = {tip.name for tip in tree.get_terminals()}
            observed_tipsets.add(tuple(sorted(current_tips)))
            if current_tips != tip_names:
                raise ValueError(f"tree tip set changes at window {i + 1}")
            membership = subset_membership(row)
            for design, mapping in designs.items():
                selected = set(mapping)
                pruned = deepcopy(tree)
                for tip in list(pruned.get_terminals()):
                    if tip.name not in selected:
                        pruned.prune(tip)
                if {tip.name for tip in pruned.get_terminals()} != selected:
                    raise ValueError(f"taxon pruning failed at window {i + 1}, {design}")
                rendered = io.StringIO()
                Phylo.write(pruned, rendered, "newick", format_branch_length="%.10g")
                line = rendered.getvalue().strip()
                for subset, keep in membership.items():
                    if keep:
                        writers[paths[design][subset]].write(line + "\n")
                        input_counts[design][subset] += 1
    finally:
        for writer in writers.values():
            writer.close()
    if missing != 12_831 or len(observed_tipsets) != 1:
        raise ValueError(f"unexpected missing tree/tip-set counts: {missing}, {len(observed_tipsets)}")

    manifest = {
        "tree_archive": str(ARCHIVE), "tree_data": str(TREES), "coordinate_data": str(COORD),
        "archive_md5": "571324702709eb1c2732a3fd993edc93", "archive_sha256": sha256(ARCHIVE),
        "tree_file_sha256": sha256(TREES), "coordinate_file_sha256": sha256(COORD),
        "taxon_assignment_source": str(MATING), "taxon_assignment_sha256": sha256(MATING),
        "total_window_records": len(coords), "valid_newick_trees": 32_472,
        "NA_tree_records": missing, "tree_tip_count": len(tip_names), "unique_tree_tip_sets": len(observed_tipsets),
        "tree_branch_lengths": "BIONJ branch lengths from PhyML GTR windows; nominal substitutions/site; each locus has 50 variable sites",
        "region_coordinate_rule": "LGI midpoint, inclusive 1,500,000–8,000,000 bp, same rule as existing Stage 1; archive study axis is 2509 reference coordinates",
        "region_counts": {"all_valid": input_counts["species_level"]["all"],
                          "region_inside": input_counts["species_level"]["region_inside"],
                          "region_excluded": input_counts["species_level"]["region_excluded"],
                          "LGII_VII": input_counts["species_level"]["LGII_VII"],
                          "LGI_outside": input_counts["species_level"]["LGI_outside"]},
        "mapping_designs": {d: {"n_tips": len(m), "n_lineages": len(set(m.values())),
                                "mapping": m} for d, m in designs.items()},
        "input_tree_counts": input_counts,
        "mating_type_policy": "18 directly resolved CJ isolates included in background design; four unresolved CJ07/CJ11/CJ12/CJ18 excluded. Species-level maps all 22 CJ samples and L6 reference to N_tetrasperma.",
        "all_qc_filters": "identical archive-valid-Newick filter for each partition and sampling design; no branch-length threshold or topology filter",
    }
    (outdir / "results" / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with (outdir / "results" / "taxon_mapping_summary.tsv").open("w") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["design", "tree_tip", "mapped_taxon", "inference_scope"])
        for design, mapping in designs.items():
            scope = "species-level" if design == "species_level" else "genetic-lineage/background summary; not a species tree"
            writer.writerows([design, tip, taxon, scope] for tip, taxon in sorted(mapping.items()))
    return manifest


def run_one_astral(astral: Path, outdir: Path, design: str, subset: str,
                   length: str = "CULength", sulength_type: int = 1,
                   seed: int = 20261121) -> dict:
    input_path = outdir / "inputs" / design / f"{subset}.trees"
    mapping_path = outdir / "inputs" / design / "taxon_mapping.tsv"
    output_dir = outdir / design / subset
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"astral_{length.lower()}_type{sulength_type}.nwk"
    log_path = output_dir / f"astral_{length.lower()}_type{sulength_type}.log"
    command = [str(astral), "-R", "-u", "2", "-t", str(min(os.cpu_count() or 1, 8)),
               "--seed", str(seed), "--length", length, "--genelength", "50",
               "--sulengthtype", str(sulength_type), "-a", str(mapping_path),
               "-i", str(input_path), "-o", str(output_path), "-v", "1"]
    started = time.time()
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, check=True)
    log_path.write_text(completed.stdout)
    run = {"design": design, "subset": subset, "length_option": length,
           "sulength_type": sulength_type, "input_tree_count": sum(1 for _ in input_path.open()),
           "command": command, "output_tree": str(output_path), "log": str(log_path),
           "elapsed_seconds": time.time() - started}
    (output_dir / f"astral_{length.lower()}_type{sulength_type}.run.json").write_text(
        json.dumps(run, indent=2) + "\n")
    return run


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--astral", type=Path, default=Path("/Users/ytabatabaee/Desktop/ASTER/bin/astral4"))
    parser.add_argument("--outdir", type=Path, default=HERE)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--skip-prepare", action="store_true")
    args = parser.parse_args()
    if not args.astral.exists():
        raise FileNotFoundError(f"ASTRAL-IV executable not found: {args.astral}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    if not args.skip_prepare:
        prepare_inputs(args.outdir)
    if args.prepare_only:
        return
    runs = []
    for design in ("species_level", "background_lineage"):
        for subset in SUBSETS:
            runs.append(run_one_astral(args.astral, args.outdir, design, subset))
    for design in ("background_pair1", "background_pair2"):
        for subset in ("all", "region_excluded"):
            runs.append(run_one_astral(args.astral, args.outdir, design, subset))
    # Validate CASTLES-II's SULength export against its internal annotation.
    for design in ("species_level", "background_lineage"):
        for subset in ("all", "region_excluded"):
            runs.append(run_one_astral(args.astral, args.outdir, design, subset,
                                       length="SULength"))
    (args.outdir / "results" / "astral_run_manifest.json").write_text(
        json.dumps({"astral_version": ASTRAL_VERSION, "runs": runs}, indent=2) + "\n")
    print(f"completed {len(runs)} ASTRAL-IV/CASTLES-II analyses")


if __name__ == "__main__":
    main()
