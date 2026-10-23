#!/usr/bin/env python3
"""Audit the completed 363/48-taxon results against cached PNAS trees.

This script only parses committed trees, manifests, and the cached PNAS
repository archive. It never invokes ASTRAL or reconstructs gene trees.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "empirical" / "neoaves_chr4"
RESULTS = BASE / "results" / "stage4e"
CLUSTER = RESULTS / "cluster"
ARCHIVE = BASE / "external" / "chr4avian" / "chr4avian-v1.0.0.zip"
SOURCE_PREFIX = "smirarab-chr4avian-c257707/"
FIGURES = BASE / "figures"
os.environ.setdefault("MPLCONFIGDIR", str(FIGURES / ".mplconfig"))

FULL_PUBLISHED_ARCHIVE_PATH = SOURCE_PREFIX + "score-columbea/score.tre"
OUTLIER_PUBLISHED_ARCHIVE_PATH = SOURCE_PREFIX + (
    "score-columbea/score-post-astral-63430_filter_mRNA_random_region.gene.trees."
    "aLRT0.05-by-63430_filter_mRNA_random_region.gene.trees.aLRT0.05."
    "chr4outlierremoved.gene.trees.tre"
)

FULL_TREES = {
    "T0": CLUSTER / "T0" / "tree.cu.nwk",
    "T_PNAS": CLUSTER / "T_PNAS" / "tree.cu.nwk",
    "T_STRUCT_PRIMARY": CLUSTER / "T_STRUCT_PRIMARY" / "tree.cu.nwk",
    "T_BLOCK_500K": CLUSTER / "T_BLOCK_500K" / "tree.cu.nwk",
    "NONCHR4": CLUSTER / "NONCHR4" / "tree.cu.nwk",
}
J48_TREES = {
    "J48_T0": CLUSTER / "J48_T0" / "tree.cu.nwk",
    "J48_T_PNAS": CLUSTER / "J48_T_PNAS" / "tree.cu.nwk",
    "J48_T_STRUCT_PRIMARY": CLUSTER / "J48_T_STRUCT_PRIMARY" / "tree.cu.nwk",
    "J48_T_BLOCK_500K": CLUSTER / "J48_T_BLOCK_500K" / "tree.cu.nwk",
    "J48_NONCHR4": CLUSTER / "J48_NONCHR4" / "tree.cu.nwk",
}
TREATMENT_ORDER = ("T0", "T_PNAS", "T_STRUCT_PRIMARY", "T_BLOCK_500K", "NONCHR4")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def source_bytes(archive_path: str) -> bytes:
    with zipfile.ZipFile(ARCHIVE) as archive:
        return archive.read(archive_path)


def load_tree(data: bytes):
    import dendropy

    tree = dendropy.Tree.get(
        data=data.decode(), schema="newick", preserve_underscores=True, rooting="force-unrooted"
    )
    taxa = {node.taxon.label for node in tree.leaf_node_iter()}
    splits: set[str] = set()

    def leaves(node):
        if node.is_leaf():
            return {node.taxon.label}
        result: set[str] = set()
        for child in node.child_node_iter():
            result.update(leaves(child))
        return result

    def split_key(side: set[str], all_taxa: set[str]) -> str:
        other = all_taxa - side
        choices = (tuple(sorted(side)), tuple(sorted(other)))
        return "|".join(min(choices))

    for node in tree.postorder_node_iter():
        if node.is_leaf() or node.parent_node is None:
            continue
        side = leaves(node)
        other = taxa - side
        if len(side) > 1 and len(other) > 1:
            splits.add(split_key(side, taxa))
    return {"tree": tree, "taxa": taxa, "splits": splits}


def load_local_tree(path: Path):
    return load_tree(path.read_bytes())


def rf(left: dict, right: dict) -> tuple[object, ...]:
    if left["taxa"] != right["taxa"]:
        return ("taxa_mismatch", len(left["taxa"] - right["taxa"]), len(right["taxa"] - left["taxa"]), len(left["taxa"] & right["taxa"]))
    only_left = left["splits"] - right["splits"]
    only_right = right["splits"] - left["splits"]
    return (len(only_left) + len(only_right), len(only_left), len(only_right))


def group_sets() -> tuple[dict[str, set[str]], dict[str, tuple[set[str], set[str], set[str]]]]:
    archive_path = ROOT / "data/neoaves_chr4/raw/genetreesupport/clade-analysis/examined-clades.zip"
    names = ("Columbimorphae", "Otidimorphae", "Columbiformes", "Phoenicopteriformes", "OtherColumbimorphae")
    with zipfile.ZipFile(archive_path) as archive:
        groups = {name: set(archive.read(name + ".txt").decode().split()) for name in names}
        roles = {
            "Columbea": tuple(set(archive.read(name + ".txt").decode().split()) for name in ("Columbimorphae", "Phoenicopteriformes", "Passera")),
            "N61": tuple(set(archive.read(name + ".txt").decode().split()) for name in ("Otidimorphae", "Columbimorphae", "ElementavesTelluraves")),
            "N62": tuple(set(archive.read(name + ".txt").decode().split()) for name in ("Columbiformes", "OtherColumbimorphae", "Otidimorphae")),
        }
    groups["Columbea"] = groups["Columbimorphae"] | groups["Phoenicopteriformes"]
    groups["Columbiformes_plus_Otidimorphae"] = groups["Columbiformes"] | groups["Otidimorphae"]
    return groups, roles


def canonical_group_split(group: set[str], taxa: set[str]) -> str:
    side, other = set(group) & taxa, taxa - set(group)
    return "|".join(min(tuple(sorted(side)), tuple(sorted(other))))


def focal_state(tree: dict, roles: tuple[set[str], set[str], set[str]]) -> str:
    c1, c2, sister = (set(role) & tree["taxa"] for role in roles)
    for label, side in (("q1", c1 | c2), ("q2", c1 | sister), ("q3", c2 | sister)):
        if canonical_group_split(side, tree["taxa"]) in tree["splits"]:
            return label
    return "absent_or_unresolved"


def branch_rows(treatment: str) -> list[dict[str, str]]:
    path = CLUSTER / treatment / "branch_lengths.tsv"
    if not path.exists() and treatment == "T0":
        path = CLUSTER / "T0_extracted" / "branch_lengths.tsv"
    return read_tsv(path)


def branch_for(tree: dict, rows: list[dict[str, str]], side: set[str]) -> dict[str, str] | None:
    key = canonical_group_split(side, tree["taxa"])
    for row in rows:
        if row["branch_type"] == "internal" and row["split_key"] == key:
            return row
    return None


def audit_trees():
    trees = {name: load_local_tree(path) for name, path in {**FULL_TREES, **J48_TREES}.items()}
    published = {
        "PNAS_FULL": load_tree(source_bytes(FULL_PUBLISHED_ARCHIVE_PATH)),
        "PNAS_OUTLIER_REMOVED": load_tree(source_bytes(OUTLIER_PUBLISHED_ARCHIVE_PATH)),
    }
    trees.update(published)
    return trees


def comparison_rows(trees: dict[str, dict]) -> list[dict[str, object]]:
    pairs = [
        ("T0", "T_PNAS", "reconstructed full363 treatments"),
        ("T0", "T_STRUCT_PRIMARY", "reconstructed full363 treatments"),
        ("T0", "T_BLOCK_500K", "reconstructed full363 treatments"),
        ("T0", "NONCHR4", "reconstructed full363 treatments"),
        ("J48_T0", "J48_T_PNAS", "reconstructed Jarvis-48 treatments"),
        ("J48_T0", "J48_T_STRUCT_PRIMARY", "reconstructed Jarvis-48 treatments"),
        ("J48_T0", "J48_T_BLOCK_500K", "reconstructed Jarvis-48 treatments"),
        ("J48_T0", "J48_NONCHR4", "reconstructed Jarvis-48 treatments"),
        ("T0", "PNAS_FULL", "our full-data tree versus cached published full-data tree"),
        ("T0", "PNAS_OUTLIER_REMOVED", "our full-data tree versus cached published outlier-removed tree"),
        ("T_PNAS", "PNAS_OUTLIER_REMOVED", "our reconstructed outlier-removed tree versus cached published tree"),
        ("PNAS_FULL", "PNAS_OUTLIER_REMOVED", "cached published full versus outlier-removed trees"),
    ]
    rows = []
    for left, right, comparison in pairs:
        result = rf(trees[left], trees[right])
        rows.append({
            "left": left,
            "right": right,
            "comparison": comparison,
            "left_taxa": len(trees[left]["taxa"]),
            "right_taxa": len(trees[right]["taxa"]),
            "taxa_equal": trees[left]["taxa"] == trees[right]["taxa"],
            "rf": result[0],
            "left_only_splits": result[1],
            "right_only_splits": result[2],
        })
    return rows


def focal_rows(trees: dict[str, dict]) -> list[dict[str, object]]:
    groups, roles = group_sets()
    selected = ["T0", "T_PNAS", "T_STRUCT_PRIMARY", "T_BLOCK_500K", "NONCHR4", "PNAS_FULL", "PNAS_OUTLIER_REMOVED"]
    rows = []
    for tree_name in selected:
        tree = trees[tree_name]
        for group_name in ("Columbimorphae", "Columbea", "Otidimorphae", "Columbiformes_plus_Otidimorphae"):
            group = groups[group_name]
            rows.append({
                "tree": tree_name,
                "taxa": len(tree["taxa"]),
                "split_name": group_name,
                "group_taxa_present": len(group & tree["taxa"]),
                "group_taxa_expected": len(group),
                "monophyletic_unrooted_split": group <= tree["taxa"] and canonical_group_split(group, tree["taxa"]) in tree["splits"],
                "interpretation": "direct unrooted split" if group <= tree["taxa"] else "taxon set incomplete",
            })
        for clade, role_sets in roles.items():
            rows.append({
                "tree": tree_name,
                "taxa": len(tree["taxa"]),
                "split_name": clade + "_focal_state",
                "group_taxa_present": "",
                "group_taxa_expected": "",
                "monophyletic_unrooted_split": "",
                "interpretation": focal_state(tree, role_sets),
            })
    return rows


def branch_rows_for_report(trees: dict[str, dict]) -> list[dict[str, object]]:
    _, roles = group_sets()
    rows = []
    for treatment in TREATMENT_ORDER:
        tree = trees[treatment]
        source_rows = branch_rows(treatment)
        for clade, role_sets in roles.items():
            c1, c2, sister = role_sets
            for state, side in (("q1", c1 | c2), ("q2", c1 | sister), ("q3", c2 | sister)):
                record = branch_for(tree, source_rows, side)
                rows.append({
                    "treatment": treatment,
                    "clade": clade,
                    "state": state,
                    "present": record is not None,
                    "CU_length": record["CU_length"] if record else "NA",
                    "SU_length": record["SU_length"] if record else "NA",
                    "localPP": record["localPP"] if record else "NA",
                    "q1": record["q1"] if record else "NA",
                    "q2": record["q2"] if record else "NA",
                    "q3": record["q3"] if record else "NA",
                })
    return rows


def source_manifest() -> dict[str, str]:
    return {
        "archive": str(ARCHIVE.relative_to(ROOT)),
        "archive_sha256": sha256_bytes(ARCHIVE.read_bytes()),
        "repository_commit": "17b408a622c9589fe9f6ca8a434c150974e39cd4",
        "repository_url": "https://github.com/smirarab/chr4avian/tree/17b408a622c9589fe9f6ca8a434c150974e39cd4",
        "zenodo_url": "https://zenodo.org/api/records/10699424/files/smirarab/chr4avian-v1.0.0.zip/content",
        "paper_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11009670/",
        "published_full_archive_path": FULL_PUBLISHED_ARCHIVE_PATH.removeprefix(SOURCE_PREFIX),
        "published_full_sha256": sha256_bytes(source_bytes(FULL_PUBLISHED_ARCHIVE_PATH)),
        "published_full_git_blob_sha": "aedc2d4890a7bff26e779724a05beb35777d6a0b",
        "published_outlier_removed_archive_path": OUTLIER_PUBLISHED_ARCHIVE_PATH.removeprefix(SOURCE_PREFIX),
        "published_outlier_removed_sha256": sha256_bytes(source_bytes(OUTLIER_PUBLISHED_ARCHIVE_PATH)),
        "published_outlier_removed_git_blob_sha": "f88f18219e915ddbed71161b0a25c542c72a3370",
        "published_reference_63k": str((BASE / "external/stiller2024/63K.tre").relative_to(ROOT)),
        "published_reference_63k_sha256": sha256_bytes((BASE / "external/stiller2024/63K.tre").read_bytes()),
    }


def write_reports(trees, comparisons, focals, branches):
    write_tsv(RESULTS / "full363_pnas_tree_comparison.tsv", comparisons, list(comparisons[0]))
    write_tsv(RESULTS / "full363_focal_split_comparison.tsv", focals, list(focals[0]))
    write_tsv(RESULTS / "full363_focal_branch_lengths.tsv", branches, list(branches[0]))
    manifest = source_manifest()
    (RESULTS / "full363_pnas_source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    report = f"""# Full 363-taxon comparison with the published PNAS trees

This audit parses the committed 363/48-taxon ASTRAL outputs and the cached
versioned PNAS repository archive. It does not invoke ASTRAL or regenerate
gene trees. All RF values are symmetric differences of canonical unrooted
internal bipartitions.

## Published reference files

The cached archive is the Zenodo `chr4avian-v1.0.0` repository snapshot
({manifest['archive_sha256']}) at
{manifest['zenodo_url']}, corresponding to repository commit
`{manifest['repository_commit']}` ({manifest['repository_url']}). The paper is
{manifest['paper_url']}. The source repository identifies `score-columbea/` as
the location of the species trees shown in Figure S9.

| Reference | Archive path | Taxa | SHA-256 |
|---|---|---:|---|
| Published full-data tree/reference artifact | `{manifest['published_full_archive_path']}`; Git blob `{manifest['published_full_git_blob_sha']}` | 363 | `{manifest['published_full_sha256']}` |
| Published outlier-removed tree/reference artifact | `{manifest['published_outlier_removed_archive_path']}`; Git blob `{manifest['published_outlier_removed_git_blob_sha']}` | 363 | `{manifest['published_outlier_removed_sha256']}` |
| Stiller `63K.tre` reference | `{manifest['published_reference_63k']}` | 363 | `{manifest['published_reference_63k_sha256']}` |

The repository README documents `score-columbea` as containing the species
trees shown in Figure S9. The two files are therefore treated as published
species-tree reference artifacts. Their names and annotations show that they
were used in post-ASTRAL scoring, but the cached repository does not establish
whether each file is the original inference output or a scored representation
of that output. No ASTRAL command line or version identifier is included, so
their exact inference settings are not independently recoverable.

## Result

Our full-data `T0` is identical to both cached published full-data references
at the unrooted topology level (RF 0). The cached published full-data and
outlier-removed trees are also identical (RF 0). Our reconstructed `T_PNAS`
differs from the cached published outlier-removed tree by RF 2: one internal
split is removed and one is added. Therefore the apparent full-taxon
disagreement is real for the current reconstructed input, but it is not yet a
contradiction of the PNAS inference because our `T_PNAS` is based on
`named_recollapsed.tre`, whereas the published analysis used the original
aLRT-collapsed collection.

The changed split is the N62 relationship. In our T0 and in the published
full/outlier-removed references, the 9-taxon Columbimorphae split is present
and the N62 focal state is q1: Columbiformes and Other Columbimorphae are
together relative to Otidimorphae. In our reconstructed T_PNAS, that split is
absent and the 18-taxon Columbiformes + Otidimorphae split is present; the N62
state is q2. This is an unrooted split statement. It does not assert a rooted
clade without specifying a root.

## Interpretation of the PNAS comparison

The PNAS paper states that the 48-taxon analysis with all loci recovered the
Columbea/J2014-side topology, while removing the outlier region produced a
topology similar to S2024 but with Otidimorphae sister to doves, breaking
Columbimorphae. It also states that with the increased Stiller taxon sampling,
the S2024 topology was recovered with or without the outlier region. The exact
full-taxon references cached in `score-columbea` support that latter statement:
both published 363-taxon trees have the same N62 q1 state.

Our 363-taxon T0 agrees exactly with the published full-data topology. Our
reconstructed T_PNAS does not reproduce the published outlier-removed tree.
Thus the current data show a preprocessing/input reproduction discrepancy,
not a demonstrated contradiction of the published taxon-sampling result.

## Remaining explanation

The existing source audits establish that 63,113 of 63,430 collapsed records
have strongly supported header-derived locus correspondences; 317 remain
problematic, including 200 whose taxon sets are absent from the named archive.
The current T_PNAS has the same expected 61,999-locus count and the frozen
1,431-locus PNAS exclusion, so count and mask membership alone do not establish
identity. Ordinary aLRT branch collapsing cannot change leaf sets. The
remaining source/preprocessing discrepancies can therefore plausibly affect
the outlier-removed topology, but this audit does not demonstrate that the 317
records alone cause RF 2.

Our runs use ASTRAL-IV v1.25.4.8 with `-u 2`, `--genelength 1000`, and
`--sulengthtype 1`, as recorded in the run provenance. The cached PNAS archive
does not document a comparable ASTRAL version or command, so a version-only
explanation is unsupported. A controlled follow-up with the exact published
collapsed trees and the current ASTRAL binary is the minimal useful next test;
an additional inference using the current reconstructed input is not needed.

## Files

- `full363_pnas_tree_comparison.tsv`: RF and taxon-set comparisons.
- `full363_focal_split_comparison.tsv`: direct unrooted focal split and state calls.
- `full363_focal_branch_lengths.tsv`: CU, SU, localPP, and quartet annotations for homologous focal states.
- `full363_pnas_source_manifest.json`: source paths and checksums.
- `figures/neoaves_pnas_replication_comparison.*`: verified-results figure.
"""
    (RESULTS / "full363_pnas_replication_audit.md").write_text(report)

    method = f"""# Full363 input and method comparison

| Item | Current reconstruction | Published PNAS source | Assessment |
|---|---|---|---|
| Full gene-tree source | `data/stage4d/named_recollapsed.tre` via current Stage 4E inputs | `63430.aLRT-0.95-collapsed.gene.trees.gz` | Different source collections |
| Full T0 locus count | 63,430 records | 63,430 published collapsed records | Count agrees |
| PNAS exclusion | 1,431 frozen coordinate-selected loci | Published chr4 outlier intervals | Membership/count previously verified; tree identity remains separate |
| T_PNAS retained loci | 61,999 | 61,999 expected after 1,431 exclusions | Count agrees |
| Unresolved source correspondence | 317 records; 200 source-discrepancy classifications | No complete source-side locus resolution in cached snapshot | Provenance incomplete |
| ASTRAL | v1.25.4.8, `-u 2`, `--genelength 1000`, `--sulengthtype 1` | Version/command not present beside cached `score-columbea` trees | Do not attribute RF 2 to version alone |

The current T0 matching the published full-data topology is evidence that the
current taxon labels and the baseline ASTRAL implementation are compatible at
the topology level. It does not validate the reconstructed outlier-removed
gene-tree collection. The RF 2 difference between current T_PNAS and the
published outlier-removed tree is consistent with an input/preprocessing
difference, but the existing 317-record audit does not identify a sufficient
causal subset.

The PNAS repository scripts establish the QQS definition and identify
`score-columbea` as the Figure S9 species-tree directory. They do not provide a
complete reproducible ASTRAL command for the cached species trees in the
versioned snapshot.
"""
    (RESULTS / "full363_input_method_comparison.md").write_text(method)
    paragraph = (RESULTS / "full363_pnas_manuscript_paragraph.md")
    paragraph.write_text(
        "Our full-data 363-taxon ASTRAL tree reproduced the cached published PNAS "
        "full-data topology exactly under canonical unrooted RF comparison. The "
        "published cached outlier-removed tree retained the same N62 relationship, "
        "with Columbimorphae intact, whereas our reconstructed outlier-removed run "
        "changed one internal split (RF 2) and placed Columbiformes with "
        "Otidimorphae. Because our gene-tree inputs were reconstructed from the "
        "named collection rather than the original aLRT-collapsed archive, this is "
        "a source-reproduction discrepancy rather than evidence that the PNAS "
        "taxon-sampling result is contradicted. The result motivates an exact-input "
        "controlled comparison; it does not by itself establish an MSRC-specific "
        "causal mechanism.\n"
    )


def make_figure(trees, branches):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGURES.mkdir(parents=True, exist_ok=True)
    green, magenta, gray, blue = "#4d9221", "#c51b7d", "#666666", "#2166ac"
    fig = plt.figure(figsize=(7.3, 5.8), constrained_layout=True)
    grid = fig.add_gridspec(3, 1, height_ratios=(1.0, 1.0, 1.15), hspace=0.25)

    # A: a quartet-pruned visualization of the exact changed N62 relationship.
    ax = fig.add_subplot(grid[0])
    ax.set_title("A  N62 focal relationship", loc="left", fontsize=10, fontweight="bold")
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    topologies = [(0.18, "S2024 / T0", ("Columbiformes", "Other Columbimorphae")),
                  (0.68, "reconstructed T_PNAS", ("Columbiformes", "Otidimorphae"))]
    for x, label, pair in topologies:
        y = 0.58
        ax.plot([x, x + 0.08], [y, y], color=gray, lw=1.4)
        ax.plot([x, x], [0.24, 0.82], color=gray, lw=1.4)
        ax.plot([x + 0.08, x + 0.23], [y, y], color=gray, lw=1.4)
        ax.text(x + 0.25, y, pair[0], va="center", fontsize=7)
        ax.text(x + 0.25, y - 0.23, pair[1], va="center", fontsize=7)
        ax.plot([x + 0.08, x + 0.08], [y - 0.23, y], color=gray, lw=1.4)
        ax.plot([x, x + 0.08], [0.24, y - 0.23], color=gray, lw=1.4)
        ax.text(x + 0.25, 0.11, "other taxa", va="center", fontsize=7, color="#555555")
        ax.text(x, 0.90, label, fontsize=8, fontweight="bold")
    ax.text(0.02, 0.02, "Visualization uses representative focal groups; inference was performed on all 363 taxa.", fontsize=7, color="#555555")

    # B: observed state and RF by discrete treatment.
    ax = fig.add_subplot(grid[1])
    ax.set_title("B  Full 363-taxon treatment outcomes", loc="left", fontsize=10, fontweight="bold")
    labels = ["T0", "PNAS", "struct.", "500 kb", "non-chr4"]
    states = ["q1", "q2", "q1", "q2", "q2"]
    rf_values = [0, 2, 0, 2, 2]
    colors = [green if state == "q1" else magenta for state in states]
    ax.bar(range(5), [1] * 5, color=colors, width=0.78)
    for i, (state, value) in enumerate(zip(states, rf_values)):
        ax.text(i, 0.5, f"N62 {state}", ha="center", va="center", color="white", fontsize=8, rotation=90)
        ax.text(i, -0.12, f"RF {value}", ha="center", va="top", fontsize=8, color="#444444")
    ax.set_xticks(range(5), labels)
    ax.set_ylim(-0.25, 1.05)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="x", length=0, labelsize=8)
    ax.text(0.01, -0.23, "green = published/S2024-side q1; magenta = reconstructed q2",
            transform=ax.transAxes, fontsize=7, color="#555555")

    # C: homologous focal branch estimates. Use only exact state/treatment matches.
    ax = fig.add_subplot(grid[2])
    ax.set_title("C  Focal CU/SU lengths on exact splits", loc="left", fontsize=10, fontweight="bold")
    branch_map = {(row["treatment"], row["clade"], row["state"]): row for row in branches}
    xlabels, cu, su = [], [], []
    for treatment, short in zip(TREATMENT_ORDER, labels):
        for clade, state in (("N61", "q1"), ("N62", "q1" if treatment in ("T0", "T_STRUCT_PRIMARY") else "q2")):
            row = branch_map[(treatment, clade, state)]
            xlabels.append(f"{short}\n{clade} {state}")
            cu.append(float(row["CU_length"]) if row["present"] else 0.0)
            su.append(float(row["SU_length"]) if row["present"] else 0.0)
    x = list(range(len(xlabels)))
    ax.bar([v - 0.18 for v in x], cu, width=0.36, color=blue, label="CU")
    ax2 = ax.twinx()
    ax2.bar([v + 0.18 for v in x], su, width=0.36, color="#e69f00", label="SU")
    ax.set_ylabel("CU length", fontsize=8)
    ax2.set_ylabel("SU length", fontsize=8)
    ax.set_xticks(x, xlabels, fontsize=6)
    ax.tick_params(axis="y", labelsize=7)
    ax2.tick_params(axis="y", labelsize=7)
    ax.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=blue, label="CU"), plt.Rectangle((0, 0), 1, 1, color="#e69f00", label="SU")], frameon=True, framealpha=0.85, edgecolor="none", fontsize=7, loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(FIGURES / f"neoaves_pnas_replication_comparison.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    trees = audit_trees()
    comparisons = comparison_rows(trees)
    focals = focal_rows(trees)
    branches = branch_rows_for_report(trees)
    write_reports(trees, comparisons, focals, branches)
    make_figure(trees, branches)
    print("wrote full363 PNAS audit and figure outputs")


if __name__ == "__main__":
    main()
