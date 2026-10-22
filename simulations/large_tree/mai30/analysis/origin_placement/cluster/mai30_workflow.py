#!/usr/bin/env python3
"""Shared, lazy-loading implementation for the Mai30 Bridges2 workflow.

The module is safe to import on a laptop: importing it never imports msrc-sim
and never starts a simulation.  The simulator is loaded only by the forward,
backward, and panel stages after the orchestration driver has verified a live
Slurm compute allocation.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable

HERE = Path(__file__).resolve().parent
MAI = HERE.parents[2]
MANIFEST = HERE / "mai30_bridges2_manifest.tsv"
EXPECTED_SIMULATOR_COMMIT = "bfba58d627f6f467af649cdb0b37f969349e0aab"
BASE_SEED = 20261022
RHOS = ("0", "0.0001", "0.001", "0.002", "0.01", "0.02")
TREE_SPECS = (
    ("outgroup.0.species.1.5.genes.1.5", "007"),
    ("outgroup.0.species.1.5.genes.1.5", "042"),
    ("outgroup.0.species.1.5.genes.1.5", "074"),
    ("outgroup.0.species.5.genes.5", "005"),
)


def rows() -> list[dict[str, str]]:
    with MANIFEST.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def selected_row(stage: str, task_id: int) -> dict[str, str]:
    candidates = [r for r in rows() if r["stage"] == stage]
    if task_id < 0 or task_id >= len(candidates):
        raise ValueError(f"task id {task_id} is outside {stage} range")
    return candidates[task_id]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return sha256_bytes(text.encode())


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent))
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def atomic_json(path: Path, value: Any) -> None:
    atomic_write(path, (json.dumps(value, indent=2, sort_keys=True) + "\n").encode())


def simulator_root() -> Path:
    value = os.environ.get("MSRC_SIM_ROOT")
    if not value:
        raise RuntimeError("MSRC_SIM_ROOT must point to the pinned msrc-sim checkout")
    root = Path(value).expanduser().resolve()
    if not (root / ".git").exists():
        raise RuntimeError(f"MSRC_SIM_ROOT is not a git checkout: {root}")
    actual = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    if actual != EXPECTED_SIMULATOR_COMMIT:
        raise RuntimeError(f"msrc-sim commit mismatch: expected {EXPECTED_SIMULATOR_COMMIT}, found {actual}")
    return root


def require_live_compute_allocation() -> None:
    """Fail closed even when an adapter is invoked without the driver."""
    required = ("SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_JOB_PARTITION", "SLURM_CPUS_PER_TASK")
    if any(not os.environ.get(name) for name in required):
        raise RuntimeError("computational adapter requires a live Slurm allocation")
    scontrol = shutil.which("scontrol")
    if not scontrol:
        raise RuntimeError("computational adapter requires scontrol")
    try:
        info = subprocess.check_output([scontrol, "show", "job", os.environ["SLURM_JOB_ID"]], text=True, stderr=subprocess.STDOUT)
        hosts = subprocess.check_output([scontrol, "show", "hostnames", os.environ["SLURM_JOB_NODELIST"]], text=True, stderr=subprocess.STDOUT).split()
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RuntimeError("Slurm could not verify the active allocation") from exc
    if "JobState=RUNNING" not in info:
        raise RuntimeError("Slurm allocation is not RUNNING")
    hostname = socket.gethostname().split(".")[0]
    if hostname not in {host.split(".")[0] for host in hosts}:
        raise RuntimeError("current host is outside the allocated compute nodes")


def load_simulator() -> dict[str, Any]:
    """Load the public pinned APIs only inside a validated compute stage."""
    root = simulator_root()
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    from msrcsim import SpeciesTree  # type: ignore
    from msrcsim.population_process import simulate_population_history  # type: ignore
    from msrcsim.rearrangement import Rearrangement  # type: ignore
    from msrcsim.structured_coalescent import simulate_genealogy  # type: ignore
    from msrcsim.wright_fisher import FrequencyHistory, FrequencyRecord  # type: ignore
    return locals()


def config_path(condition: str, replicate: str) -> Path:
    return MAI / "datasets/configs/mai30_msrc_pilot" / f"{condition}__replicate_{replicate}__msrc_no_suppression__freq_0.5.yaml"


def load_yaml(path: Path) -> dict[str, Any]:
    import yaml  # lazy; only stage adapters need it
    with path.open() as handle:
        return yaml.safe_load(handle)


def scenario_row(condition: str, replicate: str, scenario: str) -> dict[str, str]:
    source = MAI / "analysis/origin_placement/origin_scenarios.tsv"
    with source.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["condition"] == condition and row["replicate"] == replicate and row["scenario"] == scenario:
                return row
    raise RuntimeError(f"no origin scenario for {condition}/{replicate}/{scenario}")


def config_metadata(row: dict[str, str]) -> dict[str, str]:
    config = config_path(row["condition"], row["replicate"])
    scenario = scenario_row(row["condition"], row["replicate"], row["scenario"])
    if not config.exists():
        raise RuntimeError(f"missing Mai30 config: {config}")
    return {
        "config_sha256": sha256_file(config),
        "scenario_sha256": sha256_text(json.dumps(scenario, sort_keys=True)),
        "species_tree_sha256": sha256_text(load_yaml(config)["species_tree"]["newick"]),
    }


def config_key(row: dict[str, str]) -> str:
    return f"{row['condition']}__replicate_{row['replicate']}__{row['scenario']}"


def output_root() -> Path:
    return Path(os.environ.get("MAI_CLUSTER_OUT", str(HERE / "bridges2_outputs"))).expanduser().resolve()


def kind_for(row: dict[str, str]) -> str:
    return "benchmark" if row["stage"].startswith("benchmark_") else "production"


def forward_dir(row: dict[str, str], root: Path | None = None) -> Path:
    root = root or output_root()
    return root / kind_for(row) / "forward" / config_key(row)


def history_path(row: dict[str, str], history_id: int, root: Path | None = None) -> Path:
    return forward_dir(row, root) / f"history_{history_id:04d}.json.gz"


def backward_dir(row: dict[str, str], root: Path | None = None) -> Path:
    root = root or output_root()
    return root / kind_for(row) / "backward" / config_key(row) / f"rho_{row['rho']}"


def task_dir(stage: str, task_id: int, root: Path | None = None) -> Path:
    return (root or output_root()) / "tasks" / stage / str(task_id)


def task_manifest(stage: str, task_id: int, root: Path | None = None) -> Path:
    return task_dir(stage, task_id, root) / "output_manifest.json"


def seed_for(tree_index: int, scenario_index: int, history_id: int) -> int:
    return BASE_SEED + tree_index * 10_000_000 + scenario_index * 100_000 + history_id


def genealogy_seed(tree_index: int, scenario_index: int, rho_index: int, history_id: int) -> int:
    return BASE_SEED + 7_000_000_000 + tree_index * 100_000_000 + scenario_index * 10_000_000 + rho_index * 100_000 + history_id


def tree_and_scenario(row: dict[str, str]) -> tuple[int, int]:
    tree_index = next((i for i, spec in enumerate(TREE_SPECS) if spec == (row["condition"], row["replicate"])), None)
    if tree_index is None:
        raise ValueError(f"unknown Mai30 tree: {row['condition']}/{row['replicate']}")
    if row["scenario"] not in {"A", "B"}:
        raise ValueError(f"only published origin scenarios A and B are supported: {row['scenario']}")
    return tree_index, (0 if row["scenario"] == "A" else 1)


def leaf_names(newick: str) -> list[str]:
    # Labels are the tokens after '(' or ',' and before a delimiter.  Support
    # labels and branch lengths are skipped by the delimiter rule.
    return re.findall(r"(?:\(|,)\s*([^(),:;\s]+)", newick)


def validate_newick(newick: str, expected_taxa: Iterable[str]) -> None:
    if not newick.strip().endswith(";") or newick.count("(") != newick.count(")"):
        raise ValueError("gene tree is not a complete Newick record")
    observed = leaf_names(newick)
    if len(observed) != len(set(observed)):
        raise ValueError("gene tree contains duplicate taxon labels")
    if set(observed) != set(expected_taxa):
        raise ValueError(f"gene tree taxa differ from expected set: {sorted(set(observed) ^ set(expected_taxa))}")


def parse_splits(newick: str) -> tuple[set[str], set[frozenset[str]]]:
    tokens = re.findall(r"\(|\)|,|:|;|[^(),:;\s]+", newick)
    pos = 0

    def subtree() -> tuple[set[str], set[frozenset[str]]]:
        nonlocal pos
        if tokens[pos] == "(":
            pos += 1
            children = [subtree()]
            while tokens[pos] == ",":
                pos += 1
                children.append(subtree())
            if tokens[pos] != ")":
                raise ValueError("malformed Newick")
            pos += 1
            if pos < len(tokens) and tokens[pos] not in {",", ")", ":", ";"}:
                pos += 1  # optional internal support/name
            if pos < len(tokens) and tokens[pos] == ":":
                pos += 1
                if pos >= len(tokens):
                    raise ValueError("missing branch length")
                pos += 1
            leaves = set().union(*(x[0] for x in children))
            splits = set().union(*(x[1] for x in children))
            for child, _ in children:
                if 1 < len(child) < len(leaves):
                    splits.add(frozenset(child))
            return leaves, splits
        leaf = tokens[pos]
        pos += 1
        if pos < len(tokens) and tokens[pos] == ":":
            pos += 2
        return {leaf}, set()

    leaves, splits = subtree()
    if pos >= len(tokens) or tokens[pos] != ";":
        raise ValueError("Newick terminator missing")
    return leaves, splits


def quartet_topology(newick: str, quartet: tuple[str, str, str, str]) -> str:
    leaves, splits = parse_splits(newick)
    if not set(quartet) <= leaves:
        return "missing"
    q = set(quartet)
    for split in splits:
        part = split & q
        if len(part) == 2:
            pair = tuple(sorted(part))
            other = tuple(sorted(q - part))
            return f"{pair[0]}|{pair[1]}::{other[0]}|{other[1]}"
    return "unresolved"


def history_payload(row: dict[str, str], history_id: int, tree: Any, history: Any, sampled: dict[str, int], seed: int, metadata: dict[str, str]) -> dict[str, Any]:
    records = [asdict(record) if hasattr(record, "__dataclass_fields__") else dict(record.__dict__) for record in history.records]
    return {
        "schema_version": 1,
        "stage": "forward_history",
        "history_id": history_id,
        "condition": row["condition"],
        "replicate": row["replicate"],
        "scenario": row["scenario"],
        "seed": seed,
        "simulator_commit": EXPECTED_SIMULATOR_COMMIT,
        **metadata,
        "taxa": sorted(tree.taxa),
        "sampled": sampled,
        "terminal_frequencies": {taxon: float(history.terminal_frequency(taxon)) for taxon in tree.taxa},
        "structured_origin_age": history.structured_origin_age,
        "structured_origin_branch": history.structured_origin_branch,
        "records": records,
    }


def write_gzip_json(path: Path, value: Any) -> None:
    data = gzip.compress((json.dumps(value, sort_keys=True) + "\n").encode(), mtime=0)
    atomic_write(path, data)


def read_gzip_json(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt") as handle:
        return json.load(handle)


def validate_history_payload(payload: dict[str, Any], row: dict[str, str], history_id: int) -> None:
    metadata = config_metadata(row)
    required = {"schema_version", "history_id", "condition", "replicate", "scenario", "seed", "simulator_commit", "taxa", "sampled", "terminal_frequencies", "structured_origin_age", "structured_origin_branch", "records", *metadata}
    missing = required - set(payload)
    if missing:
        raise ValueError(f"history checkpoint missing keys: {sorted(missing)}")
    if payload["history_id"] != history_id or payload["simulator_commit"] != EXPECTED_SIMULATOR_COMMIT:
        raise ValueError("history checkpoint identity or simulator commit mismatch")
    for key, value in metadata.items():
        if payload[key] != value:
            raise ValueError(f"history checkpoint {key} mismatch")
    if payload["condition"] != row["condition"] or payload["replicate"] != row["replicate"] or payload["scenario"] != row["scenario"]:
        raise ValueError("history checkpoint configuration mismatch")
    taxa = set(payload["taxa"])
    if taxa != set(payload["sampled"]):
        raise ValueError("history checkpoint sampled taxa mismatch")
    if not isinstance(payload["records"], list):
        raise ValueError("history records are not a list")


def rebuild_history(payload: dict[str, Any], row: dict[str, str], tree: Any) -> Any:
    api = load_simulator()
    records = [api["FrequencyRecord"](**raw) for raw in payload["records"]]
    by_branch: dict[str, list[Any]] = {}
    for record in records:
        by_branch.setdefault(record.branch_id, []).append(record)
    return api["FrequencyHistory"](
        records=records, by_branch=by_branch,
        structured_origin_age=payload["structured_origin_age"],
        structured_origin_branch=payload["structured_origin_branch"],
    )


def file_manifest(base: Path, files: list[Path]) -> list[dict[str, Any]]:
    base = base.resolve()
    result = []
    for path in files:
        path = path.resolve()
        if not path.exists() or path.stat().st_size == 0:
            raise ValueError(f"missing or empty output: {path}")
        result.append({"path": str(path.relative_to(base)), "size": path.stat().st_size, "sha256": sha256_file(path)})
    return result


def output_manifest(path: Path, stage: str, task_id: int, row: dict[str, str], files: list[Path], expected_count: int, actual_count: int, ids: list[str] | None = None, columns: list[str] | None = None, extra: dict[str, Any] | None = None) -> None:
    base = output_root()
    value: dict[str, Any] = {
        "schema_version": 1,
        "status": "complete",
        "stage": stage,
        "task_id": task_id,
        "condition": row["condition"],
        "replicate": row["replicate"],
        "scenario": row["scenario"],
        "rho": row["rho"],
        "simulator_commit": EXPECTED_SIMULATOR_COMMIT,
        "expected_count": expected_count,
        "actual_count": actual_count,
        "files": file_manifest(base, files),
    }
    if ids is not None:
        if len(ids) != len(set(ids)):
            raise ValueError("duplicate identifiers in output manifest")
        value["ids"] = ids
    if columns is not None:
        value["columns"] = columns
    if extra:
        value.update(extra)
    atomic_json(path, value)
