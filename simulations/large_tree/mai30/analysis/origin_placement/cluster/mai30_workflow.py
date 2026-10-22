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


class CheckpointRequiresRegeneration(ValueError):
    """Raised when a legacy checkpoint lacks provenance needed for replay."""

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
    from msrcsim.speciation_split import SplitInitializationEvent  # type: ignore
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
    split = induced_quartet_split(newick, quartet)
    if split is None:
        return "unresolved"
    left, right = split
    return f"{left[0]}|{left[1]}::{right[0]}|{right[1]}"


def induced_quartet_split(newick: str, quartet: Iterable[str]) -> tuple[tuple[str, str], tuple[str, str]] | None:
    quartet = tuple(quartet)
    leaves, splits = parse_splits(newick)
    if not set(quartet) <= leaves:
        return None
    q = set(quartet)
    for split in splits:
        part = split & q
        if len(part) == 2:
            pair = tuple(sorted(part))
            other = tuple(sorted(q - part))
            return (pair, other) if pair <= other else (other, pair)
    return None


def canonical_quartet_topologies(quartet: Iterable[str]) -> tuple[str, str, str]:
    taxa = tuple(quartet)
    if len(set(taxa)) != 4:
        raise ValueError("a quartet must contain four unique taxa")
    result = []
    for left_indices in ((0, 1), (0, 2), (0, 3)):
        left = tuple(sorted((taxa[left_indices[0]], taxa[left_indices[1]])))
        right = tuple(sorted(set(taxa) - set(left)))
        result.append(f"{left[0]}|{left[1]}::{right[0]}|{right[1]}" if left <= right else f"{right[0]}|{right[1]}::{left[0]}|{left[1]}")
    return tuple(result)


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
        "split_events": [asdict(event) if hasattr(event, "__dataclass_fields__") else dict(event.__dict__) for event in getattr(history, "split_events", [])],
        "records": records,
    }


def write_gzip_json(path: Path, value: Any) -> None:
    data = gzip.compress((json.dumps(value, sort_keys=True) + "\n").encode(), mtime=0)
    atomic_write(path, data)


def read_gzip_json(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt") as handle:
        return json.load(handle)


def validate_history_payload(payload: dict[str, Any], row: dict[str, str], history_id: int) -> None:
    if "history_id" not in payload and "history" in payload:
        raise CheckpointRequiresRegeneration("legacy checkpoint uses history, not history_id; regenerate on Bridges2")
    if "sampled" not in payload and "sampled_arrangements" in payload:
        raise CheckpointRequiresRegeneration("legacy checkpoint uses sampled_arrangements; regenerate with the validated schema")
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
    scenario = scenario_row(row["condition"], row["replicate"], row["scenario"])
    if "origin_branch" in payload and payload["origin_branch"] != scenario["origin_branch"]:
        raise ValueError("history checkpoint origin branch mismatch")
    if "origin_time_from_branch_start" in payload and int(payload["origin_time_from_branch_start"]) != int(scenario["origin_time_from_branch_start"]):
        raise ValueError("history checkpoint origin time mismatch")
    taxa = set(payload["taxa"])
    if taxa != set(payload["sampled"]):
        raise ValueError("history checkpoint sampled taxa mismatch")
    if not isinstance(payload["records"], list):
        raise ValueError("history records are not a list")
    validate_record_states(payload["records"])


def boundary_observations(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """Return recorded parent/daughter boundary states without interpolation."""
    observations = []
    for record in payload["records"]:
        count = int(record["copy_count_A1"])
        total = int(record["population_chromosomes"])
        for role in (["start"] if record["is_branch_start"] else []) + (["end"] if record["is_branch_end"] else []):
            observations.append({
                "boundary_id": f"{record['branch_id']}:{role}",
                "branch_id": record["branch_id"],
                "parent_branch_id": record.get("parent_branch_id"),
                "role": role,
                "absolute_age": float(record["absolute_age"]),
                "copy_count_A1": count,
                "population_chromosomes": total,
                "arrangement_frequency": count / total,
                "polymorphic": int(0 < count < total),
                "fixed": int(count == total),
                "lost": int(count == 0),
            })
    if len({item["boundary_id"] for item in observations}) != len(observations):
        raise ValueError("duplicate population-boundary observations in checkpoint")
    return observations


def sampled_pattern(payload: dict[str, Any]) -> str:
    taxa = sorted(payload["taxa"])
    return "".join(str(int(payload["sampled"][taxon])) for taxon in taxa)


def validate_record_states(records: list[dict[str, Any]]) -> None:
    for record in records:
        required_record = {"branch_id", "absolute_age", "copy_count_A1", "copy_count_A0", "population_chromosomes", "frequency_A1", "frequency_A0", "is_branch_start", "is_branch_end"}
        if required_record - set(record):
            raise CheckpointRequiresRegeneration("checkpoint records lack complete population-boundary state; regenerate on Bridges2")
        total = int(record["population_chromosomes"])
        raw_count = float(record["copy_count_A1"])
        raw_count_a0 = float(record["copy_count_A0"])
        count = int(raw_count)
        count_a0 = int(raw_count_a0)
        if raw_count != count or raw_count_a0 != count_a0:
            raise ValueError("arrangement count is not an integer")
        if total <= 0 or count < 0 or count > total or count_a0 != total - count:
            raise ValueError("invalid arrangement/population count in history checkpoint")
        if abs(float(record["frequency_A1"]) - count / total) > 1e-10 or abs(float(record["frequency_A0"]) - count_a0 / total) > 1e-10:
            raise ValueError("checkpoint frequency does not match recorded arrangement count")


def rebuild_history(payload: dict[str, Any], row: dict[str, str], tree: Any) -> Any:
    api = load_simulator()
    records = [api["FrequencyRecord"](**raw) for raw in payload["records"]]
    by_branch: dict[str, list[Any]] = {}
    for record in records:
        by_branch.setdefault(record.branch_id, []).append(record)
    split_events = [api["SplitInitializationEvent"](**raw) for raw in payload.get("split_events", [])]
    return api["FrequencyHistory"](
        records=records, by_branch=by_branch,
        split_events=split_events,
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
