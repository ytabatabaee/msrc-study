#!/usr/bin/env python3
"""Summarize optimized origin-placement validation outputs."""
from __future__ import annotations
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=Path, required=True)
    args = parser.parse_args()
    forward = read(args.prefix.with_name(args.prefix.name + "_forward.tsv"))
    quartet = read(args.prefix.with_name(args.prefix.name + "_quartets.tsv"))
    groups = defaultdict(list)
    for row in forward:
        groups[(row["condition"], row["replicate"], row["scenario"])].append(row)
    persistence = []
    for key, rows in sorted(groups.items()):
        n = len(rows)
        mean = sum(float(r["origin_polymorphic_fraction"]) for r in rows) / n
        se = (sum((float(r["origin_polymorphic_fraction"]) - mean) ** 2 for r in rows) / (n - 1)) ** 0.5 / n ** 0.5 if n > 1 else float("nan")
        persistence.append({"condition": key[0], "replicate": key[1], "scenario": key[2], "n_histories": n,
                            "mean_origin_polymorphic_fraction": mean, "history_se": se,
                            "origin_end_polymorphic_fraction": sum(int(r["origin_end_polymorphic"]) for r in rows) / n,
                            "ever_fixed_fraction": sum(int(r["ever_fixed"]) for r in rows) / n,
                            "ever_lost_fraction": sum(int(r["ever_lost"]) for r in rows) / n})
    write(args.prefix.with_name(args.prefix.name + "_persistence.tsv"), persistence, list(persistence[0]))
    qgroups = defaultdict(lambda: defaultdict(int))
    cqgroups = defaultdict(lambda: defaultdict(int))
    for row in quartet:
        qgroups[(row["condition"], row["replicate"], row["scenario"], row["rho"], row["quartet"])][row["topology_class"]] += int(row["count"])
        cqgroups[(row["condition"], row["replicate"], row["scenario"], row["rho"], row["quartet"], row["pattern"])][row["topology_class"]] += int(row["count"])
    marginal, conditional = [], []
    for key, counts in sorted(qgroups.items()):
        total = sum(counts.values())
        ps = counts.get("species_tree", 0) / total
        pd = counts.get("discordant", 0) / total
        row = {"condition": key[0], "replicate": key[1], "scenario": key[2], "rho": key[3], "quartet": key[4],
               "p_species_tree": ps, "p_discordant": pd, "margin": ps - pd, "n_genealogies": total}
        marginal.append(row)
    for key, counts in sorted(cqgroups.items()):
        total = sum(counts.values())
        ps = counts.get("species_tree", 0) / total
        pd = counts.get("discordant", 0) / total
        conditional.append({"condition": key[0], "replicate": key[1], "scenario": key[2], "rho": key[3],
                            "quartet": key[4], "pattern": key[5], "p_species_tree": ps,
                            "p_discordant": pd, "margin": ps - pd, "n_genealogies": total})
    write(args.prefix.with_name(args.prefix.name + "_marginal_quartets.tsv"), marginal, list(marginal[0]) if marginal else ["condition"])
    write(args.prefix.with_name(args.prefix.name + "_conditional_quartets.tsv"), conditional, list(conditional[0]) if conditional else ["condition"])

if __name__ == "__main__":
    main()
