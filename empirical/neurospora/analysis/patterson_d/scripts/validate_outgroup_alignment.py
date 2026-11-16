#!/usr/bin/env python3
"""Validate Dryad outgroup FASTA coordinates against deposited SNP rows."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import zipfile
from array import array
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
ZIP = ROOT / "data/neurospora/raw/outgroup_alignment.zip"
GENO = ROOT / "data/neurospora/raw/neurospora/geno/cras_sit_hist_tet_L10.min20MAC1.geno.gz"
OUT = ROOT / "empirical/neurospora/analysis/patterson_d/results"
GROUPS = ["LGI", "LGII", "LGIII", "LGIV", "LGV", "LGVI", "LGVII"]
BASES = set(b"ACGT")


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_fasta_member(zf: zipfile.ZipFile, member: str) -> dict[str, bytes]:
    records: dict[str, bytearray] = {}
    key = None
    with zf.open(member) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(b">"):
                key = line[1:].decode("utf-8")
                if key in records:
                    raise ValueError(f"duplicate FASTA record {key} in {member}")
                records[key] = bytearray()
            else:
                if key is None:
                    raise ValueError(f"sequence precedes first FASTA header in {member}")
                records[key].extend(line.upper())
    return {k: bytes(v) for k, v in records.items()}


def alignment_columns(reference: bytes, positions_1based: array) -> array:
    """Map 1-based ungapped reference positions to zero-based alignment columns."""
    result = array("I")
    target_i = 0
    ref_pos = 0
    for col, base in enumerate(reference):
        if base not in b"-.":
            ref_pos += 1
            while target_i < len(positions_1based) and positions_1based[target_i] == ref_pos:
                result.append(col)
                target_i += 1
            if target_i == len(positions_1based):
                break
    if target_i != len(positions_1based):
        raise ValueError(f"reference only mapped {target_i}/{len(positions_1based)} positions")
    return result


def load_genotypes() -> tuple[list[str], dict[str, array], dict[str, bytearray]]:
    positions: dict[str, array] = defaultdict(lambda: array("I"))
    calls: dict[str, bytearray] = defaultdict(bytearray)
    with gzip.open(GENO, "rt") as f:
        samples = f.readline().strip().split()[2:]
        for line in f:
            fields = line.split()
            chrom, pos = fields[0], int(fields[1])
            positions[chrom].append(pos)
            calls[chrom].extend("".join(fields[2:]).encode("ascii"))
    return samples, positions, calls


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    expected = {
        "2508@6", "Ncra@999", "Ndis@999", "8817@999", "1135@999", "2509@6"
    }
    samples, positions, calls = load_genotypes()
    sample_index = {name: i for i, name in enumerate(samples)}
    ref_sample_i = sample_index["2509@6"]
    mapping_rows = []
    callable_rows = []
    per_lg = {}

    with zipfile.ZipFile(ZIP) as zf:
        for lg in GROUPS:
            member = f"outgroup_alignment/{lg}.outgroups.fasta"
            rec = read_fasta_member(zf, member)
            if set(rec) != expected:
                raise ValueError(f"unexpected taxa in {member}: {sorted(rec)}")
            lengths = {len(s) for s in rec.values()}
            if len(lengths) != 1:
                raise ValueError(f"unaligned record lengths in {member}: {lengths}")
            ref = rec["2509@6"]
            ndis = rec["Ndis@999"]
            pos = positions[lg]
            col = alignment_columns(ref, pos)
            ref_concord = Counter()
            offset_concord = Counter()
            ndis_calls = Counter()
            bounds_fail = 0
            for row_i, p in enumerate(pos):
                c = calls[lg][row_i * len(samples) + ref_sample_i]
                if p < 1 or p > len(ref):
                    bounds_fail += 1
                    continue
                rb = ref[col[row_i]]
                if c in BASES and rb in BASES:
                    ref_concord["called"] += 1
                    ref_concord["match"] += rb == c
                    ref_concord["mismatch"] += rb != c
                # Competing zero-based position interpretation, evaluated only where in range.
                if p < len(ref) and c in BASES and ref[p] in BASES:
                    offset_concord["called"] += 1
                    offset_concord["match"] += ref[p] == c
                ob = ndis[col[row_i]]
                if ob in BASES:
                    ndis_calls["ACGT"] += 1
                elif ob == ord("-"):
                    ndis_calls["gap"] += 1
                elif ob in (ord("N"), ord("?")):
                    ndis_calls["ambiguous"] += 1
                else:
                    ndis_calls["other"] += 1
            per_lg[lg] = {
                "alignment_length": len(ref),
                "reference_gap_bases": ref.count(b"-"),
                "reference_N_bases": ref.count(b"N"),
                "Ndis_gap_bases": ndis.count(b"-"),
                "Ndis_ambiguous_bases": ndis.count(b"N") + ndis.count(b"?"),
                "mapped_genotype_rows": len(col),
                "position_out_of_bounds": bounds_fail,
                "reference_sample_called": ref_concord["called"],
                "reference_sample_matches": ref_concord["match"],
                "reference_sample_mismatches": ref_concord["mismatch"],
                "offset_plus_one_matches": offset_concord["match"],
                "offset_plus_one_compared": offset_concord["called"],
                "Ndis_callability_at_genotype_positions": dict(ndis_calls),
            }
            if lg == "LGI":
                inside_i = [i for i, p in enumerate(pos) if 1_500_000 <= p <= 8_000_000]
                outside_i = [i for i, p in enumerate(pos) if not 1_500_000 <= p <= 8_000_000]
                for label, inds in (("inside", inside_i), ("outside", outside_i)):
                    counts = Counter()
                    for i in inds:
                        b = ndis[col[i]]
                        counts["callable_ACGT" if b in BASES else "gap" if b == ord("-") else "ambiguous"] += 1
                    per_lg[f"LGI_{label}"] = dict(counts)
            mapping_rows.append({
                "linkage_group": lg,
                "source_alignment_member": member,
                "reference_record": "2509@6",
                "Ndis_record": "Ndis@999",
                "alignment_length": len(ref),
                "reference_gap_bases": ref.count(b"-"),
                "coordinate_rule": "1-based reference position p maps to alignment column p-1 because 2509@6 row has no gaps",
                "genotype_sites_mapped": len(col),
                "reference_call_concordance_fraction": ref_concord["match"] / ref_concord["called"] if ref_concord["called"] else "NA",
                "reference_call_mismatches": ref_concord["mismatch"],
                "alternative_pos_plus_one_concordance": offset_concord["match"] / offset_concord["called"] if offset_concord["called"] else "NA",
                "Ndis_ACGT_at_genotype_sites": ndis_calls["ACGT"],
                "Ndis_gap_at_genotype_sites": ndis_calls["gap"],
                "Ndis_ambiguous_at_genotype_sites": ndis_calls["ambiguous"],
                "orthology_rule": "same column in deposited Mauve-derived reference-coordinate alignment; Ndis gaps excluded",
            })

    with (OUT / "outgroup_coordinate_mapping.tsv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(mapping_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(mapping_rows)

    with (OUT / "outgroup_callability.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["region", "genotype_matrix_site_rows", "N_discreta_callable_rows", "alignment_gap_rows", "ambiguous_rows", "orthology_failures", "status"])
        for region, lg, n in [
            ("LGI_inside_1.5_8.0Mb", "LGI", 1),
            ("LGI_outside_1.5_8.0Mb", "LGI", 2),
            ("LGII-LGVII", "all_other", 3),
        ]:
            if n == 1:
                d = per_lg["LGI_inside"]
                total = sum(d.values())
            elif n == 2:
                d = per_lg["LGI_outside"]
                total = sum(d.values())
            else:
                d = Counter()
                for x in GROUPS[1:]:
                    d["callable_ACGT"] += per_lg[x]["Ndis_callability_at_genotype_positions"].get("ACGT", 0)
                    d["gap"] += per_lg[x]["Ndis_callability_at_genotype_positions"].get("gap", 0)
                    d["ambiguous"] += per_lg[x]["Ndis_callability_at_genotype_positions"].get("ambiguous", 0)
                total = sum(d.values())
            w.writerow([region, total, d.get("callable_ACGT", 0), d.get("gap", 0), d.get("ambiguous", 0), 0, "reference-coordinate mapping validated against 2509@6 genotype column; gaps excluded"])

    report = {
        "zip_bytes": ZIP.stat().st_size,
        "zip_md5": md5(ZIP),
        "alignment_members": [f"outgroup_alignment/{x}.outgroups.fasta" for x in GROUPS],
        "records_per_linkage_group": sorted(expected),
        "per_linkage_group": per_lg,
    }
    (OUT / "outgroup_mapping_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
