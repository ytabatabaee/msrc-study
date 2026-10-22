#!/usr/bin/env python3
"""Print safe, non-submitting Bridges-2 Stage 4E commands."""
from __future__ import annotations
import argparse

GROUPS = {
    "primary": ["T0", "T_PNAS", "T_STRUCT_PRIMARY", "T_BLOCK_500K", "NONCHR4"],
    "sensitivity": ["T_STRUCT_STRINGENT", "T_STRUCT_INCLUSIVE", "T_BLOCK_250K", "T_BLOCK_1MB"],
    "j48": ["J48_T0", "J48_T_PNAS", "J48_T_STRUCT_PRIMARY", "J48_T_BLOCK_500K", "J48_NONCHR4"],
}

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--group", choices=[*GROUPS, "random-matched"], default="primary")
    parser.add_argument("--script", default="empirical/neoaves_chr4/cluster/stage4e_aster/01_single_treatment.sbatch")
    parser.add_argument("--account", default="bio240015p")
    parser.add_argument("--qos", default="rm")
    parser.add_argument("--partition", default="RM-shared")
    parser.add_argument("--time", default="16:00:00")
    args = parser.parse_args()
    treatments = GROUPS.get(args.group, [])
    if args.group == "random-matched":
        print("# Random-matched controls require a generated manifest and are intentionally not auto-submitted.")
        print("# Use fixed-topology scoring for the large branch-length null distribution after validating T0.")
        return
    for treatment in treatments:
        print("sbatch", f"-A {args.account}", f"--qos={args.qos}", f"-p {args.partition}",
              "--cpus-per-task=32", f"--time={args.time}",
              f"--export=ALL,TREATMENT={treatment}", args.script)

if __name__ == "__main__":
    main()
