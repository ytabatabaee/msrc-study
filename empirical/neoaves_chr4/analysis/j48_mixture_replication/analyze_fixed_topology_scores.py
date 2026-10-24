#!/usr/bin/env python3
"""Record the fixed-candidate quartet-score method and pending status."""
from __future__ import annotations
import argparse, json
from pathlib import Path

TEXT = """# Fixed-topology quartet-score method

For fixed candidate trees T0 and NONCHR4, define Q(T) as the sum over the
candidate's induced quartets of the number of input gene-tree records that
resolve each quartet in the same unrooted topology as T. A record contributes
only when all four taxa are present and its induced quartet is resolved;
missing taxa and polytomies are excluded symmetrically. Report Delta-Q =
Q(T0) - Q(NONCHR4), raw eligible quartet counts, and a normalized difference.

This is a fixed-candidate concordance score. It is not ASTRAL localPP and is
not automatically identical to an ASTRAL internal optimization objective.

## Current status

No score is reported yet. The committed outputs contain inferred trees and
branch-local annotations only; they do not contain fixed-candidate scores.
The optional ASTRAL `-C -c` workflow must first be run on a small known test
and its help/output semantics inspected on Bridges-2. If that interface
reports a constrained objective compatible with the definition above, it can
be compared with the direct scorer. Otherwise, use a cluster-side direct
scorer on the returned mixture gene-tree records.

The scorer must keep T0 and NONCHR4 fixed at every epsilon, use the same
records for both candidates, exclude unresolved or incomplete quartets
symmetrically, and report eligible counts with every difference. No smooth
curve or confidence interval is justified before independent inputs are
scored.
"""

def main():
    p = argparse.ArgumentParser(); p.add_argument("--output", type=Path, default=Path("empirical/neoaves_chr4/results/stage4e/j48_mixture_fixed_topology_score_method.md")); args = p.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True); args.output.write_text(TEXT); print(json.dumps({"status": "pending_cluster_scoring", "output": str(args.output)}))

if __name__ == "__main__": main()
