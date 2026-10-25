# J48 mixture replication workflow

This directory analyzes five independent with-replacement
replicates at each of epsilon 0.02, 0.05, 0.075, 0.10, and 0.15. The 25
inputs contain exactly 62,945 records and use deterministic seeds 76001–76405
in five blocks. The frozen manifest records source, locus-manifest, and taxon
mapping checksums, requested and realized weights, and record counts.

`audit_localpp.py` audits the 13 completed mixture outputs. It confirms that
very small localPP values are valid branch-local annotations; they do not
represent whole-tree probabilities or fixed-candidate quartet scores.

`prepare_manifest.py` creates the manifest locally and realizes one input on a
cluster node with `--realize`. `cluster/run_replicates.sbatch` runs one
single-call `topology_su` ASTRAL-IV inference per job and extracts both CU and
SU annotations. It is restart-safe and refuses to overwrite an existing
incomplete result. All 25 results are now present under
`results/stage4e/J48_MIX_REPLICATES/`. `validate_results.py` verifies the
returned runs, and `analyze_replicates.py` produces replicate-level tables,
Wilson intervals, focal-state summaries, branch-length tables, and figures.
Postprocessing never rewrites ASTRAL outputs.

The fixed-topology score method is documented in
`results/stage4e/j48_mixture_fixed_topology_score_method.md`. Existing branch
annotations are insufficient for Delta-Q. The executable direct scorer defines
Q on the fixed candidate-disagreement quartet universe and has a brute-force
reference implementation for small fixtures; it has not been run on the full
mixtures. A constrained ASTRAL output must still be validated against its
executable semantics before being interpreted as the same score.

Run the deterministic tests with the configured project Python interpreter:

```bash
python analysis/j48_mixture_replication/tests/test_replication_workflow.py
```

The runbook contains the exact remote paths and SSH/tar transfer commands. It
does not submit jobs automatically.
