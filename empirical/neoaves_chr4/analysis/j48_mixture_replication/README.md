# J48 mixture replication workflow

This directory prepares and analyzes five independent with-replacement
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
incomplete result. `validate_results.py` and `analyze_replicates.py` are used
after results return. No replication result is present until Bridges-2 output
has been transferred and validated.

The fixed-topology score method is documented in
`results/stage4e/j48_mixture_fixed_topology_score_method.md`. Existing branch
annotations are insufficient for Delta-Q. A constrained ASTRAL output must be
validated against its executable semantics or a direct fixed-candidate scorer
must be run on the returned gene-tree records.

Run the deterministic tests with the configured project Python interpreter:

```bash
python analysis/j48_mixture_replication/tests/test_replication_workflow.py
```

The runbook contains the exact remote paths and SSH/tar transfer commands. It
does not submit jobs automatically.
