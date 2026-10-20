# Stage 13 Phlag cross-system benchmark

Run the benchmark with:

```bash
/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria/bin/python \
  empirical/cross_system_boundary_detection/phlag/scripts/14_phlag_cross_system_benchmark.py
```

The script records the upstream Phlag commit, audits taxon compatibility, copies
the unchanged Atlantic cod Stage 11 raw outputs, and writes validity-aware
diagnostics and comparisons. It does not collapse individual or haplotype tips
into species or population labels. The upstream bundled example is retained in
`raw_results/` as the finite-output software check.
