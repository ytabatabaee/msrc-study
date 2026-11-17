# House-mouse chr17 Patterson D final validation

This additive finalization preserves the frozen pipeline and prior results. `results/` distinguishes corrected fixed-REF+ALT D, superseded fixed-REF-only values, exact q3−q2 topology, full-contig vs physical-span vs exact-tree-matched contrasts, completed filter sensitivity, and incomplete sample-level sensitivity. Quartet and phasing audits are in `../results/final_validation/`; the concise current-state audit is `results/current_state_audit.md`.

Reproduction from repository root (environment: `.venv-malaria310`):

```sh
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/polarization_validation/tests/test_polarization.py
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/exact_quartet_comparison/tests/test_quartet_scoring.py
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/final_validation/tests/test_scalar_reference.py
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/final_validation/scripts/build_final_validation.py
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/final_validation/scripts/summarize_robustness.py
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/final_validation/scripts/summarize_exact_tree_robustness.py
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/final_validation/scripts/block_uncertainty_validation.py
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/final_validation/scripts/summarize_spatial.py
```

Primary corrected scan reproduction is documented in `../README.md`; exact quartet reproduction is in `../exact_quartet_comparison/README.md`. Optional robustness retry command: `.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/polarization_validation/scripts/robustness_scan.py`; it resumes complete chunk checkpoints; this run verified 65/65 chunks and generated the transversion and 80/90/100% call-rate results. No full 55.7-GB VCF download was performed. The separate sample-level scan was attempted but the endpoint truncated BGZF at byte offset 48,567,409,620 before outputs completed; no sample-level estimates are reported.

Paper readiness: **CONDITIONAL GO** for a carefully framed population-level allele-sharing analysis with the filter sensitivity scan complete and sample-level influence flagged as unresolved; **not** ready as phased t-haplotype evidence.
