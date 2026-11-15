# Chr17 Patterson D ancestral-state validation

This is an additive validation of the frozen 55-sample, three-design analysis in `../`. It preserves the original files and applies fixed-outgroup polarization from the actual nucleotide calls: all eight *M. spretus* genotypes must be complete and fixed for either REF or ALT. Fixed-ALT sites use derived REF frequencies (`1 - ALT frequency`). Polymorphic/missing outgroups remain excluded. The scan reads the indexed AllMouse PASS VCF over HTTP range requests and does not download the full VCF.

## Run

From the repository root, using `.venv-malaria310`:

```sh
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/polarization_validation/scripts/run_analysis.py
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/polarization_validation/tests/test_polarization.py
```

The scan checkpoints compact 5-kb polarization-class contributions after each 500-kb indexed fetch in `results/site_checkpoint.json`. Re-running resumes completed chunks. The corrected analysis writes only into this directory.

`results/corrected_window_D.tsv` and `results/corrected_multiresolution_D.tsv` contain pooled ABBA/BABA counts and D. `results/polarization_class_contributions.tsv` retains fixed-REF, fixed-ALT, polymorphic-outgroup, and missing-outgroup/missing-population counts. `results/original_vs_corrected_D.tsv` preserves paired original and corrected values. `results/polarization_audit.md` gives the numerical audit and inversion-union summaries.

Outgroup and ingroup calls are unphased diploid VCF genotypes. The statistic is population allele-frequency asymmetry and does not isolate the t-bearing chromosome. See `../README.md` and the exact-quartet comparison report for shared limitations.


## Optional robustness scan status

`results/robustness_checkpoint.json` records the optional transversion-only and complete-call-rate sensitivity scan over the matched 5–37.055 Mb region. It is partial (59 of 65 indexed chunks completed) because the remote indexed VCF endpoint repeatedly returned a truncated BGZF block at a high file offset. No complete robustness result tables are produced; the primary fixed-outgroup corrected scan is complete.
