# Fire-ant chromosome 16 nucleotide Patterson's D

This directory contains the completed Stage 1 feasibility audit and Stage 2 nucleotide-level Patterson's D analysis. The estimator uses individual VCF calls; TWISST and local-tree quartet frequencies are kept as separate genealogy measurements and never substituted for nucleotide site patterns.

## Reproduce Stage 1

From repository root, run:

```sh
python3 empirical/fire_ants_chr16/analysis/patterson_d/scripts/audit_stage1.py
python3 empirical/fire_ants_chr16/analysis/patterson_d/scripts/audit_stage1.py --scan-vcf
python3 empirical/fire_ants_chr16/analysis/patterson_d/scripts/freeze_quartets.py
python3 empirical/fire_ants_chr16/analysis/patterson_d/scripts/run_stage2.py
python3 empirical/fire_ants_chr16/analysis/patterson_d/scripts/run_stage2.py --transversions-only
python3 empirical/fire_ants_chr16/analysis/patterson_d/scripts/run_stage2.py --strict-calls
python3 empirical/fire_ants_chr16/analysis/patterson_d/scripts/make_figures.py
python3 -m unittest discover -s empirical/fire_ants_chr16/analysis/patterson_d/tests -v
```

The VCF is large and ignored by Git. It was downloaded from the Wurm Lab source linked in `results/nucleotide_source_inventory.tsv`; SHA256 is recorded in the manifest. The `--scan-vcf` pass streams the complete file because no remote or local tabix index was available. It assigns records only by exact BUSCO gene-body interval to one of the existing four-gene windows.

## Estimator and limitations

The estimator requires DP>=10, treats homozygous pseudo-diploid calls as one haploid allele, excludes heterozygotes, and requires at least 50% called per ingroup group. Outgroups must be unanimous among calls (at least four, or all available when a frozen subset has fewer than four); REF-fixed outgroups make ALT derived and ALT-fixed outgroups make REF derived. All contributions are pooled before D is calculated. The unit tests include pure ABBA/BABA, balanced patterns, fixed-ALT polarization, missing/heterozygous calls, and REF/ALT label reversal.

The primary and transversion-only outputs are in `results/`; Figures A–F are in PNG, PDF and SVG under `figures/`. Results are restricted to variants in four BUSCO gene bodies per author-defined window. The genotype release is variant-only, so invariant-site callability, full physical-window coverage, per-read mapping quality, and reference-bias checks are unavailable. `results/ascertainment_audit.md` describes these limits. No commit or push was made.

See `results/data_feasibility_audit.md` and `results/quartet_design_audit.md` for findings and limits.
