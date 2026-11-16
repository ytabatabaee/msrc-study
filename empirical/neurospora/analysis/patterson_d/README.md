# Neurospora nucleotide Patterson's D

This Stage 1 follow-up analyzes the actual deposited Martin & Van Belleghem L10 haploid variable-site genotype matrix with the locally supplied, checksum-validated *N. discreta* alignment. It estimates frequency-weighted nucleotide ABBA/BABA contributions; local-tree topology weights are never substituted for site patterns.

## Inputs and reproducibility

Raw archive/alignment files remain under ignored `data/neurospora/raw/`. The alignment is Dryad file 11939, expected size 63,382,228 bytes and MD5 `95a3eb66abc7d0c44cda5cb7953b2a4f`. The genotype archive is the Zenodo 4946041 dataset, verified at 15,577,951 bytes and MD5 `571324702709eb1c2732a3fd993edc93`.

From repository root:

```sh
python3 empirical/neurospora/analysis/patterson_d/scripts/audit_inputs.py
python3 -m unittest discover -s empirical/neurospora/analysis/patterson_d/tests -v
python3 empirical/neurospora/analysis/patterson_d/scripts/run_patterson_d.py
python3 empirical/neurospora/analysis/patterson_d/scripts/make_figures.py
```

The coordinate check maps 1-based genotype POS to 0-based alignment column after confirming the 2509@6 reference row has no gaps; every source reference allele was compared at the mapped column. Main results are in `results/window_patterson_d.tsv`, `inside_outside_summary.tsv`, `block_uncertainty.tsv`, `spatial_null_validation.tsv`, and `matched_nucleotide_genealogy.tsv`. Figures A–F are generated in PNG, PDF, and SVG.

## Key qualifications

Frozen rooted design: L10 mat-a (10 verified isolates), L10 mat-A (8), *N. hispaniola* 8817@999, rooted by *N. discreta* Ndis@999. CJ07, CJ11, CJ12 and CJ18 are excluded as unresolved. Exact outgroup-bearing quartet topology matching is impossible because the deposited local trees lack *N. discreta*. The genotype matrix contains filtered variable sites, not invariant/callable genome sequence. Consequently the D results are conditional on ascertainment and do not by themselves establish introgression or MSRC.
