# Anopheles 2La joint genealogical and population-genetic tracks

This directory adds 50-kb tracks aligned by the frozen Stage-2R 2L window grid. It retains the published/frozen Patterson's D outputs and writes the allele-state-aware recalculation separately.

## Inputs and populations

Genotypes are read from MalariaGEN Ag3.10, sample set `fontaine-2015-rebuild`, using the configured `gs://vo_agam_release_master_us_central1` source. Frozen homozygous sample IDs come from `data/anopheles_2la/metadata/sample_manifest.tsv` and `data/anopheles_2la/processed/stage1r_primary_homozygote_samples.tsv`. The four arrangement groups are G_s (gambiae standard), G_i (gambiae inverted), C_s (coluzzii standard), and C_i (coluzzii inverted). The quartet design table preserves all eight existing test/outgroup combinations; A_merus and B_merus are the primary analyses and the others are sensitivities.

## Estimators

Ag3 stores genotype values as per-site allele indices: 0 is REF and 1–3 are alternate nucleotides. The frozen Patterson's D helper treated every called value above 1 as multiallelic, which can incorrectly discard a genuinely biallelic site. The new `D_patterson_corrected` and `f_dM` use the count of observed nucleotide states across the quartet at each site and retain sites with at most two observed states. The existing D outputs are preserved and compared window by window in `results/anopheles_2la_corrected_vs_frozen_D.tsv`.

D and f_dM use the frozen quartet populations and no additional `gamb_colu` mask, matching the original D analysis settings. Sites must have no more than 25% missing diploid samples, at least two callable chromosomes in each quartet group, and an outgroup that is complete and fixed. f_dM uses sitewise frequency-based numerator and denominator terms, saved as per-window components and physical-block components; it is not reconstructed from pooled ABBA/BABA totals. MAC >=2 is retained as a sensitivity filter alongside the primary MAC 0 setting.

For pi, dXY and Hudson's F_ST, all Ag3 non-N reference positions are queried, including invariant sites, and the documented `gamb_colu` accessibility mask is applied. Per-group pi at each callable site is the unbiased probability that two sampled chromosomes differ. Pairwise dXY is the per-site cross-population mismatch probability, averaged over sites callable in both groups; each window records callable sites and the sum of effective allele-pair comparisons. Hudson's F_ST is the ratio of summed `dXY - (pi_X + pi_Y)/2` to summed dXY over sites with at least two called chromosomes in both populations. It can be negative under sampling noise. Primary comparisons are G_s/G_i and C_s/C_i; cross-species arrangement comparisons are secondary.

Inside/outside summaries use the original 2La bounds (20,524,058–42,165,532, 1-based inclusive) and exclude boundary-overlap windows. Physical 500-kb block delete-one jackknife intervals account for local linkage/autocorrelation better than treating 50-kb windows as independent. Population arrangement groups have geographic structure; all contrasts remain susceptible to population/geographic confounding and are not by themselves evidence that recombination suppression caused the observed patterns.

## Reproduction

From the repository root, activate `.venv-malaria310`, set `SSL_CERT_FILE` and `REQUESTS_CA_BUNDLE` to `python -m certifi`, then run:

```sh
NUMBA_DISABLE_JIT=1 python empirical/anopheles_2la/analysis/joint_genomic_tracks/scripts/scan_anopheles_tracks.py --chunk-bp 1000000
python empirical/anopheles_2la/analysis/joint_genomic_tracks/scripts/summarize_anopheles_tracks.py
python empirical/anopheles_2la/analysis/joint_genomic_tracks/scripts/plot_anopheles_tracks.py
python empirical/anopheles_2la/analysis/joint_genomic_tracks/scripts/integrate_anopheles_tracks.py
```

The scanner checkpoints by chunk and resumes completed chunks by default. The pilot uses `--pilot`; `--force` reruns the selected scan without reusing its checkpoints.
