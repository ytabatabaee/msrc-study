# MSRC pairwise coalescence and Anopheles 2La: initial results

Date: 2026-11-21. The mathematical implementation, numerical tests, and
first callable-base empirical comparison are complete. Other empirical
systems were screened for denominator validity; they do not yet permit
absolute per-base prediction checks.

## Pair-time calculation and validation

The exact three-state matrix-exponential operator is documented in
`../derivation.md` and implemented in `../scripts/pairwise_msrc.py`. It
preserves the legacy `msrc-sim` frequency-dependent coalescence hazards,
effective background switching, finite-copy frequency floor, finite structured
origin, and ordinary MSC after the origin. Eight focused unit tests passed,
including comparisons of exact expectations with 2,500 direct event-time
genealogies per initial pair class at fixed frequency. A separate stochastic
validation used 100 independent forward histories in each of 12 Moran/Wright–
Fisher, neutral/selected, and initial-frequency cells. The history-paired
analytic-versus-one-genealogy difference interval included zero in 35/36
pair-state cells. Small effective population sizes make those stochastic
validation intervals intentionally broad; this demonstrates implementation
agreement rather than biological calibration.

For the tested fixed-frequency grid (B frequency 0.05–0.95, normalized
switching \(2N_{e,in}\rho=0–100\), duration
\(\tau/(2N_{e,in})=0.1–8\), and \(N_{e,in}/N_{e,out}=0.5–10\)), R was
monotone decreasing with increasing switching in every p-duration-Ne-ratio
combination examined. Equivalently, reducing switching increased the relative
divergence ratio in this grid. The result does not establish universal
monotonicity under arbitrary histories, demography, or linked genealogies.
Duration and p strongly change the magnitude. See
`../results/theoretical_parameter_grid.tsv` and
`../figures/pairwise_msrc_parameter_sensitivity.pdf`.

## Anopheles callable-base comparison

We reused the frozen 50-kb AgamP4 2L table at
`empirical/anopheles_2la/analysis/joint_genomic_tracks/results/anopheles_2la_dxy_fst_windows.tsv`.
For both primary within-species comparisons, dXY and Hudson FST have identical
callable-site counts, and πX/πY were calculated over this same pairwise mask.
We form each window's \(\Delta=d_{XY}-(\pi_X+\pi_Y)/2\) and
\(R=2d_{XY}/(\pi_X+\pi_Y)\) from that matched table; we do not use the
separately generated standalone π tracks with different callability. All
statistics include callable invariant sites. Boundary-overlap windows were
excluded. The fixed region is 2L:20,524,058–42,165,532 (1-based inclusive);
the windows use 50-kb bins and uncertainty resamples physical 500-kb blocks.

| Species comparison | Region | πX | πY | dXY | Δ | R | pooled Hudson FST |
|---|---:|---:|---:|---:|---:|---:|---:|
| Gambiae G_s vs G_i | outside | 0.014163 | 0.013313 | 0.014101 | 0.000363 (0.000302–0.000423) | 1.0264 (1.0216–1.0316) | 0.0257 |
| Gambiae G_s vs G_i | inside | 0.015468 | 0.015006 | 0.024177 | 0.008940 (0.008318–0.009587) | 1.5867 (1.5352–1.6414) | 0.3698 |
| Coluzzii C_s vs C_i | outside | 0.012235 | 0.013148 | 0.013410 | 0.000719 (0.000662–0.000776) | 1.0566 (1.0522–1.0618) | 0.0536 |
| Coluzzii C_s vs C_i | inside | 0.012811 | 0.014335 | 0.024009 | 0.010435 (0.009833–0.011081) | 1.7688 (1.7074–1.8357) | 0.4346 |

Parentheses are 95% physical-block bootstrap intervals for Δ and R. Independent
inside/outside block bootstrap gives Δ contrasts of +0.00858
(0.00795–0.00922) for Gambiae and +0.00972 (0.00911–0.01037) for Coluzzii;
R contrasts are +0.561 (0.508–0.614) and +0.713 (0.651–0.779), respectively.
The initial regional estimates supplied in the request are reproduced.
Within-background diversity also rises inside 2La relative to outside: Gambiae
πX +9%, πY +13%; Coluzzii πX +5%, πY +9%. dXY rises by 71% and 79%,
respectively. Thus the data show the requested joint pattern: sharply higher
between-background divergence without a corresponding loss of both
within-background diversities.

At the 50-kb scale, within-region matched Δ and R correlate with frozen
arrangement-distance A(w): Spearman ρ=0.881/0.929 for Gambiae Δ/R inside 2La,
and 0.820/0.886 for Coluzzii. The 500-kb block intervals exclude zero. By
contrast, inside-2La correlations with q_arrangement are weak (ρ≈0.05/0.02 for
Gambiae and 0.13/0.11 for Coluzzii; block intervals include zero). Thus the
large correlation obtained by pooling the whole chromosome arm would largely
reflect the common inversion boundary and must not be read as a within-region
genealogy calibration. Correlations with mean tree pair distance inside are
modest (Δ/R ρ≈0.27/0.35 Gambiae, 0.21/0.29 Coluzzii). The detailed
inside/outside block-bootstrap estimates are in
`anopheles_genealogy_track_correlations.tsv`. These are exploratory
correspondences: A(w) is NJ patristic distance in substitution units, and the
local tree metrics are not calibrated pairwise coalescence times. They cannot
independently calibrate the model's time scale.

## Interpretation, identifiability, and limits

The legacy MSRC model can generate excess between-background pair time while
retaining or increasing within-background diversity for suitable frequencies,
switching, and duration; the fixed-history curves demonstrate existence, not
a fitted explanation of 2La. The empirical target alone does not identify
\(N_e,\rho,\tau,p(t),\mu\): absolute statistics confound mutation rate and
time scale, while R removes that scale but is only one ratio of three
dependent diversity quantities. Arrangement-stratified sample counts are not
natural p(t). Centimorgans/Mb are not converted to the simulator's effective
switching rate. Outside 2La is not assumed to be a single panmictic
population; geography, species structure, and introgression remain competing
explanations for nonzero baseline Δ.

The local gene-tree signal is concordant spatially with the pairwise distance
signal but is not an independent estimate of coalescence time. A joint
likelihood would require a calibrated genealogy time observable and a
geography/introgression demographic null. Without those, we cannot conclude
that one empirically identified parameter set explains all three signals.

As a transparent sensitivity check, we compared the three inside statistics
to an equal-weight average of outside \(\pi_X,\pi_Y,d_{XY}\), then ranked the
grid by unweighted squared distance in those three normalized moments. The
best grid point for Gambiae has \(N_{e,in}/N_{e,out}=2,p=0.5,
2N_{e,in}\rho=1,\tau/(2N_{e,in})=0.5\); predicted ratios
\((\pi_A,\pi_B,d_{AB})\) are (1.137, 1.137, 1.810), versus observed (1.116,
1.083, 1.744). For Coluzzii the best point has ratio 2, p=0.5, normalized
switching 0.3, and duration 0.5; predictions (1.050,1.050,1.931) are close
to observed (0.991, 1.109, 1.857) but do not reproduce the within-group
asymmetry. The top 25 candidates span multiple p, switching, duration, and Ne
ratios in both species. This is a grid-ranking diagnostic, not a likelihood,
confidence region, or parameter estimate: equal weights ignore their spatial
sampling covariance, the outside average is not a fitted demographic null,
and no genealogy-time term is included. The current analysis therefore shows
that a common-order parameter regime can reproduce the simultaneous scale of
within diversity and excess dXY, but does not validate a unique shared
biological parameterization or explain the calibrated local-genealogy signal.
See `anopheles_grid_moment_fit_top25.tsv` and
`../figures/anopheles_observed_vs_grid.pdf`.

Fire ants, Neurospora L10, and house mouse currently provide SNP-conditioned
diversity panels without validated callable invariant-site denominators. For
these, absolute π and dXY comparisons remain blocked; absent sites are not
called invariant. Atlantic cod's window divergence tracks are time-scaled
SNAPP population-tree divergences, not pairwise gene-lineage coalescence times,
and the available nucleotide summaries are SNP-ascertained. They are assessed
separately and no common six-dataset parameterization is forced.

## Outputs and next steps

- Derivation: `../derivation.md`.
- Reusable implementation/tests: `../scripts/` and `../tests/`.
- Raw history-level numerical validation: `frequency_history_validation.tsv`;
  summary and manifest are adjacent.
- Per-window matched Anopheles tracks and regional/block outputs: files named
  `anopheles_*` in this directory.
- Figures: `../figures/anopheles_matched_statistics.pdf` and
  `../figures/pairwise_msrc_parameter_sensitivity.pdf` (PNG counterparts also
  provided).

Next steps: conduct a model adequacy/sensitivity fit against spatially
stratified Anopheles targets using explicit population history and geography
priors; calibrate or replace local NJ branch summaries with matched sequence
genealogies in time units; assess geographic/introgression alternatives; and
only add absolute fire-ant, mouse, or Neurospora comparisons after callable
invariant-site masks pass their existing source audits. This initial work does
not modify the frozen empirical source analyses.
