# Joint genomic tracks

This directory holds reusable nucleotide-statistic estimators and their tests. Dataset-specific scripts, results, and figures stay in each species' own `analysis/joint_genomic_tracks/` directory, alongside that species' existing Patterson's D analysis. It keeps `D_top` (a
frequency difference between alternative gene-tree topologies) distinct from
`D_Patterson` (an ABBA–BABA site statistic).

## f_dM

`fdm.py` implements the frequency-based, modified f_d site numerator and
piecewise denominator used by Dsuite `Dinvestigate`, pools terms before taking
a ratio, orients alleles against a fixed ancestral outgroup, and provides
physical-window aggregation and delete-one-block jackknife utilities. The cod scanner and validator are kept under
`empirical/atlantic_cod/analysis/joint_genomic_tracks/scripts/`; they reuse the
frozen quartet manifest, existing 250-kb tree-window inventory, genotype
parser, and published region intervals. Cod results and figures remain under
`empirical/atlantic_cod/analysis/joint_genomic_tracks/`.

Run the estimator unit tests with:

```sh
python3 -m unittest empirical.scripts.joint_genomic_tracks.tests.test_fdm -v
```

The test suite checks the four piecewise branches against a transcription of
the official Dsuite `D.cpp` implementation, orientation reversal, missing and
undefined sites, ratio-of-sums pooling, window aggregation, and block
jackknifing. This is source-code equivalence validation, not an execution of
the Dsuite binary. The Dsuite CLI comparison remains unrun where the binary is
not installed.

## Directory layout

```text
empirical/scripts/joint_genomic_tracks/           shared estimators and tests
empirical/<dataset>/analysis/joint_genomic_tracks/ dataset scripts, results, figures
```

## Scope and data rules

Outputs are additive. Existing Patterson's D results, frozen designs, and
figures remain authoritative and unchanged. Variant-only VCFs support these
segregating-site ratios under the stated filters, but do not by themselves
provide callable invariant bases; do not use them for unbiased per-base
`d_XY`, `pi`, or `F_ST` without suitable invariant-site/callability data.
Every dataset must pass an input and sample-design audit before an estimator is
run. Unsupported statistics should be reported as unavailable, never imputed.
