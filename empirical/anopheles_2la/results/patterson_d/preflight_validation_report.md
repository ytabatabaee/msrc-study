# Patterson's D final preflight validation

Status: **READY_FOR_FULL_SCAN**. The 988-window scan was not run.

## Numerical correctness

The retained scalar estimator and the production NumPy implementation were
compared on a deterministic 500-site genotype fixture containing missing,
heterozygous, multiallelic, fixed, and intermediate-frequency calls. ABBA and
BABA sums agreed within `1e-8` and D agreed within `1e-10`, with identical
exclusion counts. The independent algebraic fixture in
`independent_estimator_validation.tsv` covers the limiting values -1, 0, and
1, allele-frequency weighting, missing calls, and REF/ALT complementation.

Dsuite is not installed in the available environments, so no Dsuite comparison
is claimed. The validation is against an independently calculated reference
implementation of the stated estimator plus the retained scalar code. Dsuite
or another VCF-level tool can be added as a separate cross-check later.

## Real-data semantics

The completed pilot used the authenticated Ag3.10 `snp_calls` array and
verified the returned sample IDs and order for every query. Each frozen window
was requested as a 1-based inclusive 50,000-bp interval, and the returned
position coordinate was checked for uniqueness. The pilot therefore examined
50,000 coordinate-aligned records per window; invariant records have zero
allele-frequency contributions and were not silently converted into gene-tree
or quartet observations. The retained pilot TSV does not preserve the raw
genotype arrays, so exact invariant-versus-polymorphic counts cannot be
reconstructed from the aggregate output alone. The full scan should retain
those optional diagnostic counts if they are needed for the manuscript.

Outgroup polarization used complete diploid calls fixed for one allele. REF
was never assumed to be ancestral. Multiallelic calls, unsuitable missingness,
missing or polymorphic outgroups, and sensitivity MAC failures were excluded
with explicit counts. The eight-window pilot examined 3,199,224 primary
positions and had 30,283 informative positions; exclusions were 503,680
multiallelic, 349,315 site-missing, 228,489 outgroup-missing, 42,039
outgroup-polymorphic, and 7,145 insufficient-call sites. No duplicate
positions were returned. A site fixed for the derived allele in ingroups is
handled by the allele-frequency estimator and requires no additional filter.

## Performance and memory

Vectorization reduced the synthetic 50-kb-equivalent computation from 9.863339
to 0.453567 seconds (21.7x) while preserving the sums and D. The vectorized
fixture used about 96.10 MiB peak RSS versus 72.35 MiB for the scalar
subprocess. This is a bounded synthetic measurement, not a promise for remote
query time. A simple CPU-only extrapolation is approximately two hours for
988 windows × 8 designs × 2 filters at the measured independent-window rate;
the production code batches contiguous windows, so network access must be
measured from the first full chunk separately.

The historical pilot rows retain the original `peak_rss_mb` field. On macOS
that raw value was KiB, so the corrected interpretation is 742.4 MiB. New
runs use the platform-aware `peak_rss_mib` field. Historical rows were not
rewritten.

## Checkpoint and statistical audit

Checkpoint validation now rejects unknown, duplicate, truncated, coordinate-
mismatched, and design-mismatched rows; writes TSV and JSON checkpoints
atomically; skips a query when all rows for its region are complete; and
validates the full output table before declaring completion. The full scan
remains restartable and does not reset a completed pilot to design-prepared.

Jackknife calculations use only informative physical blocks and apply the same
block deletion to inside and outside contributions for the difference. Results
with too few blocks or zero-denominator deletion replicates return missing
uncertainty rather than fabricated SEs. Boundary-overlap windows are excluded
from the primary inside/outside summary; a physical block may still contain
both region classes, with region-specific sums retained.

## Execution gate

The existing SLURM wrapper requests 8 CPUs and 14 GB, or 1.75 GB per core,
below the 2,000 MB/core RM-shared limit. It creates logs when executed from
the documented working directory and runs only with an explicit `--run-full`.
No job was submitted. The next execution should be preceded by an authenticated
first-chunk preflight and monitoring of remote data-access time.
