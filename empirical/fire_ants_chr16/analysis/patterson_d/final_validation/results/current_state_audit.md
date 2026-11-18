# Current-state audit

## Completed before this package

Stage 1 reconciles 213 author-defined four-BUSCO windows (117 chr1, 42 chr16A,
2 chr16B, 52 supergene) to 852 source gene bodies and VCF scaffold coordinates.
Stage 2 streams 2,650,251 VCF records, assigning biallelic SNPs only to exact
BUSCO intervals, and computes primary, transversion-only, strict-call, geographic,
outgroup-subset, spatial-block, and exact induced-quartet results. Tests cover
the site estimator and key source mappings. Frozen Q1/Q2 memberships were defined
from metadata and the independent chr1–15 tree before D was estimated.

The primary estimator is frequency weighted, pools ABBA/BABA before forming D,
uses DP≥10 and ≥50% call fraction per ingroup, and requires unanimous sampled
outgroup alleles (at least four calls, except smaller frozen subsets). It
supports fixed REF and fixed ALT polarization. Existing transversion, geography,
and outgroup-geography sensitivities preserve the primary direction.

## Key gaps addressed here

The previous strict DP≥15 / 80% comparison is uninformative (Q1 zero, Q2 two
supergene informative SNPs; no useful comparison-region denominator). It cannot
be called a robustness pass. A single additive VCF pass now evaluates the
predeclared 4×4 DP/callability grid for the frozen primary Q1 and Q2. It does not
alter primary outputs or frozen membership. See `filter_attrition.tsv` and
`genotype_filter_sensitivity.tsv`.

The release contains GT:DP only; no GQ, allele depths, or read-level mapping
metrics are supplied. Variant ascertainment is restricted to 852 selected BUSCO
gene bodies, not callable sequence. Scaffold intervals are not continuous
chromosome-16 SNP coordinates. The mapped inversion annotations are already
source-derived and are used as midpoint/window classes only.

Interval audit from the frozen coordinate map finds 37 overlapping gene pairs,
10 across different author windows. The pipeline converts all interval hits to
unique window IDs before pooling, so a SNP in overlapping genes within one
window is counted once per window; a SNP spanning two author windows is assigned
to both, consistent with the existing mapping rule. The latter is a potential
cross-window dependence, not duplicated independent evidence. The 10 cross-window
interval pairs are identified, but their effect on regional pooled totals was
not separately quantified in this package; no claim of zero impact is made.

## Remaining limits

Individual-sample influence, gene/window influence, matched-information
subsampling, and branch-length comparisons were not part of the original run.
The present variant-only file does not permit direct tests of read mapping bias,
allele balance, invariant-site callability, or copy-number/multimapping effects.
In(16)3 has zero source-midpoint-assigned nucleotide windows; its D is
unavailable. In(16)2 has only four assigned windows. Physical-block results are
the defensible uncertainty unit; the circular shift of the 96-window mask
ignores irregular physical spacing and gene density and remains descriptive.

The exact induced-quartet track uses the frozen sampled groups and is already
matched by author window. Its combinations are summarized within each tree and
are not independent replicates. Nucleotide sites and trees derive from related
BUSCO loci, so their concordance is not independent data replication.
