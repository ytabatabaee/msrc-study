# Atlantic cod f_dM — Phase 1 component

This additive analysis estimates Dsuite-style `f_dM` in the four frozen
Atlantic cod quartet designs (LG01, LG02, LG07, LG12), using exact frozen
sample IDs and the existing 250-kb tree-window inventory. The four in-region
SNAPP VCFs and shared outside-SNAPP VCF are streamed from the already present
local downloads; no download or alteration of the original Patterson's D
results is required. SNAPP VCFs are used because they contain the frozen
individuals; the Dsuite VCF is not a drop-in source for every exact frozen
quartet.

Run from the repository root:

```sh
python3 empirical/atlantic_cod/analysis/joint_genomic_tracks/scripts/scan_cod_fdm.py
python3 empirical/atlantic_cod/analysis/joint_genomic_tracks/scripts/validate_cod_fdm.py
```

The scan filters PASS (or `.`) biallelic single-base records, requires complete
ingroup calls, and retains sites with a fixed called outgroup. It polarizes
using that outgroup state. Region summaries pool site numerators and
denominators, with 500-kb physical-block delete-one jackknife intervals.
The window-level Patterson D and f_dM are shown in separate figures, each with one statistic and one y-axis per panel: `figures/atlantic_cod_patterson_D_tracks` and `figures/atlantic_cod_fdm_tracks` (PNG/PDF/SVG). Observed windows are shown as points; gaps are left unconnected. The regional comparison is in `figures/atlantic_cod_fdm_regions` (PNG/PDF/SVG).

## Results

The inside estimates differ in sign among linkage groups: LG01 +0.423 (95%
block-jackknife interval 0.393–0.452), LG02 +0.322 (0.266–0.377), LG07 −0.357
(−0.377 to −0.337), and LG12 +0.369 (0.331–0.407). The corresponding outside
estimates are near zero and all intervals include zero: LG01 −0.024 (−0.087 to
0.038), LG02 −0.015 (−0.043 to 0.013), LG07 −0.007 (−0.045 to 0.031), and
LG12 +0.023 (−0.032 to 0.077). The `f_dM` signs are conditional on the frozen
P1/P2/P3 ordering and do not alone establish introgression or its direction.
These ratios are not literal admixture proportions.

A 421-row cross-check against the existing window-level `D_Patterson` values
passes, with maximum absolute difference 4.97e-11. This verifies that the
new scanner reproduces the old ABBA–BABA statistic on matched windows under
its filters. It does not independently validate f_dM against a running Dsuite
binary. Source hashes and validation status are recorded in
`results/atlantic_cod_fdm_provenance.json`.

Because the input VCFs are variant-only, this component does not report
per-base `d_XY`, nucleotide diversity, or Hudson's `F_ST`. Those need callable
invariant-site denominators or equivalent validated callability data.
