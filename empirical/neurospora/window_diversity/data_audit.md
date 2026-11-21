# Neurospora L10 — first priority, SNP-panel results complete — data audit

Date: 2026-11-20. Existing MSRC gene-tree coordinates and sample manifests have been audited first. This package now includes clearly labeled SNP-panel analyses where verified genotype calls are available. None of the reported SNP-panel quantities is per callable base: all-sites denominators remain unavailable, so per-callable-nucleotide FST, dXY and pi are not claimed.

## Dataset audit

The local source is `data/neurospora/raw/neurospora/geno/cras_sit_hist_tet_L10.min20MAC1.geno.gz`: 2,265,369 filtered variable rows, 22 haploid CJ columns plus four reference samples. Existing local genealogy comprises 45,303 50-variable-site windows. Directly labeled or pair-resolved sample IDs yield 8 mat-A (CJ01, CJ03, CJ05, CJ09, CJ13, CJ15, CJ19, CJ23) and 10 mat-a (CJ02, CJ04, CJ06, CJ08, CJ10, CJ14, CJ16, CJ20, CJ22, CJ24). CJ07, CJ11, CJ12 and CJ18 remain unresolved and were excluded. This is a verified 8-vs-10 subset, not the archive's full 10-vs-12 design. Focal isolates are haploid homokaryons.

The requested Corcoran et al. Dryad all-sites VCFs and consensus sequences (Dryad 10.5061/dryad.162mh) are not local; the listed file-stream retrieval returned HTTP 403 in this environment. Dryad consensus coordinates use the 2509 mat-a reference, while this calculation uses linkage-group coordinates in the exact local genotype/tree-window inventory; no unverified contig lift was assumed. The SNP-panel analysis retained 44,986 of 45,303 tree windows with at least one site fully called in all 18 samples. It reports Hudson FST and unbiased pi/dXY conditional on complete observed variable SNPs. Inside LGI 1.5–8.0 Mb versus outside LGI, differences were +0.8053 FST (95% 100-kb block interval +0.7184 to +0.9875), +0.2326 dXY (+0.2018 to +0.2742), −0.05774 mat-A pi (−0.06722 to −0.03834), and −0.05894 mat-a pi (−0.06922 to −0.03906), all per observed SNP. This is qualitatively consistent with arrangement divergence but is ascertained and not a per-base result. Neurospora has existing B(w) and q1/q2/q3 tracks but no manuscript-defined arrangement-informed A(w) for these mating-type groups; the matched pair branch contrast is not A(w).

See `empirical/neurospora/window_diversity/results/l10_results_report.md`, `empirical/neurospora/window_diversity/scripts/analyze_l10_snp_panel.py`, per-window files under `empirical/neurospora/window_diversity/results/`, and `empirical/neurospora/window_diversity/figures/l10_snp_panel_tracks.{pdf,png}`.

## Estimators and shared limitations

For each SNP panel, `pi = n/(n-1) * (1 - sum(p_k^2))`, `dXY = 1 - sum(p_Xk*p_Yk)`, and Hudson FST is `sum[dXY - (pi_X + pi_Y)/2] / sum[dXY]` over observed SNPs in a window. Per-window pi and dXY are means over retained polymorphic records. Because sources lack invariant callable sequence, these are conditional per-observed-SNP statistics; SNP-only denominators are never presented as callable genome lengths. Window tracks are joined by the existing native tree-window IDs/coordinates. Suppressed-region boundaries are indicated from each system's existing source annotations. Spatial uncertainty uses physical block resampling rather than treating adjacent windows as independent; contrasts remain descriptive and do not replace a spatial-null hypothesis test.

All reproducible scripts, per-window outputs, joined track-correlation tables, dataset-specific reports and figures are in this directory. Next steps are: (1) obtain/reconcile Dryad all-sites L10 data for a verified 10-vs-12 callable-base analysis; (2) recover the published independent fire-ant 8+8 VCF or callable data and validate the local 78-vs-30 cohort result at read level; (3) reconcile mouse gVCF reference blocks across the frozen samples, compute common per-window callable lengths, and repeat geographic sensitivity analyses. Do not describe the present SNP-panel values as absolute nucleotide-diversity estimates.

## Next step

Recover the Dryad all-sites VCF and consensus files, validate the 2509-reference coordinate lift, and replace SNP-panel denominators only after callable positions are reconciled.
