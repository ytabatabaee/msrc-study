# Current-state audit

The completed full scan has 988 unique 2L windows, 8 frozen quartet/outgroup designs, and two filters (primary MAC=0 and MAC≥2): 15,808 rows, all finite. Source interval: 2La 20,524,058–42,165,532, 1-based inclusive; source reference: AgamP4 / Ag3.10. Results cover the 2L grid; structural-boundary-overlap windows are retained in tracks but excluded from the primary inside/outside grouping. Existing independent scalar validation covers four selected windows (24 design/filter comparisons), maximum absolute discrepancy 8.53e-14. Existing primary analysis includes block jackknife and spatial-resolution summaries.

Documentation inconsistency: the workflow README still says the full scan has not been submitted, while frozen scan outputs, manifest, checkpoint, and full report record completion. Treat the README sentence as stale; the original provenance remains intact. The README was not overwritten in this package.

The full scan used Ag3.10 at `gs://vo_agam_release_master_us_central1`; no alternative release was substituted. The validation runtime had no importable `malariagen_data`, and this work did not query raw genotypes or download large data.
