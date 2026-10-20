# Anopheles Stage 12 spatial distance analysis

Processed 978 Stage-4R local NJ trees. Branch lengths are patristic distances in the published NJ substitution-distance scale; no molecular clock or coalescence-time interpretation is made.

The primary arrangement contrast averages within-species 2La/2L+a pair means and equalizes species before forming `A_distance = opposite - same`. Heterokaryotypes are excluded from that contrast because haplotype phase is not established, but all mapped tips contribute to arrangement-blind summaries. The individual-pair weighted sensitivity is retained per window.

The frozen 2La interval is 20.524058--42.165532 Mb. Boundary-overlap windows remain in the track and are not silently interpolated. The observed inside/outside contrast spatial-null result is recorded in `stage12_spatial_null.tsv`; it is an alignment test against the prespecified interval, not a test that a change point exists.

A limitation is that the NJ trees are local substitution-distance summaries and may be affected by informative-site density. The per-window post-filter SNP count is included for scale diagnostics.
