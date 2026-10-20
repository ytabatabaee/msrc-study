# Quartet-definition audit

The primary topology track uses `empirical/house_mouse_t_complex/results/stage2_fixed_quartet_scan.tsv` rows with `treatment=ALL_TIPS`, preserving the frozen four-group quartet definitions and the original `q_species`, `q_t_alt`, and `q_other` support fractions. The historical column `delta_species_alt` is q_species - q_t_alt. Stage 5 writes `D = q_T - q_S = -delta_species_alt`; historical files are unchanged.

The supplementary TTT track uses `TTT_delta` from `empirical/house_mouse_t_complex/results/stage2_arrangement_pattern_spatial_summary.tsv` and applies the same sign transformation. TTT is a status-pattern-specific sampling scheme and is not interchangeable with ALL_TIPS.

Resolved support sums were checked at every one of 4,046 windows. The mean resolved sum was 0.999706794483; the maximum absolute deviation of all support including unresolved from one was 1e-10. Distance values were copied from the frozen Stage 4 table without recomputation.
