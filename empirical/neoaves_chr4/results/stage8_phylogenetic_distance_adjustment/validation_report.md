# Validation report

Parsed 6247 chromosome-4 local ML trees with at least the six exact structural taxa; 0 source records had parser errors. The Stage 7 eligibility rules reproduce candidate regions CR_01, CR_10, CR_11. The Stage 8 residual table contains 514566 pair-region-window rows.

Raw Stage 7 A-distance reproduction is written to `raw_stage7_reproduction.tsv`; it is used as a regression audit and does not overwrite Stage 7. Pair-specific baselines are calculated from windows outside the focal region (broad) or outside all canonical regions (strict).

The outside baseline is estimated from the same chromosome-4 local ML distance series used to form residuals. Therefore outside residuals are centered by construction and inside-versus-outside summaries are descriptive, not independent tests of the baseline model.
