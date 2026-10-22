# Full 363-Taxon Neoaves T0 ASTRAL-IV Validation

Generated UTC: 2026-10-22T16:10:39.447939+00:00

This report validates committed output only. No ASTRAL-IV subprocess was invoked.

## Provenance

- Recorded ASTRAL-IV host: `r305.ib.bridges2.psc.edu`; job `49555313`; CPUs `32`.
- Recorded input: `/ocean/projects/cis260126p/shared/neoaves_stage4e/repo/empirical/neoaves_chr4/data/stage4e/T0.tre`; SHA-256 `c0a0ec9b432405da548a6cdb37460ccf1c0d9d1b5a5211d74d5f397c70e0f943`.
- ASTRAL binary SHA-256: `38d2aad2879605b3bd3d26bb7dbc978dce9cec510d684f64da4a58ead836872a`; `-u 2`, `--genelength 1000.0`, `--sulengthtype 1`.
- Both original runs returned code 0; the stderr logs report optimization score `30019936712278` for CU and SU.
- The committed run metadata does not include process end timestamps, so a walltime cannot be reconstructed from repository files alone.

## Tree and annotation audit

- Taxa: **363** in each original and extracted tree; duplicate labels: none.
- Unique nontrivial unrooted internal splits: **360** in each tree (expected `363 - 3 = 360`).
- Original CU/SU unrooted RF symmetric difference: **0**; shared splits: **360**.
- Original CU tree: 360 internal CU estimates; terminal CU estimates: 0.
- Original SU tree: 360 internal SU estimates and 363 unique terminal SU estimates.
- Extracted tree: 360 internal CU/SU estimates, 363 unique terminal SU estimates, and 0 terminal CU estimates.
- Maximum extracted-vs-original CU difference: `0.0`.
- Maximum extracted-vs-original SU difference: `0.0`.
- Maximum CU annotation difference between independent original runs: `0.0`.
- Maximum SU annotation difference between independent original runs: `0.0`.
- Maximum support/q annotation differences between original runs: `0.0`.

## Internal-length summaries

- CU range: `0.0` to `9.05635`; median `2.140465`.
- SU range: `0.000112919` to `0.0790686`; median `0.00514644`.
- All localPP values finite/nonnegative: `True`.
- The CU and SU trees have the same topology; only the edge-length representation differs.

## Focal relationships

Focal roles and taxon sets come from the frozen Stage 3B role mapping and `focal_split_states()`; no approximate branch was selected.
- An independent reference species tree for focal agreement was not available in the repository. The existing `empirical/neoaves_chr4/results/stage4d/trees/FULL363_ORIGINAL_REPRO.nwk` path is empty, and the large `.tre` file is a gene-tree collection rather than a reference species tree; no focal agreement claim is made.

| Relationship | T0 state | Reference state | CU | SU | localPP | q1 | q2 | q3 | exact reference split | |
|---|---|---|---:|---:|---:|---:|---:|---:|---|
| Columbea | q2 | unavailable | 0.0567416 | 0.00198804 | 1.0 | 0.370114 | 0.302169 | 0.327717 | unavailable |
| N61 | q1 | unavailable | 0.0256619 | 0.000368452 | 1.0 | 0.350229 | 0.343002 | 0.306769 | unavailable |
| N62 | q1 | unavailable | 0.0492232 | 0.00307844 | 1.0 | 0.36536 | 0.328157 | 0.306483 | unavailable |

## PASS/FAIL checks

| Check | Result |
|---|---|
| taxa_363 | PASS |
| unique_internal_splits_360 | PASS |
| original_cu_su_unrooted_rf_zero | PASS |
| source_cu_lengths_finite_nonnegative | PASS |
| source_su_lengths_finite_nonnegative | PASS |
| extracted_cu_terminal_lengths_missing | PASS |
| extracted_su_terminal_lengths_present | PASS |
| extracted_source_hash_matches | PASS |
| extracted_cu_topology_matches_original | PASS |
| extracted_su_topology_matches_original | PASS |
| focal_mappings_exact_in_t0 | PASS |

## Interpretation and scope

T0 establishes the full-taxon reference topology and branch-length annotations for later treatment comparisons. It does not by itself establish a chromosome-4 structural effect or causal rearrangement mechanism. Future analyses must compare stable unrooted splits and distinguish topology changes from CU and SU length changes on matching splits.
