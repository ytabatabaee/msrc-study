# Patterson's D full scan

- Completed frozen windows: 988/988.
- Result rows: 15808 (8 designs × 2 filters × 988 windows).
- Valid finite-D rows: 15808; undefined D rows: 0.
- Sum of unique chunk query time: 464.7 s; computation time: 2940.8 s; recorded elapsed time: 3405.5 s.
- Maximum recorded peak RSS: 917.2 MiB.
- Selected-window scalar validation: 24 comparisons; PASS.

## Pooled inside/outside estimates (primary filter)

| design | region | D | informative-site contributions |
|---|---|---:|---:|
| A_merus | inside | -0.73493 | 180425 |
| A_merus | outside | -0.01341 | 200142 |
| A_merus | left_flank | -0.01814 | 48245 |
| A_merus | right_flank | -0.01274 | 53040 |
| A_merus | inside − outside, 500-kb block jackknife | -0.72152 [-0.75258, -0.69046] | 99 blocks |
| B_merus | inside | 0.80885 | 185700 |
| B_merus | outside | -0.00216 | 160736 |
| B_merus | left_flank | -0.00825 | 39442 |
| B_merus | right_flank | 0.02021 | 41174 |
| B_merus | inside − outside, 500-kb block jackknife | 0.81101 [0.78613, 0.83588] | 99 blocks |
| C_merus | inside | -0.70276 | 182297 |
| C_merus | outside | 0.00060 | 188237 |
| C_merus | left_flank | -0.00864 | 46139 |
| C_merus | right_flank | -0.03571 | 50448 |
| C_merus | inside − outside, 500-kb block jackknife | -0.70336 [-0.74344, -0.66329] | 99 blocks |
| D_merus | inside | 0.81868 | 213555 |
| D_merus | outside | 0.01239 | 203504 |
| D_merus | left_flank | -0.00097 | 48003 |
| D_merus | right_flank | -0.00004 | 55839 |
| D_merus | inside − outside, 500-kb block jackknife | 0.80628 [0.77822, 0.83434] | 99 blocks |

Boundary-overlap windows are retained in the full window table and excluded from inside/outside estimands. Block-size sensitivity includes 250 kb, 500 kb and 1 Mb. Spatial aggregation uses pooled ABBA/BABA counts. Full methods, exclusions, topology comparison and interpretation are in the companion TSV/Markdown outputs.
