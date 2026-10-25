# Patterson's D pilot

The eight frozen Stage-1R windows were queried independently as eight 50-kb Ag3 regions (400,000 bp total). Each design has a primary no-MAC result and a MAC>=2 sensitivity result. This is a feasibility and estimator-validation pilot; it is not an inside/outside inferential analysis.

## Access and mapping

- Windows: 8; designs: 8; rows: 128.
- Sample order was checked against the returned Ag3 `sample_id` coordinate for every query.
- Outgroup polarization required complete diploid calls fixed for one allele; missing or polymorphic outgroups were excluded.
- The primary designs A_merus and B_merus are exploratory population site-pattern contrasts under the frozen baseline, not automatic introgression tests. C/D and melas are sensitivity contrasts with smaller P2 or outgroup groups.

## Primary results

| test | window | D | ABBA | BABA | informative sites |
|---|---|---:|---:|---:|---:|
| A_merus | outside_left_distal_10Mb | -0.143225 | 18.5002 | 24.6855 | 452 |
| A_melas | outside_left_distal_10Mb | -0.0653028 | 22.1011 | 25.1893 | 468 |
| B_merus | outside_left_distal_10Mb | -0.0917372 | 21.5352 | 25.8855 | 372 |
| B_melas | outside_left_distal_10Mb | -0.0252035 | 24.8171 | 26.1004 | 383 |
| C_merus | outside_left_distal_10Mb | 0.122634 | 27.5595 | 21.5384 | 415 |
| C_melas | outside_left_distal_10Mb | 0.136136 | 28.4653 | 21.6437 | 414 |
| D_merus | outside_left_distal_10Mb | 0.161696 | 25.647 | 18.5074 | 458 |
| D_melas | outside_left_distal_10Mb | 0.177704 | 28.3619 | 19.8028 | 466 |
| A_merus | outside_left_flank_50kb | -0.142498 | 25.4736 | 33.9399 | 512 |
| A_melas | outside_left_flank_50kb | -0.0691104 | 31.4275 | 36.094 | 572 |
| B_merus | outside_left_flank_50kb | 0.125793 | 31.4482 | 24.4203 | 433 |
| B_melas | outside_left_flank_50kb | 0.157319 | 39.0048 | 28.4006 | 492 |
| C_merus | outside_left_flank_50kb | -0.187316 | 22.1113 | 32.3041 | 484 |
| C_melas | outside_left_flank_50kb | -0.130833 | 27.4866 | 35.7615 | 542 |
| D_merus | outside_left_flank_50kb | 0.0784969 | 27.4869 | 23.4857 | 550 |
| D_melas | outside_left_flank_50kb | 0.119053 | 34.9501 | 27.5136 | 606 |
| A_merus | inside_left_boundary_50kb | -0.917894 | 4.19407 | 97.968 | 329 |
| A_melas | inside_left_boundary_50kb | -0.722792 | 10.2648 | 63.7938 | 326 |
| B_merus | inside_left_boundary_50kb | 0.821703 | 80.5226 | 7.88106 | 345 |
| B_melas | inside_left_boundary_50kb | 0.929591 | 153.068 | 5.5853 | 456 |
| C_merus | inside_left_boundary_50kb | -0.903725 | 5.09064 | 100.662 | 346 |
| C_melas | inside_left_boundary_50kb | -0.664614 | 13.3169 | 66.0956 | 354 |
| D_merus | inside_left_boundary_50kb | 0.860974 | 80.6067 | 6.02182 | 363 |
| D_melas | inside_left_boundary_50kb | 0.939756 | 155.295 | 4.82304 | 473 |
| A_merus | inside_deep_q1 | -0.789438 | 9.04399 | 76.8593 | 274 |
| A_melas | inside_deep_q1 | -0.632734 | 20.5452 | 91.3365 | 347 |
| B_merus | inside_deep_q1 | 0.812608 | 102.588 | 10.6058 | 463 |
| B_melas | inside_deep_q1 | 0.850239 | 129.576 | 10.488 | 536 |
| C_merus | inside_deep_q1 | -0.720451 | 13.0898 | 80.5598 | 312 |
| C_melas | inside_deep_q1 | -0.700026 | 16.8583 | 95.5401 | 372 |
| D_merus | inside_deep_q1 | 0.843273 | 101.499 | 8.63011 | 486 |
| D_melas | inside_deep_q1 | 0.838155 | 123.637 | 10.886 | 563 |
| A_merus | inside_deep_center | -0.687826 | 12.2101 | 66.016 | 439 |
| A_melas | inside_deep_center | -0.571278 | 21.7217 | 79.6107 | 574 |
| B_merus | inside_deep_center | 0.599154 | 69.9926 | 17.5444 | 469 |
| B_melas | inside_deep_center | 0.607497 | 94.7247 | 23.129 | 627 |
| C_merus | inside_deep_center | -0.61165 | 16.0113 | 66.4465 | 458 |
| C_melas | inside_deep_center | -0.540358 | 24.936 | 83.5658 | 619 |
| D_merus | inside_deep_center | 0.656917 | 75.0149 | 15.5326 | 592 |
| D_melas | inside_deep_center | 0.642465 | 102.754 | 22.3676 | 749 |
| A_merus | inside_deep_q3 | -0.446754 | 18.3192 | 47.9051 | 434 |
| A_melas | inside_deep_q3 | -0.393931 | 21.4171 | 49.2583 | 480 |
| B_merus | inside_deep_q3 | 0.674884 | 76.5851 | 14.8661 | 513 |
| B_melas | inside_deep_q3 | 0.675392 | 87.131 | 16.8816 | 538 |
| C_merus | inside_deep_q3 | -0.433495 | 20.8049 | 52.645 | 506 |
| C_melas | inside_deep_q3 | -0.365615 | 24.3784 | 52.4784 | 559 |
| D_merus | inside_deep_q3 | 0.662836 | 79.0302 | 16.0245 | 576 |
| D_melas | inside_deep_q3 | 0.68668 | 88.3259 | 16.4076 | 601 |
| A_merus | inside_right_boundary_50kb | -0.854451 | 7.74432 | 98.6708 | 414 |
| A_melas | inside_right_boundary_50kb | -0.736907 | 11.8673 | 78.3462 | 431 |
| B_merus | inside_right_boundary_50kb | 0.842762 | 129.193 | 11.0237 | 454 |
| B_melas | inside_right_boundary_50kb | 0.900501 | 171.522 | 8.97989 | 553 |
| C_merus | inside_right_boundary_50kb | -0.810449 | 9.93028 | 94.8463 | 417 |
| C_melas | inside_right_boundary_50kb | -0.585893 | 20.9655 | 80.291 | 452 |
| D_merus | inside_right_boundary_50kb | 0.873059 | 131.231 | 8.89383 | 482 |
| D_melas | inside_right_boundary_50kb | 0.908651 | 184.522 | 8.83128 | 593 |
| A_merus | outside_right_flank_50kb | -0.316894 | 16.7559 | 32.302 | 477 |
| A_melas | outside_right_flank_50kb | -0.243582 | 21.6422 | 35.5807 | 534 |
| B_merus | outside_right_flank_50kb | 0.479806 | 47.2094 | 16.5955 | 358 |
| B_melas | outside_right_flank_50kb | 0.478541 | 57.8882 | 20.4163 | 413 |
| C_merus | outside_right_flank_50kb | -0.0626005 | 30.9763 | 35.1136 | 490 |
| C_melas | outside_right_flank_50kb | 0.0535947 | 39.9975 | 35.9283 | 547 |
| D_merus | outside_right_flank_50kb | 0.585635 | 57.763 | 15.0949 | 454 |
| D_melas | outside_right_flank_50kb | 0.625219 | 73.5735 | 16.9663 | 516 |

## Runtime and exclusions

- Maximum recorded resident memory: 742.4 MiB.
- Sum of independent regional query times: 88.6 s.
- Exclusion counts are recorded as JSON in `pilot_results.tsv`; they are design/filter diagnostics and are not evidence for introgression.

Patterson's D and Stage 2R `D_topology` remain separate estimands. No nonzero value is interpreted as proof of introgression or MSRC causation.
