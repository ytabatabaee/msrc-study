# Atlantic cod Stage 8b report

Stage 8b preserves the original study's rooted Atlantic-cod baseline topology and within-cod calibrated branch lengths after outgroup pruning. No root optimization was performed.

## Time-tree audit

- Original induced Atlantic-cod root age: 0.06537923.
- Root-to-tip range: 0.06537923 to 0.06537923.
- Maximum root-to-tip deviation: 1.1999998e-11.
- Root source: original Matschiner et al. baseline tree; outgroups were pruned and no MRCA data selected the root.

## Projection comparison

- Outside projection root age: 0.059036327.
- All-window projection root age: 0.065656219.
- Outside versus original: RMSE 0.010782669, MAE 0.0084991254, correlation 0.92779813.
- All minus outside median absolute internal-node shift: 0.005152591; maximum absolute shift: 0.011771344.
- Largest all-minus-outside clade shifts: {Gadmor_avc_spc+Gadmor_lfo_spc}: 0.01177; {Gadmor_bor_spc+Gadmor_icc_spc+Gadmor_ico_spc+Gadmor_kie_spc+Gadmor_lfc_spc+Gadmor_low_spc}: 0.008073; {Gadmor_avc_spc+Gadmor_bor_spc+Gadmor_icc_spc+Gadmor_ico_spc+Gadmor_kie_spc+Gadmor_lfc_spc+Gadmor_lfo_spc+Gadmor_low_spc}: 0.008073; {Gadmor_avc_spc+Gadmor_avo_spc+Gadmor_bor_spc+Gadmor_icc_spc+Gadmor_ico_spc+Gadmor_kie_spc+Gadmor_lfc_spc+Gadmor_lfo_spc+Gadmor_low_spc}: 0.008073; {Gadmor_avc_spc+Gadmor_avo_spc+Gadmor_bat_spc+Gadmor_bor_spc+Gadmor_icc_spc+Gadmor_ico_spc+Gadmor_kie_spc+Gadmor_lfc_spc+Gadmor_lfo_spc+Gadmor_low_spc+Gadmor_twc_spc+Gadmor_two_spc}: 0.00662.

## Fit diagnostics

| tree version | treatment | SSE | RMSE | MAE | correlation | max absolute residual |
|---|---|---:|---:|---:|---:|---:|
| original_reference | outside_reference_comparison | 0.007673552382 | 0.01078266874 | 0.008499125371 | 0.9277981309 | 0.02968247372 |
| outside_projection | outside | 0.0002349362685 | 0.001886700939 | 0.001303141127 | 0.9895049507 | 0.005764952476 |
| all_projection | all | 0.004444560632 | 0.00820620666 | 0.00680954344 | 0.8588463608 | 0.0164753001 |

The outside projection is the validation of aggregate local MRCA times against the independently calibrated source scale. The all-window comparison tests sensitivity to inversion/supergene-linked windows while holding topology and root fixed. This analysis does not estimate MSRC parameters.

## Comparison with current Stage 8

The fixed original-tree analysis gives the same root ages and the same direction of all-window root expansion as the current ASTRAL-topology Stage 8 result, because the selected root split is concordant with the original baseline arrangement. Internal-node shifts are similar in scale but not identical: Stage 8 reports a median absolute shift of 0.00413352 and maximum shift of 0.012448, whereas Stage 8b reports 0.00515259 and 0.0117713. The fixed-root analysis therefore does not materially change the broad conclusion that including all windows changes the apparent temporal scale, while it provides the independent rooted biological reference requested here.
