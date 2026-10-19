# Numerical validation

`validate_against_msrc_sim.py` compares the symmetric theorem with direct
`msrc-sim` genealogies. It forces every one of the 16 terminal patterns,
uses `p=1/2`, maps `rho=2m), and uses `lambda=1/N). The script accepts
the requested `m) grid, several sister-branch lengths, and multiple root
durations. It writes the all-pattern comparison table and two figures:
`quartet_probabilities_vs_m.png` and `msc_limit_vs_m.png`.

A checked direct run used 50 loci per pattern for
`m={0,0.01,0.1,1}` and sister-branch lengths `{0.1,0.3}`, with
`lambda=0.1) and a long shared root extension. The largest absolute
theory/simulation difference was 0.173 across the small 50-locus cells; the
largest absolute deviation divided by the binomial standard error was 2.46.
The matrix tests are the high-precision validation: the closed two-lineage
formula agrees with `scipy.linalg.expm`, and all conditional probabilities
conserve mass. Larger replicate counts should be used for publication figures.
