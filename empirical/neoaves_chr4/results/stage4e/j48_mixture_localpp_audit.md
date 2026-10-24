# J48 mixture localPP audit

The raw annotated Newick files were reparsed for every focal branch in all 13 completed mixtures. CU and SU trees were checked against the same frozen taxa and canonical unrooted splits; raw q1/q2/q3 and pp1/pp2/pp3 fields were compared.

## Interpretation correction

The earlier audit incorrectly compared a branch-local q winner with a biological focal-state label. Those labels have different definitions. The corrected table reports the frozen biological focal state separately from the branch-local q winner, branch-local posterior winner, q1/pp1 favored status, near-ties, and low localPP.

For an ASTRAL annotated branch, q1 and pp1 refer to the resolution represented by that inferred branch. At epsilon 0 Columbea, the biological state is q2 but branch-local q1 is largest and localPP is 1; this is locally favored and is not a conflict. At epsilon 0.10 N62, the biological state is q1 but q3 is slightly larger than q1 (approximately 0.340876, 0.313992, 0.345133) and localPP is approximately 2.93e-5. At epsilon 0.10 Columbea, q3 likewise exceeds q1 and localPP is approximately 1.07e-11. These are genuine locally disfavored inferred branches retained in the globally selected 48-taxon topology.

## Interpretation limits

Absent focal splits are missing observations and are not assigned zero support. q1/q2/q3 are branch-local labels and must not be concatenated across different biological splits. localPP is support for an inferred local branch, not the probability that the entire 48-taxon species tree is correct. The completed outputs contain no fixed-candidate Delta-Q trajectory; that requires scoring both predeclared candidate trees against each mixture.

Rows with localPP < 0.001: 12. Machine-readable details: `j48_mixture_localpp_audit.tsv`.
