# J48 mixture localPP audit

The raw annotated Newick files were reparsed for every focal branch in all 13 completed mixtures. CU and SU trees were checked against the same frozen taxa and canonical unrooted splits; raw q1/q2/q3 and pp1/pp2/pp3 fields were compared.

## Finding

The unexpectedly small localPP values are annotation-consistent and are not evidence of a failed ASTRAL run. ASTRAL's localPP and q1/q2/q3 describe the quartet induced by the inferred branch. The globally selected 48-taxon topology can retain a branch whose local quartet posterior is weak when the full topology optimization trades off several branches. At epsilon 0.10, for example, the N62 branch is in the T0 tree while its raw q frequencies are approximately (0.340876, 0.313992, 0.345133), so q3 is slightly larger than q1; its localPP is approximately 2.93e-5. This is a local conflict within an otherwise valid globally optimized tree, not a fixed-topology score.

## Interpretation limits

Absent focal splits are missing observations and are not assigned zero support. q1/q2/q3 are branch-local labels and must not be concatenated across different inferred splits. localPP is support for an inferred local branch, not the probability that the entire 48-taxon species tree is correct. The completed outputs contain no fixed-candidate Delta-Q trajectory; that requires scoring both predeclared candidate trees against each mixture or validating the constrained ASTRAL interface.

Rows with localPP < 0.001: 12. Machine-readable details: `j48_mixture_localpp_audit.tsv`.
