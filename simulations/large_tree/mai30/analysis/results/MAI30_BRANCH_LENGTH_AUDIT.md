# Mai30 ASTRAL branch-length audit

The large replicate-042 error is real as a numerical estimation issue, but it
does not indicate a failed demographic conversion. The scaled species tree
uses generation lengths divided by `2 * N_scaled`, which equals `g/N_H` under
the validated SimPhy convention. For the 042 MSC control, the 27 true and
inferred internal splits have true lengths from 1.88 to 136.09 CU. The mean
absolute error is 32.74 CU and the signed bias is −32.74 CU.

ASTER/ASTRAL-IV v1.25.4.8 was run with `--length CULength`. Its ASTRAL-IV
source computes quartet based CULength as

\[
 -\log\{1.5-1.5p\},
\]

when the focal quartet frequency `p` is above one third. With 500 gene trees,
very long branches have `p` indistinguishable from one after finite sampling.
The estimator therefore cannot resolve the difference between a long branch
and an even longer branch; its estimates remain finite while the true CU
length can be tens or hundreds of units. This explains the large negative
bias for replicate 042. It is a saturation and finite-sample issue, not
evidence that ASTRAL interpreted the scaled generations as substitution units.

The audit retains every true internal edge. It reports three separate cases in
`analysis/results/mai30_branch_length_edges.tsv`:

- all true internal edges, including splits absent from the inferred topology;
- common true/inferred splits, which measure topology-correct branch error;
- the complete inferred topology comparison in
  `analysis/results/mai30_scaled_wf_astral.tsv`.

No difficult long branches were discarded. Terminal edges and the unrooted
root edge are excluded because the target is internal species-tree splits.
ASTRAL support annotations are retained in its output Newick files, but no
support value is converted into a branch length.

For replicate 042's MSC control, all 27 internal splits are common because the
inferred topology is exact. The branch-length diagnostic is shown in
`analysis/figures/mai30_branch_length_audit.png`; it shows the estimated CU
values flattening relative to the true long branches. Replicate 007's five
common MSC splits have true lengths below 0.042 CU and mean absolute error
0.028 CU, illustrating the separate short-branch regime.

Consequently, the epsilon analysis reports signed and absolute branch errors
with true-length strata and does not use the replicate-042 aggregate MAE alone
as evidence for an MSRC effect. A true-topology evaluation is reported only
when the inferred topology contains that split; it cannot recover an
unobserved branch length from an unrooted ASTRAL output.
