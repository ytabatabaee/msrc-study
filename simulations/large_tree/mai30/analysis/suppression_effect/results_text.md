# Results

Within-MSRC suppression effects were small and inconsistent across the four
pilot replicates. For replicate 007, κ=0.5 increased ASTRAL normalized RF
relative to κ=1 by 0.0741, while κ=0.1 had no RF difference. For replicates
042 and 005, all ASTRAL topologies were exactly recovered under every MSRC
treatment. Replicate 074 also had exact recovery under all MSRC treatments.

The extreme-ILS replicate 007 recovered 9 of 27 true internal splits under
κ=1, 7 under κ=0.5, and 9 under κ=0.1. The ordinary MSC control recovered 5
of 27. These are realized finite-sample outcomes. The incorrectly inferred
splits and their identities are retained in `split_recovery.tsv`; the shared
difficult branches should not be interpreted as evidence that MSRC reduces
ILS.

The five-history bootstrap illustrates substantial variation. For replicate
007, bootstrap ASTRAL RF ranges were 0.667--0.889 for κ=1, 0.630--0.815 for
κ=0.5, and 0.630--0.889 for κ=0.1. The ordinary MSC pseudo-block range was
0.741--0.963. These few bootstrap replicates are diagnostic rather than a
formal confidence interval study.

Branch-length error varies substantially even when topology is unchanged.
The ASTRAL table and `branch_length_comparison.tsv` report signed, absolute,
and relative errors only on matching splits. An earlier version used the raw
SimPhy tree and therefore had a factor-1,000 generation-versus-scaled-time
error; the corrected replicate-042 true lengths span approximately
1.88--136.09 CU. Remaining large errors arise mainly from the finite-sample
CULength fit and are present across treatments; they cannot be assigned to
suppression without a larger, independently replicated study.
