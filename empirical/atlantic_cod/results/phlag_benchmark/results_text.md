# Results

Phlag ran successfully on all four linkage groups. It analyzed 9 focal internal branches per linkage group and produced 36 raw output files. The basal three-taxon clade N009 was excluded by the upstream reader's basal-bifurcation collapse.

Under the default HMM state call, no eligible focal branch flagged a window in any linkage group. The Phlag union therefore had precision, recall, F1, and interval IoU equal to zero for LG01, LG02, LG07, and LG12. The branch-specific best and median summaries are in `phlag_branch_summary.tsv`; they are also zero.

The fixed comparator rules produced the following pooled window-level results: topology-only F1=0.7195 and mean LG IoU=0.4270; time-only F1=0.7516 and mean LG IoU=0.4322; joint F1=0.7939 and mean LG IoU=0.4402. LG07 was detected by the time-only sign rule with F1=0.5324 and IoU=0.3035, while Phlag had F1=0 and IoU=0.

These numerical comparisons do not establish superiority of the comparator detectors over Phlag because the SNAPP MCC tree input violates the ordinary gene-tree interpretation assumed by Phlag. They show that the original Phlag run completed and was highly conservative on this input representation.
