# Forward persistence validation

The completed local validation contains 20 independent Wright--Fisher
histories for Mai replicate 007 under scenarios A--C. The histories use the
authoritative `msrc-sim` forward process and are shared across switching-rate
comparisons. Results are in `validated_20_tree007_forward.tsv` and the
summarized `validated_20_tree007_persistence.tsv`.

The planned 1,000-history all-tree validation was attempted but stopped before
completion because long forward histories dominate runtime. It is not reported
as complete. The older-root D scenario requires a cluster-scale checkpointed
run before it can support a persistence claim.
