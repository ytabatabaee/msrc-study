# Validated Mai30 quartet-reversal experiment

This run extends the earlier 20-history pilot with the optimized authoritative
`msrc-sim` implementation at commit
`bfba58d627f6f467af649cdb0b37f969349e0aab`. It uses the four frozen Mai
species trees, scenarios A (the existing internal origin) and B (the
deterministic deeper internal origin), and the prespecified switching grid
`0, 0.0001, 0.001, 0.002, 0.01, 0.02`.

Each forward Wright–Fisher history is simulated once and reused across all
six switching values. One genealogy is generated per history, switching value,
and selected quartet. The independent-history unit is therefore the history,
not the individual quartet row. Checkpoints are atomic compressed JSON files
under `validated_100/`; they include the frequency records, sampled terminal
arrangements, event counts, seeds, and simulator commit.

The deterministic panel contains 12 quartets per tree/scenario, stratified by
structural exposure and induced internal-branch length. It is selected from
species-tree geometry and forward-history design metadata only; observed gene
tree topologies are not used for selection. The full 30-taxon collection and
ASTRAL analysis are intentionally outside this pilot.

The target is 100 histories for each of the eight configurations. The current
run is partial: scenarios A and B for replicate 007 are complete; A for 042,
074, and 005 is incomplete; B for those three trees has not been run. The
current counts are recorded in `validated_100/forward_histories_manifest.tsv`.

Resume or run a disjoint configuration with:

```bash
/Users/ytabatabaee/Desktop/msrc-sim/.venv/bin/python \
  simulations/large_tree/mai30/analysis/origin_placement/run_validated_experiment.py \
  --histories 100 --tree-indices 0 1 2 3 \
  --scenarios A_existing_internal B_deep_internal \
  --output-dir simulations/large_tree/mai30/analysis/origin_placement/validated_100
```

The 1,000-history expansion remains gated until all initial configurations
have enough complete histories to support an independent validation decision.
Raw checkpoints are retained locally for restart but excluded from the commit.
