# Production recommendations

Do not launch all 900 production collections yet. The completed pilot shows
that the four selected species-tree replicates are sufficient to verify the
workflow, but not to estimate a systematic suppression effect. Replicate 007
is particularly informative for debugging because all treatments remain
difficult for ASTRAL, while the apparent treatment differences are comparable
to history-resampling variation.

Use a staged, parameter-blind expansion:

1. Select a fixed small set of species-tree replicates from each of the three
   no-outgroup Mai conditions before inspecting suppression outcomes. Include
   low, intermediate, and high ILS representatives using the original Mai
   summaries only.
2. Run all three MSRC treatments for every selected replicate, preserving 25
   histories × 20 loci initially.
3. Run the ordinary MSC controls separately with an explicit statement that
   their pseudo-block bootstrap is not a frequency-history bootstrap.
4. Reassess variance and runtime after the condition-0.15 benchmark on
   Bridges2. If uncertainty is dominated by between-history variation, add
   independent histories before adding more loci per history.
5. Launch the remaining manifest tasks only after the staged results show that
   the resource allocation and estimands are stable.

The existing production manifest and SLURM arrays remain ready for this staged
workflow. Structural or favorable suppression outcomes must not be used to
select the staged replicates.
