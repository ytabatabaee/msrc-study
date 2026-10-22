# Mai30 MSRC–MSC ASTRAL suppression-effect report

## Findings

The completed pilot does not show systematic evidence that stronger
recombination suppression changes ASTRAL topology accuracy. The within-MSRC
contrast is zero for many collections and changes direction in the difficult
replicate. In replicate 007, moderate suppression is worse than κ=1 by one
normalized-RF step (`+0.0741`), while strong suppression matches κ=1.

κ=1 can outperform ordinary MSC because it is not ordinary MSC: finite
background switching, arrangement inheritance, frequency-dependent
coalescence, and the finite structured interval remain active. The ordinary
MSC control has none of these processes. The comparison therefore identifies
the effect of replacing MSC with the MSRC hierarchy, not a pure suppression
effect.

Branch-length estimation is sensitive even when topology is identical. This
variation is consistent with finite quartet-frequency sampling, CULength
estimation, and long or poorly supported branches. The pilot cannot separate
these mechanisms reliably. Only common splits are used for the reported error
metrics.

The history-level bootstrap is deliberately small: five resamples per
collection. It demonstrates that replicate 007’s ASTRAL RF varies materially
when complete histories are resampled. It is not sufficient for a formal
test. The full production manifest remains unlaunched, and condition 0.15
requires a cluster resource benchmark.

## Answers to the scientific questions

1. κ=1 sometimes differs from MSC because finite structured backgrounds remain
   in the MSRC process; κ=1 only sets the cross-arrangement fraction.
2. The pilot gives no systematic suppression effect on ASTRAL topology.
3. Yes, branch-length errors vary among treatments even with unchanged
   topology, but the cause is unresolved at pilot scale.
4. History resampling in replicate 007 produces large RF variation, showing
   that finite sampling can explain a substantial part of the observed spread.
5. The outcomes are qualitatively compatible with the finite-M theory’s
   background-dependent genealogy law, but the theory does not predict
   30-taxon ASTRAL accuracy.
6. Demonstrating systematic bias requires more independently selected Mai
   species-tree replicates across all three conditions, more history-level
   bootstrap support, and explicit panmictic-limit controls.
7. These results belong in supplementary validation material, not as a main
   manuscript claim of MSRC-induced ASTRAL bias.
