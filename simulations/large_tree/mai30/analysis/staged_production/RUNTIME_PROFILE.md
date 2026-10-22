# Mai30 runtime profile

The local condition-0.15 probe (`replicate 001`, one history, 20 loci) was
stopped after 562.8 seconds without a validated history.  It is recorded as
incomplete and is not counted as a simulation result.

The interrupted log points into the backward structured-genealogy loop while
querying `FrequencyHistory.frequency_at`, after the forward history had been
constructed.  This makes the main bottleneck the number of backward
generation/event steps in the long, finite-population demographic history, not
Newick I/O or ASTRAL.  The condition has larger recovered population sizes
and can have much longer ancestral intervals than the completed pilot
conditions.  The generation-by-generation Wright--Fisher process therefore
remains expensive after the `c=0.001` scaling.

`scripts/profile_runtime.py` provides an opt-in benchmark that invokes the
authoritative resumable collection runner and records wall time, exit status,
and output tails.  It does not alter the simulator, substitute a diffusion,
or write completion markers itself.  Bridges2 should run this benchmark for
one representative history from each condition before the staged array.

The completed pilot collections used approximately 96--2,130 seconds per
500-tree collection and 71--174 MB peak RSS.  Those values cannot be safely
extrapolated to condition 0.15 from the incomplete local probe.  Resource
requests and walltime must therefore be confirmed by the cluster benchmark.
