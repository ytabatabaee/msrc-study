# Staged resource estimates

The staged design contains 30 selected species-tree replicates, three MSRC
treatments, 25 independent histories, and 20 loci per history: 90 collections
and 45,000 gene trees.  The manifest is prepared but no staged collection or
ASTRAL array has been submitted from this workspace.

The existing pilot gives a lower-bound planning range of roughly 96--2,130
seconds for a completed 500-tree collection and 71--174 MB peak RSS.  The
condition-0.15 local probe did not complete after 562.8 seconds, so its final
cost is unknown.  A Bridges2 one-history benchmark for each condition is a
decision gate before production.  Until it is run, any total walltime estimate
is provisional and should not be used to claim feasibility.

The requested SLURM scripts use one collection per array task, six concurrent
collection tasks, four concurrent ASTRAL tasks, eight GB per simulation task,
and four GB per ASTRAL task.  These are conservative starting requests based
on the pilot, not a substitute for the cluster benchmark.
