# Methods

The experiment uses the four-taxon species tree `((1,2),(3,4))` and imports
`simulate_population_history`, `simulate_genealogy`, and
`simulate_msc_genealogy` directly from the authoritative `msrc-sim` checkout.
The MSRC population process is the simulator's finite-M Moran process with
`M=2Ne`; `rho` is passed as the effective background-switching intensity. The
root rearrangement begins at the older end of the root branch, so the
structured interval is the root extension and the corrected simulator enters
the ordinary panmictic ancestor at the explicit rearrangement origin.

Each forward frequency history is an independent Monte Carlo unit. The
default pilot uses one genealogy per history; `--genes-per-history` can be
increased, but nested genealogies are never counted as independent histories.
The ordinary MSC comparator uses one independent genealogy per history slot.

The topology is the simulator's `topology_index`. True pairwise MRCA times,
mean pairwise time, `R_T`, internal quartet branch length, and total tree
length are calculated from the realized rooted Newick genealogy. No sequence
simulation or substitution branch-length inference is used.

The exact Moran quartet theory is evaluated through
`full_species_tree_quartet_probabilities`, including its `P(Q,Z)` recursion.
Theory comparison is restricted to topology and joint-law quantities because
the current exact operator does not derive a genealogical-time distribution.

Reproducible commands:

```text
/Users/ytabatabaee/opt/anaconda3/bin/python simulations/topology_time_separation/scripts/signal_separation.py --output simulations/topology_time_separation/results/discovery --histories 100
/Users/ytabatabaee/opt/anaconda3/bin/python simulations/topology_time_separation/scripts/signal_separation.py --validation --output simulations/topology_time_separation/results/validation --histories 1000 --genes-per-history 1 --seed 20261031
```

The configured project venv lacks SciPy, which is required by the authoritative
simulator and sparse theory operator. The recorded runs therefore used the
available Anaconda interpreter at `/Users/ytabatabaee/opt/anaconda3/bin/python`;
the simulator commit is recorded in each stage manifest.
