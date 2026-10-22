# Methods

Four fixed Mai30 pilot species trees were analyzed. For each tree, four origin
placements were defined deterministically: the original internal origin, the
oldest eligible non-root internal origin, a root-extension origin 50 scaled
generations before root speciation, and a doubled-root-extension origin near
the older boundary. All placements preserve the simulator's origin-age
convention.

The full set of 27,405 four-taxon subsets was used for structural exposure.
Forward histories were simulated with `msrcsim.population_process` using the
scaled Wright--Fisher process, neutral selection, and initial frequency 0.5.
Terminal arrangement patterns were sampled from the resulting histories.

The attempted genealogy discovery grid used rho values from 0 through 1,
including the original 0.002, 0.01, and 0.02 values. Local backward
simulation was stopped after runtime stress in 30-taxon ancestral simulation;
the executed genealogy result files therefore contain explicit unexecuted
status rows rather than extrapolated quartet probabilities. No ASTRAL
collection was generated for an origin-placement scenario.
