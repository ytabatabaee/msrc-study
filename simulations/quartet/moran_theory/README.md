# Exact finite-M Moran-integrated shared-population quartet kernel

This directory implements the single shared ancestral population calculation
for the hierarchical model in `msrc-sim`: a forward finite-population Moran
frequency history and a backward structured genealogy conditional on that
history. It does not implement the full species-tree recursion or sister
branches.

The main implementation is `moran_quartet_operator.py`. It constructs the
finite Moran generator, the 19-state conditional quartet generator, and the
sparse dual operator

`A4 = Q_M kron I_19 + blockdiag(G4(0)^T, ..., G4(M)^T)`.

The transpose is essential. A realized Moran path is ordered old to young,
whereas a genealogy sampled at the young boundary is traversed young to old.
The operator is a dual/Feynman--Kac average over these two directions; it is
not being claimed as the microscopic ancestral dual of the Moran reproduction
events.

Run unit tests from the repository root with:

```text
pytest simulations/quartet/moran_theory/tests
```

The standalone validation imports the authoritative `/Users/ytabatabaee/Desktop/msrc-sim`
checkout by default. A compact smoke run is:

```text
python simulations/quartet/moran_theory/validate_against_msrc_sim.py --smoke --replicates 1000
```

Increase `--replicates` and use the default grid for publication-scale
Monte Carlo summaries. Outputs are written to `results/`.

The full rooted species-tree recursion is in `full_species_tree.py`, with its
derivation in `FULL_SPECIES_TREE_DERIVATION.md`. Its validation driver is:

```text
python simulations/quartet/moran_theory/validate_full_species_tree.py --replicates 500
```

The recursion preserves the common parent count through each speciation and
multiplies left/right sister contributions only conditionally on that count.
The derivation documents a finite-M boundary distinction: literal zero-
duration `msrc-sim` speciation still performs binomial daughter initialization,
so it is distinct from the separate collapsed shared-population kernel.
