# Reachable unconditional quartet distributions in the two-regime model

This audit concerns the finite 38-state two-regime model, not the original
finite-M Moran model. A regime is a two-state CTMC whose states have fixed
background frequencies `p_low` and `p_high`. The regime generator is combined
with the transposed conditional genealogy generators. The model is therefore an
exact marginal genealogy for its stated regime process, with a maximum
four-lineage operator dimension of 38.

## Implementation audit

`two_regime_species_tree.py` uses the forward regime generator `R` in
`R \otimes I` and the transposed backward genealogy blocks in the diagonal
part. The initial regime distribution is applied at the older end of the root
branch. Row-stochastic inheritance kernels are applied at speciation before
descendant contributions are multiplied, so daughter regime contributions
remain correlated through their common parental regime. The early terms retain
terminal pattern mass when a topology is resolved before all four lineages have
entered the ancestral root. These identities are tested by joint-law
normalization and marginalization tests.

The effective parameter `rho` is an ancestry-background switching rate for one
locus. It is not a chromosome-wide recombination breakpoint rate.

## Exact identities

If `p_low = p_high = p`, every genealogy block is the same. The regime process
then factors from the genealogy and, because `R` is conservative, integrates to
one. The result is exactly the fixed-frequency genealogy at `p`, independent of
the regime switching rates and initial regime distribution.

If the regime rates are zero and inheritance is identity, the result is the
initial-regime-weighted mixture of the two fixed-frequency calculations. This
does not imply that arbitrary nonidentity inheritance is a mixture of fixed
regime models, because inheritance can change the regime seen by descendants.

An equality `q2=q3` follows whenever the complete parameterized model is
invariant under a taxon permutation that exchanges the two discordant quartet
resolutions, including branch lengths, population parameters, terminal
sampling, and inheritance kernels. Complementing background labels alone does
not exchange taxa and therefore does not prove this equality.

## Bounded reachability search

The reproducible search in `analyze_two_regime_reachability.py` evaluates 100
random parameter points for each rooted shape. It varies both regime
parameters and row-stochastic inheritance kernels. It is a diagnostic search,
not a proof of coverage or of an optimization bound.

The recorded witness has an unbalanced tree and gives approximately

```text
q = (0.35417018, 0.26868940, 0.37714042)
max(q2,q3)-q1 = 0.02297024
|q2-q3| = 0.10845102
```

Thus an unconditional anomaly and substantial unconditional discordant
asymmetry are both possible in this model. The witness uses nonidentity
inheritance kernels; the bounded identity-only search did not find a positive
anomaly, but that observation is not a theorem.

At the representative interior point, finite differences in `(p_low,
p_high,rho)` have Jacobian rank 2 for `(q1,q2)`, with singular values reported
in `summary.tsv`. This establishes local two-dimensional variation of the
simplex image at that point, not full-simplex coverage. No open-region or
full-simplex theorem is established here. The finite search also cannot prove
that a simplex vertex is attainable or unattainable.

The lower panels of the generated simplex figure show the conditional law for
the anomaly witness and the unconditional search points. Conditional point
areas are proportional to pattern probability; they are not confused with
unconditional points.

## Comparison with the Moran model

`model_comparison.tsv` contains one matched, identity-inheritance comparison
with the exact finite-`M` Moran solver. It is a benchmark, not an equivalence
claim. The two models have different frequency processes: the Moran model
retains all count states and their history, while the regime model has only two
fixed-frequency regimes. Similar quartet vectors at one parameter setting do
not establish distributional equivalence.

## Conclusions and limitations

The two-regime model is mathematically coherent and its unconditional formulas
are exact for that model. It can produce both unconditional anomalies and
discordant asymmetry. The strongest rigorous statements currently available
are the factorization and symmetry propositions above. Full-simplex coverage,
global bounds on the anomaly margin, and equivalence to the finite-`M` Moran
model remain unresolved. The model is suitable as a clearly labeled reduced
stochastic alternative, not as a replacement for the original Moran theory.
