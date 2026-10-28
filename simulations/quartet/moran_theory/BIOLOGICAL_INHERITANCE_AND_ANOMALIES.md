# Biologically constrained inheritance and quartet anomalies

This audit separates three mechanisms: identity inheritance, symmetric
regime inheritance, and a founding-count approximation. It also compares the
two-regime quartet law with a standard two-component network MSC mixture.

## Inheritance models

Identity inheritance uses `H=I`. Symmetric inheritance uses

\[
H_\epsilon=\begin{pmatrix}1-\epsilon&\epsilon\\
\epsilon&1-\epsilon\end{pmatrix}.
\]

The founding model is specified before looking at quartet probabilities. A
parent in regime `e` has frequency `p_e`; a daughter founder count is sampled
as `X~Binomial(M,p_e)`, and the daughter is classified low when
`X/M <= (p_low+p_high)/2`. The resulting two-by-two probabilities are used as
the daughter regime kernel. This is a coarse two-regime approximation to a
founding process, not a derivation of the exact finite-M Moran inheritance
law.

## Targeted search

The search evaluates 864 deterministic cases: two rooted shapes, three
frequency pairs, three switching rates, three regime-transition settings, two
initial regime weights, two copy counts, and four inheritance models. The
complete table is `results/biological_inheritance/search.tsv`.

The largest anomaly margins in this grid were negative for every constrained
model:

| inheritance | balanced max `max(q2,q3)-q1` | unbalanced max |
|---|---:|---:|
| identity | -0.08810 | -0.01746 |
| symmetric, epsilon=0.02 | -0.08858 | -0.01817 |
| symmetric, epsilon=0.10 | -0.08981 | -0.01992 |
| founding classification | -0.09012 | -0.01864 |

These are numerical findings on the stated grid, not impossibility theorems.
Identity inheritance can nevertheless produce small unconditional discordant
asymmetry in the unbalanced tree; the maximum in the grid was 0.01251. Thus
`q2 != q3` and `q2/q3 > q1` are distinct claims.

The previously reported anomaly remains reproducible under the unrestricted
inheritance model. Its exact witness is

\[
q=(0.35417018,0.26868940,0.37714042),
\]

with `q3-q1=0.02297024` and `|q2-q3|=0.10845102`. It uses nonidentity,
taxon-shared regime kernels. Dense matrix-exponential recomputation agrees to
machine precision. The witness therefore establishes a model-level numerical
possibility, but it does not establish that biologically constrained
inheritance produces anomalies.

## Introgression comparison

For a standard four-taxon network MSC mixture, let the base tree favor `Q1`
and let an introgressed parental edge favor `Q2` or `Q3`. With inheritance
weight `gamma` and internal branch lengths `t_base,t_intro`, the quartet law is

\[
q=(1-\gamma) m_{Q1}(t_{base})+\gamma m_{Qj}(t_{intro}),
\]

where the favored topology has probability `1-2 exp(-t)/3` and each other
topology has probability `exp(-t)/3`. This is the standard quartet-mixture
calculation used for network MSC quartet concordance factors; see
[Solís-Lemus and Ané (2016)](https://doi.org/10.1371/journal.pgen.1005896)
and the discussion of quartet CFs under hybridization in the accompanying
[PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC4780787/).

The unrestricted MSRC witness is matched exactly by an introgression mixture
favoring `Q3` with

```text
gamma   = 0.5941079083
t_base  = 0.2364818452
t_intro = 0.2015585930
max error = 0
```

This is a constructive nontrivial equality of quartet vectors, but it compares
introgression with the unrestricted-kernel witness. It does not show that
identity or founding inheritance has an anomalous quartet. A second exact fit
is recorded for an identity-inheritance vector,
`(0.4329387874,0.2836441797,0.2834170329)`, using
`gamma=0.05800342` and an introgressed `Q2` component; this demonstrates that exact quartet-vector overlap does
not require an anomalous MSRC vector.

The equality is an existence result. No positive-area overlap theorem or
full-simplex coverage theorem has been established. A rank-two parameter
Jacobian would not prove either claim.

## Moran comparison and interpretation

The original finite-M Moran comparison in
`results/two_regime_reachability/model_comparison.tsv` gives, for one matched
identity-inheritance setting,

```text
two-regime       (0.36371581, 0.31822689, 0.31805730)
finite-M Moran   (0.36468612, 0.31765694, 0.31765694)
```

The two models are not equivalent. The comparison does not demonstrate an
unconditional Moran anomaly. No anomaly was found in the targeted Moran cases
examined here.

The rigorous conclusions are limited. The 38-state model is coherent and can
produce unconditional anomalies when arbitrary inheritance kernels are
allowed. Under identity, small symmetric perturbations, and the prespecified
founding classification, this search found no anomaly. A proof that identity
inheritance always preserves species-tree dominance remains open. Conditional
anomalies and unconditional discordant asymmetry should not be conflated.
