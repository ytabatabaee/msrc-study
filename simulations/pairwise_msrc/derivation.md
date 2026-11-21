# Pairwise coalescence expectations and mutation-statistic predictions under legacy MSRC

## Scope and simulator contract

This analysis uses the single-locus structured-coalescent implementation in
`ytabatabaee/msrc-sim` at commit
`1658db065a3b42a2fc1d0bebda874eecc4da6d3b`. `structured_coalescent.py`
already simulates event times for pairs and larger samples but does not expose
conditional expectations for AA, AB, and BB pairs. This document and
`scripts/pairwise_msrc.py` supply an exact expectation operator for two
lineages, while retaining the generator's existing rates and boundary rules.

The legacy process assigns each lineage one background label and uses an
effective, frequency-dependent lineage-switch rate. It is a marginal
single-locus model. The linked-genealogy SMC extension in
`docs/msrc_smc_mathematical_model.md` remains unfinished: marginal switching
alone does not specify correlated linked genealogies or segment/phase state.
No SMC/ARG claims are made here.

Let B have forward frequency \(p(t)\), and A frequency \(1-p(t)\). In the
backward process the per-lineage transitions are
\[
A\to B: a(t)=\rho(t)p(t),\qquad B\to A: b(t)=\rho(t)[1-p(t)].
\]
The structured coalescence hazards are exactly those in the validated
simulator:
\[
c_A(t)=\frac{1}{2N_e(t)\max(1-p(t),1/[2N_e(t)])},\qquad
c_B(t)=\frac{1}{2N_e(t)\max(p(t),1/[2N_e(t)])}.
\]
The floor is the simulator's finite-copy regularization. For an explicitly
conditioned pair in an absent background it allows coalescence at rate 1 in
generation-equivalent units; this is a numerical/model boundary convention,
not a literal biological population of zero copies.

## Three transient states and exact propagator

For pair states ordered as AA, AB, BB, the row generator among transient
states is
\[
Q(t)=\begin{pmatrix}
-(2a+c_A)&2a&0\\
b&-(a+b)&a\\
0&2b&-(2b+c_B)
\end{pmatrix}.
\]
The missing row mass is absorption by coalescence. For an initial state \(i\),
write the transient probability column as \(u_i(t)\), with \(u_i(0)=e_i\).
Then \(u_i'(t)=Q(t)^\mathsf{T}u_i(t)\), and survival is
\(S_i(t)=\mathbf 1^\mathsf{T}u_i(t)\). The expected time to the structured
origin boundary \(\tau\) is the area under survival:
\[
E[T_i; T_i\wedge\tau]=\int_0^\tau S_i(t)\,dt.
\]

For a constant segment of length \(d\), the augmented exponential
\[
\exp\left(d\begin{pmatrix}Q^\mathsf{T}&0\\\mathbf1^\mathsf{T}&0\end{pmatrix}\right)
\]
propagates the three transient probabilities and accumulates their integrated
survival exactly. Products of these propagators handle piecewise-constant
frequency, \(N_e\), and \(\rho\) histories without time-step approximation.

At the finite structured-origin boundary labels disappear. If the ordinary
ancestral MSC after that boundary has expected residual waiting time \(W_0\),
then
\[
E[T_i]=\int_0^\tau S_i(t)\,dt + S_i(\tau)W_0.
\]
The boundary age is not added again: it is already represented in the
structured survival integral. For ordinary ancestral segments \((d_j,N_j)\),
with \(N_j=\infty\) representing no coalescence, the tail having size
\(N_{anc}\),
\[
W_0=\sum_j S_j\frac{1-e^{-d_j/(2N_j)}}{1/(2N_j)}
       + S_{end}(2N_{anc}),
\]
where a no-coalescence segment contributes \(S_jd_j\). These formulas yield
\(E[T_{AA}],E[T_{AB}],E[T_{BB}]\) by starting at each basis state. Variable
frequency histories are integrated conditionally on the full forward history;
uncertainty is then summarized over independent frequency histories, rather
than treating multiple genealogies that share one frequency history as
independent replicates.

## Mutation statistics and the ordinary null

Under neutral infinite sites, a common per-generation mutation rate \(\mu\),
random pair sampling within each background, and no additional within-group
structure beyond the modeled genealogy,
\[
E[\pi_A]=2\mu E[T_{AA}],\quad E[\pi_B]=2\mu E[T_{BB}],\quad
E[d_{AB}]=2\mu E[T_{AB}].
\]
Thus
\[
E[\Delta_{AB}]=2\mu\left(E[T_{AB}]-\tfrac12(E[T_{AA}]+E[T_{BB}])\right),
\]
and the ratio of expectations is
\[
R_{AB}^{(mom)}=\frac{2E[d_{AB}]}{E[\pi_A]+E[\pi_B]}
=\frac{2E[T_{AB}]}{E[T_{AA}]+E[T_{BB}]}.
\]
This is the large-window moment ratio, not generally
\(E[2d_{AB}/(\pi_A+\pi_B)]\). A matched ordinary panmictic coalescent has
equal within- and between-group pair means, hence \(E[\Delta]=0\) and
\(R^{(mom)}=1\). A geographic or otherwise structured ordinary-coalescent
null need not have those values and must be matched to the study's sampling
design before interpreting positive excess as an MSRC-specific effect.

With identical callable sites and consistent weighting, Hudson's ratio form
obeys
\[
F_{ST}=\frac{d_{AB}-(\pi_A+\pi_B)/2}{d_{AB}}=1-\frac{1}{R_{AB}}.
\]
Consequently \(F_{ST}\), \(R\), and \(\Delta\) reuse the same diversity
information; they are not three independent confirmations. For finite windows,
ratios of random quantities also differ from ratios of expectations.
Mutation-rate equality, neutral mutation, correct sample assignments,
haploid/diploid sampling treatment, and consistent callability are assumptions.
Background-conditional samples must be genuine background chromosomes; a
heterozygous carrier cannot be treated as a pure background without resolved
phase.

## Special cases and numerical validation

When \(\rho=0\), the states do not switch. For AA and BB, the expectation is
\((1-e^{-c_i\tau})/c_i+e^{-c_i\tau}W_0\); AB cannot coalesce before the
labels disappear, so it has \(\tau+W_0\). At balanced frequency, exchanging
background labels leaves the expectations unchanged. Strong switching
typically removes dependence on the initial pair state over a sufficiently
long structured phase. At \(\tau=0\), all three means equal the ordinary MSC
mean.

Eight `unittest` checks cover these limits, symmetry, time-scale invariance,
the finite-copy floor, the algebraic FST/R identity, and direct event-time
simulation against `simulate_genealogy`. The direct fixed-history check draws
2,500 independent genealogies per pair class; every exact mean is within 4.5
Monte Carlo standard errors. The frequency-path validation runs Moran and
Wright–Fisher histories under neutral and selected paths (100 independent
histories/cell, initial B frequencies 0.25/0.5/0.75, \(2N_e=20\),
\(\tau=2N_e\), \(2N_e\rho=1\)); each history has one direct genealogy per
initial pair state. The analytic-history versus direct-genealogy difference
interval includes zero in 35/36 cells. This is a numerical validation at small
population size, not empirical calibration.

`theoretical_parameter_grid.tsv` varies \(p\), \(2N_{e,in}\rho\),
\(\tau/(2N_{e,in})\), and \(N_{e,in}/N_{e,out}\) over broad fixed-history
ranges. In this grid, R decreases monotonically as switching increases for
every tested combination, so reducing switching (the direction expected
under stronger effective suppression) increases R. This is a result for this
specific fixed-frequency legacy generator and grid, not a universal
monotonicity theorem. Suppression also changes R with duration, frequency, and
the relative effective sizes of the structured and ancestral phases; it can
yield large divergence while retaining high within-background diversity, but
the fit is not unique. Parameter identifiability must be assessed jointly.

## Anopheles 2La comparison and identifiability

The empirical input is the frozen 50-kb AgamP4 2L pairwise table in
`empirical/anopheles_2la/analysis/joint_genomic_tracks/results/`. Its 988
windows include 432 inside 2La, 554 outside, and two boundary-overlap windows
excluded from regional contrasts. For each primary comparison, the source
records the same callable-site count for dXY and Hudson FST; the stored
cross-population \(\pi_X,\pi_Y\) were also evaluated on that pair's eligible
sites. Standalone π scans have different masks and were not substituted.
Statistics are per callable nucleotide, invariant positions included. The
coordinate system is 1-based inclusive AgamP4; the interval is
2L:20,524,058–42,165,532. Physical 500-kb blocks are the spatial uncertainty
units. This interval estimates local spatial variation only, not uncertainty
over populations, geography, or sample ascertainment.

The arrangement groups were selected homozygotes, so their sample fractions
are not estimates of natural \(p(t)\). A-distance and local tree metrics are
NJ patristic distances in substitution units, not calibrated coalescence
times; quartet support is a topology probability. Their overlap with diversity
patterns is examined descriptively, with block resampling, and not equated to
pair-time estimates. Outside 2La is a biological comparison, not an assumed
panmictic null: species, geography, introgression, and arrangement-associated
sampling can all contribute to its small nonzero divergence excess.

The plotted recombination literature gives regional rate comparisons
(Stump et al. 2007); heterokaryotype recombination is strongly reduced in the
inversion and flanks but not necessarily absent. Centimorgans per megabase are
not converted into \(\rho\). At most these data inform sensitivity about the
direction and possible scale of suppression in an effective switching rate.

From \(\pi\) and \(d_{XY}\) alone, mutation rate and time scale appear as
products, and even after fixing a mutation rate, \(N_e\), structured duration,
frequency trajectory, and switching are confounded. R removes a common
mutation-rate/time scale but constrains only a ratio of pair-time means. The
regional means provide three overlapping statistics per group and cannot
uniquely recover \(N_e,\rho,\tau,p(t)\); extra independent demographic,
mutation, geography, and genealogy data are needed. The current fixed-p
parameter grid is sensitivity analysis, not a fitted posterior. A fully
identified joint fit would need a calibrated genealogy observable and a
geographic/introgression null, neither supplied by the existing NJ branch
units.

## Applicability to the other empirical systems

Existing fire-ant, Neurospora L10, and house-mouse t-complex window analyses
are SNP-panel conditional: absent positions are not callable invariant
bases. Their absolute \(\pi\), \(d_{XY}\), and callable-base rate comparisons
remain blocked by unavailable or incomplete masks, so no absolute
theory-to-data prediction test is reported. Their conditional SNP statistics
are not substituted into the infinite-sites per-base equations. For Atlantic
cod, local population tree divergences are calibrated SNAPP population-tree
times, not gene-lineage pair coalescence times; the existing nucleotide
summaries are SNP-ascertained. No shared parameterization is imposed across
these systems.

## Source

Stump AD et al. (2007), "The reproductive isolation of the *Anopheles gambiae*
species complex revealed by the analysis of introgression loci," including
recombination estimates. [PubMed record](https://pubmed.ncbi.nlm.nih.gov/18092999/).
