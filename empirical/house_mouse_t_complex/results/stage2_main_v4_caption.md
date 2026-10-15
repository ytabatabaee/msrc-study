# Manuscript main v4 figure caption

The house-mouse t-complex comparison separates three haplotype treatments,
arrangement-conditioned quartet distributions, and balanced resampling. (A)
Three four-group ASTRAL cladograms use the same 4,046 linked 5-kb gene trees:
standard-only tips recover `Q_SPECIES`, t-haplotype-only tips also recover
`Q_SPECIES`, and the mixed treatment (standard + t) recovers `Q_T_ALT`. The
estimated internal branch lengths are CU = 0.056283, 0.202004, and 0.049097,
respectively, with localPP = 0.999985, 1.000000, and 1.000000. The biological
quartet frequencies `(q_species, q_t_alt, q_other)` are
(0.369910, 0.343656, 0.286434), (0.455385, 0.267675, 0.276940), and
(0.323323, 0.365364, 0.311313), respectively. ASTRAL inference is unrooted;
*M. spretus* is used only to orient the display. Named quartet values are
mapped to biological splits and are not fixed aliases for ASTRAL q1/q2/q3.
The terminal branch layout is schematic, while CU is the estimated internal
branch length; localPP is ASTRAL's local support statistic rather than an
independently replicated genome-wide confidence measure.

(B) Exact quartet composition for all eight domesticus/musculus/castaneus
arrangement patterns, with S denoting a standard haplotype and T denoting a
reconstructed t-haplotype sequence. Homogeneous SSS and TTT favor
`Q_SPECIES`; STS and TTS favor `Q_T_ALT`; TST favors `Q_OTHER`.

(C) Pattern weights and weighted contributions to
`D = q_species - q_t_alt`. The eight contributions reconstruct the frozen
all-haplotype quartet distribution within floating-point tolerance. Mixed-state
quartets comprise 64.6% of induced quartets, with STS and TTS providing the
strongest negative contributions.

(D) The plotted values are deterministic balanced resampling results from
1,000 replicates, not full-data ASTRAL trees. Each replicate contains 12
selected tips for B0 standard-matched and B1 t-matched, and 21 selected tips
for B2 mixed-balanced; these sample sizes differ by design. The linked 5-kb
windows are not independent biological loci. The full-data t-only result and
the resampling result agree on `Q_SPECIES`, but the historical crosses and
local genealogies come from different experimental and sampling contexts, so
the comparison does not estimate recombination suppression or establish the
MSRC mechanism.
