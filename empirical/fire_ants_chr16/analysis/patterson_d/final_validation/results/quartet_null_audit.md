# Rooted quartet and null audit

The original Q1/Q2 designs remain frozen. Membership comes from the explicit
Supplementary Data 1 arrangement field cross-checked against the source
individual-to-group table, exact VCF/tree-tip IDs, and sample metadata. The
selected individuals are male haploid candidates; worker pools and
random-allele-treated samples are excluded. Individual details are in
`sample_provenance_validation.tsv`.

Q1 is invicta/macdonaghi SB (P1), invicta/macdonaghi Sb (P2), richteri Sb (P3),
and *S. saevissima* (O), with 94/40/32/6 individuals. Q2 is richteri SB, richteri
Sb, invicta/macdonaghi Sb, and *S. saevissima*, with 28/32/40/6. The independent
chr1–15 ASTRAL background tree supports the within-population P1/P2 group and
the rooted ingroup/outgroup relationship; the original frozen audit reports
PASS for both checks. Thus the designs are correctly rooted for the specified
sampled groups and were not selected by observed D.

This rooted design is not equivalent to a demographic null of exchangeable
P1/P2 lineages. Arrangement, geographic origin, species ancestry (notably the
invicta/macdonaghi grouping), and selection on alternative backgrounds can
structure allele frequencies. These features complicate interpreting expected
D=0 as a universal population-history null. Geographic and outgroup subset
analyses exist and preserve the primary direction, but they do not remove
ancestry or arrangement-associated confounding.

The VCF estimator is a separate question: DP, call fraction, outgroup
unanimity, ancestral polarization, and pooled frequency products are validated
by the original tests and the additive parser tests. Neither correct rooting nor
correct arithmetic establishes introgression, MSRC, or the causal history of
recombination suppression.
