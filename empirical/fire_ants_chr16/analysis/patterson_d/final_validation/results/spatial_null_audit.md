# Spatial uncertainty and null audit

The existing block uncertainty table was retained and audited rather than
rerunning its bootstrap. It pools ABBA/BABA within physical bins and resamples
physical blocks; the 250-kb, 500-kb and 1-Mb results for the primary
supergene-versus-chr16A+B contrast remain positive. These are linked-genomic
block intervals, not independent-SNP standard errors. The 96 author windows
are irregularly spaced, and windows without informative data cannot add
independent information; sparse sampled windows represent a discontinuous
BUSCO track rather than contiguous sequence.

The original circular-shift p=0.0208 rotates a 52-window mask over a 96-window
author sequence. Rotations preserve window order/count but not physical
coordinates, gene density, gaps, or the observed window sampling process.
Exchangeability is therefore not established. Retain this value only as a
descriptive rotation statistic, not a conventional spatial significance test.

Chr16B has two windows and cannot support a separate regional uncertainty claim.
Chr1 is an independently sampled chromosome control, useful for background
contrast but not a matched local negative control for chromosome-16 structure,
gene density, or mapping context. Physical-block uncertainty is the principal
defensible uncertainty summary where there are enough observed blocks.
