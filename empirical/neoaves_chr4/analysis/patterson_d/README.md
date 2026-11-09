# Neoaves chromosome-4 nucleotide Patterson's D

This directory contains an independent ABBA--BABA analysis of the original
Stiller et al. intergenic nucleotide alignments. It is deliberately separate
from the existing local-tree quartet statistic `D = q2 - q1`.

The frozen PNAS outlier intervals are the GalGal4 window-start intervals used
by the published `delta_quartet.R` predicate:

* 25,030,000--32,670,000
* 33,510,000--34,470,000
* 44,130,000--56,810,000

The local Stage-4D manifest supplies the exact Stiller locus IDs and GalGal6
coordinates. The original alignment archive is kept outside Git under
`empirical/neoaves_chr4/external/stiller2024/`.

## Reproduction

```bash
# The raw archive is downloaded only if it is absent and is checksum-verified.
bash empirical/neoaves_chr4/analysis/patterson_d/scripts/download_alignments.sh
<configured-python> empirical/neoaves_chr4/analysis/patterson_d/scripts/audit_inputs.py
<configured-python> empirical/neoaves_chr4/analysis/patterson_d/scripts/patterson_d.py all
<configured-python> -m unittest discover -s empirical/neoaves_chr4/analysis/patterson_d/tests -p 'test_*.py'
```

`patterson_d.py` reads the tar archive directly, so extraction is not needed.
It writes locus-level counts, pooled region summaries, block-jackknife
uncertainty, and SVG/PNG/PDF figures. All D values pool ABBA and BABA site
counts before taking their ratio.

The taxon designs are frozen in `config/frozen_quartets.tsv` before any D
values are calculated. The primary outgroup is *Gallus gallus*; the alternative
outgroup is *Meleagris gallopavo*, another galliform taxon present in the
Stiller taxon map. The report records coverage and whether each alternative is
usable for each design.

