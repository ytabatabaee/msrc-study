# Neoaves Patterson's D / genealogy integration

This is an independent integration of the existing nucleotide ABBA--BABA
counts with the frozen chromosome-4 local-genealogy support and distance
tracks. It does not reprocess the 3.6-GB alignment archive and does not alter
the original Patterson's D or Stage 6--10 outputs.

The exact Stage-4D `locus_id` is the join key. Integer `Gene` fields in the
QQS and distance tables are first mapped through the Stage-4D manifest; no
positional join is used. N61 uses the frozen quartet
`Tauraco_erythrolophus, Columba_livia, Ciconia_maguari, Gallus_gallus`.
For this quartet, ABBA supports QQS q3 and BABA supports QQS q2, so the
nucleotide contrast is comparable to `q3-q2`, not to the existing `q2-q1`.

The PNAS intervals are retained using the repository's documented published
GalGal4 predicate on the GalGal6 locus axis. The report explicitly records
that this is an assembly-coordinate caveat rather than an exact lift.

Run from the repository root:

```bash
python3 empirical/neoaves_chr4/analysis/patterson_d_genealogy_integration/scripts/integrate.py
python3 -m unittest discover -s empirical/neoaves_chr4/analysis/patterson_d_genealogy_integration/tests -p 'test_*.py'
```

The script writes matched tracks at 100 kb, 250 kb, 500 kb, and 1 Mb,
quartet-specific branch summaries from the local named gene trees, block and
circular-shift correspondence summaries, and SVG/PDF/PNG figures. The
downloaded alignment archive is not required for this integration.
