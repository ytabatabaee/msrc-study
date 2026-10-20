# Phlag input compatibility

The Stage 6 values are normalized QQS summaries for three source-defined Neoaves focal designs. Phlag's `--read-qqs-path` format is indexed by labeled species-tree edges and by gene-tree position. The stored Stage 6 rows do not contain a universal Phlag edge identity, and the three designs use different quadripartitions, so they are not silently passed to Phlag as though they were one Phlag sequence.

A matched Columbea pilot is prepared from the original local trees using the fixed representatives Columba_livia, Phoenicopterus_ruber, Tauraco_erythrolophus, and Gallus_gallus. The source 63K species tree is pruned to these tips and its focal internal edge is labeled. This preserves original local-tree observations while making the taxon set exact. Any resulting method comparison is limited to this pilot sampling design.
