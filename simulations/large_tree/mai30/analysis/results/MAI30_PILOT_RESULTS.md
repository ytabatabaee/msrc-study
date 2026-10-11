# Mai30 pilot results

Status: **not executed**. The proposed MSRC pilot is prepared but remains
computationally gated.

The authenticated D1 archive is locally checksum-verified. Four 30-leaf
replicates were selected from measured original-gene-tree discordance, and
ordinary MSC controls generated 500 complete trees for each selected replicate
at four root settings. Those controls passed the output contract and showed
small distributional differences from the original SimPhy collections.

The current forward-generation Wright–Fisher implementation is too expensive
at the original Mai population sizes for the planned 500 genes per history.
Runtime measurements and the explicitly scaled, unexecuted pilot configs are
documented in:

- `MAI30_VALIDATION_RESULTS.md`
- `RUNTIME_BENCHMARK.md`
- `datasets/configs/mai30_msrc_pilot/`

The next step requires an explicit decision to use the labeled `c=0.001`
coalescent-timescale-preserving approximation or to add a validated simulator
acceleration. No full MSRC grid has been launched.
