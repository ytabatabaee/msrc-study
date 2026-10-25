# Checkpoint/resume validation

The checkpoint key is `(window_id, test_id, filter_label)`. Existing rows are
validated against the frozen window coordinates and design identities before
they are accepted. Unknown keys, duplicates, truncated rows, coordinate
mismatches, and design mismatches raise an error.

TSV and JSON checkpoint writes use a temporary sibling followed by `os.replace`,
so an interrupted write cannot replace a valid completed table. Before a
regional query, the runner computes missing keys. A complete region is skipped;
when only some keys are absent, only those rows are calculated. The final
completion and provenance counts are computed from the complete validated table,
not just rows made in the current process.

The deterministic tests cover an empty checkpoint, a partially complete
region, invalid coordinates, duplicate-key rejection, and a second restart with
all expected keys complete. They also verify that the eight pilot regions are
exactly eight independent 50-kb requests totaling 400,000 bp. No real pilot
rows were overwritten during this audit.

The remaining limitation is that the test backend is in-process rather than a
live network mock; the production skip decision is pure key arithmetic and is
independent of MalariaGEN. A full run should preserve the `.partial` file, if
present, for diagnosis rather than treating it as a completed result.
