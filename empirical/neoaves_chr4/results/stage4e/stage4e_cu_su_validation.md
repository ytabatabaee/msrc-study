# Stage 4E CU/SU validation

The verified ASTER source interface uses `--length CULength` for coalescent-unit branch lengths and `--length SULength` for CASTLES substitution-per-site branch lengths. ASTRAL-IV v1.25.4.8 reports both internal annotations on the single SULength-annotated output, so production needs one inference per treatment. SU estimation uses `--genelength` (default 1000) and `--sulengthtype` (1: species-tree branch length; 2: mean gene-tree branch length). The local smoke test is generated from a four-taxon input with nonzero gene-tree branches and validates the legacy modes; extraction tests validate the single-run schema. Terminal CU remains unestimated.

The committed T0 output has now been independently validated in
`t0_full363_validation.md`. It contains 363 taxa, 360 unique unrooted internal
splits, zero CU/SU RF distance, finite nonnegative internal CU/SU estimates,
and 363 unique terminal SU estimates. The extracted single-run artifacts
match the original CU and SU runs exactly to the recorded precision, with
terminal CU left missing. The source stderr logs report return code 0 and
optimization score `30019936712278` for both runs. The committed run metadata
does not include process end timestamps, so runtime is not reconstructed.

T_PNAS and the remaining treatment and random-control comparisons remain
pending. Missing CU/SU estimates are retained as missing and never converted
to zero.
