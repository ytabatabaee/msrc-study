# Stage 4E CU/SU validation

The verified ASTER source interface uses `--length CULength` for coalescent-unit branch lengths and `--length SULength` for CASTLES substitution-per-site branch lengths. SU estimation uses `--genelength` (default 1000) and `--sulengthtype` (1: species-tree branch length; 2: mean gene-tree branch length). The local smoke test is generated from a four-taxon input with nonzero gene-tree branches and validates both output modes.

No 363-species ASTER output is present in the repository at collection time; all empirical treatment and random-control comparisons therefore remain pending. Missing CU/SU estimates are retained as missing and never converted to zero.
