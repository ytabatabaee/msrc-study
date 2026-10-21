# Stage 8c methods

Stage 8c runs the original TCMM implementation from
`https://github.com/shayesteh99/TCMM`, commit
`e6064b46ed0359ab7358d5110521fce1b1a825f7`, using its
`multiple_tree_matching.py` entry point. The external source is not copied into
this repository. Raw trees, branch files, rate files, optimizer logs, and
command lines are retained under `raw/`.

The query is the rooted 12-population baseline topology. The external unary
root stem in the source induced tree is set to zero in the TCMM query copy,
matching the Stage 8b induced baseline. This is a fixed-topology normalization,
not a rerooting operation. TCMM is run with lambda 0 and 0.001. Lambda 0 is the
primary no-query-branch-length setting; 0.001 is the published default
regularization setting.

Reference trees are the 250 outside and 426 all-window published 250-kb
MCC-like trees. Equal weighting is used because the available reference rows
are one tree per genomic window. No per-window posterior samples were
available, and no SNAPP or gene-tree inference was rerun.

TCMM outputs are evaluated in two ways. First, raw branch lengths and
root-to-tip distances are audited without modification. Second, descendant-set
matched node metric heights are summarized for comparison with Stage 8b. A raw
TCMM metric tree is not called a calibrated time tree unless it is ultrametric
and its units are supported by the input. Pairwise fit diagnostics distinguish
the TCMM tree-metric output from Stage 7 pairwise MRCA-time residuals.

Recorded execution command:

```text
MPLCONFIGDIR=/tmp/msrc-paper-matplotlib /Users/ytabatabaee/opt/anaconda3/bin/python empirical/atlantic_cod/scripts/08c_tcmm_baseline_projection.py --tcmm-repo /tmp/tcmm-stage8c --tcmm-python /Users/ytabatabaee/opt/anaconda3/bin/python
```

The manifest records input hashes, TCMM commit, interpreter, lambda values,
reference counts, and limitations.
