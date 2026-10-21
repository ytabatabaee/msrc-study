# Diagnosis of current boundary detectors

The frozen cross-system detector uses
`empirical/shared/spatial_boundary_detector.py`. It robustly standardizes each
feature with a median/MAD scale, keeps observed coordinates and gaps, and uses
exact dynamic programming for a fixed number of piecewise-constant segments.
The selected interval count is penalized. The Stage 13b shared-background
comparison uses a different objective and scaling from the Stage 13 independent
segment comparator; `stage13b/results/detector_objective_audit.tsv` records the
difference.

The results show why topology-only, continuous-only, and joint fits diverge:

* The blocks have different estimands and source units. Topology is a simplex
  or contrast of fractions; continuous features are SNAPP time summaries or
  substitution distances. Robust standardization fixes a numerical scale but
  does not make their noise models or biological directions equal.
* A joint SSE objective rewards a segment that explains total squared residual
  reduction. If topology and length changes occur in different locations or
  have opposite directions, their sum can choose an interval that is good for
  neither modality. The Atlantic cod LG12 candidate audit demonstrates this:
  topology improvement and time improvement favored different candidates.
* A single-interval model cannot represent separated or heterogeneous
  intervals. The house-mouse Stage 13 results explicitly show that the
  multi-interval reference is not represented by one outer hull.
* Interval penalties and endpoint censoring affect modalities differently when
  their missing windows differ. The code requires complete rows for the chosen
  feature block, so joint fits can use a smaller support than either marginal
  fit.
* Branch-length tracks can be driven by a terminal taxon, SNP density, or
  substitution-scale artifact. The mouse Stage 6 diagnostic found a
  taxon-specific downstream discontinuity in `B(w)`. A joint detector can
  correctly optimize that signal and still be biologically misleading.

The current results support a methodological diagnosis, not a universal
ranking. Stage 13 means topology versus continuous union IoU/F1 were cod
`0.627/0.678` versus `0.784/0.846`, fire ants `0.530/0.693` versus
`0/0`, Anopheles `0.995/0.997` versus `0.997/0.999`, and mouse
`0.143/0.269` versus `0.117/0.200`; the joint detector did not improve over
the better single modality. These comparisons use frozen structural regions
only for post-fit evaluation.

The Atlantic cod Stage 12b result adds a more specific warning: the original
joint objective can select an interval with incompatible modality directions,
while a prespecified robust equalized objective behaves better on LG12. This
is evidence for objective sensitivity, not permission to tune on structural
coordinates.

Phlag remains a frozen comparator. House-mouse Stage 13d and Neoaves Stage
11C/Stage 11B should be evaluated on matched domains, with matched four-taxon
and full multi-taxon inputs kept as separate analyses.
