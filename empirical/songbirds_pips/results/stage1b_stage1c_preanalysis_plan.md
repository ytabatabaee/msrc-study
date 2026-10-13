# Stage 1C preanalysis plan placeholder

Stage 1C is not justified by the present Stage 1B run because the focal raw-read pilot did not execute and the published shared-PIP signal was not reproduced.

The intended Stage 1C question remains: is genotype-across-species structure at `shared_CM027530.1` stronger than in matched non-PIP regions?

A final control-window rule must be frozen only after Stage 1B produces usable focal-region data, stable 231-pair SFS estimation, and a complete or nearly complete dXY matrix. Until then, no background windows should be selected, no control NJ trees should be inspected, and no Stage 1C result should be generated.

Provisional independent eligibility constraints, to be revisited after successful Stage 1B:

- same scaffold `CM027530.1`;
- 300 kb windows outside all known PIP intervals on the scaffold;
- exclude low-mappability or low-callability windows using thresholds defined from coverage/QC, before inspecting any control distance tree;
- target 10-50 controls only if computationally manageable;
- primary statistic: focal genotype-vs-species distance contrast compared with the distribution from eligible controls.
