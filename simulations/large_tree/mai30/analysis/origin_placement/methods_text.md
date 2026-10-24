# Methods

The authoritative `msrc-sim` Wright--Fisher process generated forward
frequency histories. The optimized grouped backward simulator was used only
for the small exploratory genealogy pilot. Switching intensities were
interpreted as effective background-switching intensities; rho=0.1 and 1.0
were excluded from that pilot because their event counts were computationally
stressful. Selected quartets were prespecified by structural exposure bins,
not by observed topology outcomes. Independent forward histories, rather than
individual genealogies, are the uncertainty unit.
## Validated A/B quartet pilot

The validated extension uses the optimized `msrc-sim` structured coalescent at
commit `bfba58d627f6f467af649cdb0b37f969349e0aab`. It reuses each independent
forward Wright–Fisher frequency history across six effective background-
switching intensities and records one genealogy per history/rho/quartet.
Scenario A is the existing internal origin and scenario B is the deterministic
deeper internal origin. Quartet panels were selected from tree geometry and
forward design metadata, not from observed topologies. History-level intervals
and the multiple-comparison reversal gate are used; the pilot does not claim an
ASTRAL result.
