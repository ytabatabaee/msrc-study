# MSRC control-limit audit

The staged comparison uses the three MSRC switching fractions
`kappa={1.0,0.5,0.1}`.  The pilot uses baseline intensity 0.02, so the
configured effective background-switching intensity is `rho=0.02*kappa`.
Here `kappa=1.0` means that this intensity is unsuppressed relative to the
other MSRC treatments.  It is **not** an ordinary MSC control: backgrounds are still
sampled at the terminal boundary, coalescence remains background- and
frequency-dependent during the structured interval, and the genealogy uses
the finite rearrangement-origin boundary before ordinary panmictic ancestry.

The ordinary MSC collections in the pilot are retained as a separate external
baseline.  They should not be described as the `kappa=1` condition.

An exact ordinary-MSC limit is obtained by using the simulator's
`simulate_msc_genealogy` path and omitting the structured background process;
there is no claim that the existing `kappa=1` parameterization reaches that
limit.  A high-switching numerical limit experiment is prepared separately
from the staged manifest and is not used to relabel the pilot.

The current contrast is therefore causal only within the hierarchical MSRC
model: changing the effective switching fraction while holding the Mai tree,
population mapping, terminal frequency, structured interval, and sampling
design fixed.
