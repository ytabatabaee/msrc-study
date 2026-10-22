# Control audit

The four pilot treatments use the same four fixed Mai D1 species-tree
backbones, the same 30 ingroup taxa, the same scaled population-size mapping,
the same scaled branch lengths, and 500 gene trees per collection. MSRC
collections contain 25 independent forward Wright–Fisher frequency histories
with 20 loci per history. The ordinary MSC controls contain 500 directly
simulated gene trees; they do not contain shared forward frequency histories.

ASTRAL-IV was run with the same command settings for every collection:
`-t 2 -u 1 --length CULength`. All input collections and outputs passed the
complete-genealogy contract. Branch errors in the analysis are calculated only
on splits present in both the inferred and true species trees.

The primary suppression contrast is within MSRC:

\[
\Delta(\kappa)=M_{\mathrm{MSRC},\kappa}-M_{\mathrm{MSRC},1}.
\]

The ordinary MSC is reported separately. It is not a paired baseline for this
contrast because it has no arrangement backgrounds and was generated from a
different stochastic process.

## Why κ=1 is not ordinary MSC

In these configurations `kappa` is the effective cross-arrangement fraction
in the background-switching rate, with baseline rate 0.02. At κ=1,
background switching is present but finite; the sampled arrangement labels,
finite-frequency history, same-background coalescence rules, and structured
root interval remain. Ordinary MSC has no background process and uses direct
panmictic coalescence. Thus κ=1 can differ from MSC in either direction.

The configurations use `ordinary_msc_outside: true`, so unresolved lineages
enter an ordinary ancestor after the modeled structured interval. This does
not remove the structured dynamics before that boundary. No existing control
in this pilot is an exact panmictic limit of the full MSRC model; that would
require a separate limit experiment with background switching sufficiently
fast relative to structured coalescence and matched ancestry.

The exact finite-M quartet theory supports background-dependent topology laws,
but it does not derive ASTRAL behavior for this 30-taxon large-tree setting.
The ASTRAL comparisons are therefore simulation results.
