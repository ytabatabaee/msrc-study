# Multi-resolution report

The native 5-kb track is retained for local detail. The 25-kb track preserves
more short-range variation but has the largest median step size in D(w)
(0.1905) and is visibly noisy. The 100-kb
track is smoother (0.1200) but has fewer bins
and lower mean observed-window coverage (0.816).
The 250-kb track is the smoothest but attenuates the local spatial patterns.
The 50-kb track is therefore the initial main-paper candidate: its D(w) step
size is 0.1420, with mean observed-window
coverage 0.866. This choice is a
figure-readability judgment and does not use inversion coordinates or rerun
boundary detection.

All resolutions use pooled ALL_TIPS quartet counts, means of valid distance
observations, common bin boundaries, and explicit observed coverage. The
original native 5-kb and Stage 5 250-kb outputs are preserved separately.
