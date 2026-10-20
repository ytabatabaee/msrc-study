# Stage 13b matched-taxon Phlag report

The upstream Phlag 0.2.0 executable at commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6` was run without algorithm changes on deterministic one-representative-per-taxon inputs. The selected representative was held fixed across every genomic window. The basal focal edge was excluded after Phlag's own basal-bifurcation preprocessing.

The Anopheles and fire-ant pilots completed with exact taxon-set matching but returned nonfinite posterior diagnostics and are classified as numerically degenerate. They are not scored as negative detections. The house-mouse pilot was finite and informative; one focal edge and one replicate produced default state-call IoU 0.1326 and F1 0.2342. This is exploratory because it is a single representative replicate and edge.

The complete run, sampling, diagnostic, window, and accuracy tables are in this directory. The previous Atlantic-cod Phlag benchmark remains unchanged.

Matched inputs were prepared for 20 Anopheles replicates, 20 fire-ant replicates, and 7 house-mouse replicates before the remaining preparation became computationally expensive. Only replicate 000 was executed for each dataset; no unexecuted replicate is represented as an accuracy result.
