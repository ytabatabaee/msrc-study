# Stage 13b methods

Stage 12b Model B was reproduced from the frozen implementation and outputs before cross-system comparison. Its shared-background objective uses one outside mean, robust median/MAD/IQR scaling with an SD floor, the original penalty, censored endpoints, and deterministic tie-breaking. Stage 13 independent-segment output is retained as a separate formulation.

Primary genomic IoU is calculated from merged interval sets using half-open arithmetic internally; disconnected predictions are not replaced by their hull. Censored observation endpoints remain labeled as censored. The current spatial-null table records the prespecified ordered-block design and explicitly withholds p-values where full detector recomputation was not completed.

Phlag uses upstream version 0.2.0 at commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6`, with one fixed representative per taxon across windows. SB/Sb labels were renamed to neutral unique lineage labels because Phlag/DendroPy treats case variants as duplicate taxa. The basal focal edge is excluded after Phlag’s own basal-bifurcation collapse.

# Stage 13b results

The frozen Atlantic cod Stage 12b IoUs were reproduced exactly: LG01 0.9676, LG02 0.9444, LG07 0.9566, and LG12 0.9147. The mean is 0.9458. Stage 13 differs because it permits independent flank means and uses a different scaling/penalty formulation; these objectives are not interchangeable.

Matched-taxon Phlag pilots were completed on Anopheles, fire ants, and house mouse. Anopheles and fire-ant executions returned finite exit status but nonfinite posterior diagnostics, so they are classified as numerically degenerate rather than negative detections. The house-mouse pilot returned finite posterior probabilities and a nonconstant state path.

The full 20-replicate matched-taxon benchmark remains a follow-up because each upstream Phlag run takes roughly one to two minutes and the two non-mouse pilots are numerically degenerate under the original branch-length inputs.
