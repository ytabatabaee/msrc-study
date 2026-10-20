# Fire-ant chromosome-16 patristic distance track (Stage 9)

This analysis uses the 213 published RAxML local trees and the verified Stage-6C individual-to-group mapping. Branch-length distances are patristic distances in substitution units. They are not coalescence times and no molecular-clock assumption is made.

The primary contrast is group-balanced: each of the four opposite-haplotype group pairs is averaged within pair and then equally across pairs; each of the two same-haplotype cross-species pairs is treated analogously. The supergene/outside comparison uses 52 and 44 chromosome-16 windows with complete focal groups.

## Results

- `A_distance`: supergene mean=0.0017179535494069613, outside mean=-0.003446625419768419, difference=0.00516457896917538.
- `A_relative`: supergene mean=0.24172471843367566, outside mean=-0.3029444963355287, difference=0.5446692147692044.
- `A_distance_individual_weighted`: supergene mean=0.0012692479905132047, outside mean=-0.00345070180343527, difference=0.004719949793948475.
- `A_relative_individual_weighted`: supergene mean=0.17305036168519752, outside mean=-0.33626011939360995, difference=0.5093104810788075.
- `A_distance_outgroup_normalized`: supergene mean=0.18646732232632265, outside mean=-0.1705563458511423, difference=0.35702366817746495.

The individual-pair-weighted values are a sensitivity analysis and can be influenced by unequal group sample sizes. The outgroup-normalized sensitivity divides the group-balanced distance contrast by the mean distance from focal individuals to nonfocal individuals observed in the same local tree.

## Spatial null

- `A_distance` / `exact_physical_order_circular_shift`: observed supergene-minus-outside contrast=0.00516457896918; p=0.0208333333333 over 96 alignments or coordinate states.
- `A_distance` / `unique_membership_state_coordinate_sensitivity`: observed supergene-minus-outside contrast=0.00516457896918; p=0.0309278350515 over 97 alignments or coordinate states.
- `A_distance` / `continuous_coordinate_length_weighted_sensitivity`: observed supergene-minus-outside contrast=0.00516457896918; p=0.0241243813869 over 97 alignments or coordinate states.
- `A_relative` / `exact_physical_order_circular_shift`: observed supergene-minus-outside contrast=0.544669214769; p=0.0104166666667 over 96 alignments or coordinate states.
- `A_relative` / `unique_membership_state_coordinate_sensitivity`: observed supergene-minus-outside contrast=0.544669214769; p=0.020618556701 over 97 alignments or coordinate states.
- `A_relative` / `continuous_coordinate_length_weighted_sensitivity`: observed supergene-minus-outside contrast=0.544669214769; p=0.0603060648818 over 97 alignments or coordinate states.
- `A_distance_individual_weighted` / `exact_physical_order_circular_shift`: observed supergene-minus-outside contrast=0.00471994979395; p=0.0208333333333 over 96 alignments or coordinate states.
- `A_distance_individual_weighted` / `unique_membership_state_coordinate_sensitivity`: observed supergene-minus-outside contrast=0.00471994979395; p=0.0309278350515 over 97 alignments or coordinate states.
- `A_distance_individual_weighted` / `continuous_coordinate_length_weighted_sensitivity`: observed supergene-minus-outside contrast=0.00471994979395; p=0.0241243813869 over 97 alignments or coordinate states.
- `A_relative_individual_weighted` / `exact_physical_order_circular_shift`: observed supergene-minus-outside contrast=0.509310481079; p=0.0104166666667 over 96 alignments or coordinate states.
- `A_relative_individual_weighted` / `unique_membership_state_coordinate_sensitivity`: observed supergene-minus-outside contrast=0.509310481079; p=0.020618556701 over 97 alignments or coordinate states.
- `A_relative_individual_weighted` / `continuous_coordinate_length_weighted_sensitivity`: observed supergene-minus-outside contrast=0.509310481079; p=0.0603060648818 over 97 alignments or coordinate states.
- `A_distance_outgroup_normalized` / `exact_physical_order_circular_shift`: observed supergene-minus-outside contrast=0.357023668177; p=0.0104166666667 over 96 alignments or coordinate states.
- `A_distance_outgroup_normalized` / `unique_membership_state_coordinate_sensitivity`: observed supergene-minus-outside contrast=0.357023668177; p=0.020618556701 over 97 alignments or coordinate states.
- `A_distance_outgroup_normalized` / `continuous_coordinate_length_weighted_sensitivity`: observed supergene-minus-outside contrast=0.357023668177; p=0.0603060648818 over 97 alignments or coordinate states.

## Limitations

The source trees are substitution-unit RAxML trees with irregular four-BUSCO window spacing. The track is descriptive and spatially tested against the frozen author-labeled supergene window mask; the interval is an observed BUSCO-window span rather than a nucleotide-resolved breakpoint annotation.
