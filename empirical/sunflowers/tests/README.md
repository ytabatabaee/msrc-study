# Sunflower Stage 1 tests

The focused suite covers fixed individual-quartet extraction, the three induced quartet resolutions and internal-branch distance calculation, unresolved quartets, pairwise-distance summaries, genomic interval overlap, conservative adjacency, and run detection that does not bridge missing intervals.

Run with:

```sh
python3 -m unittest discover -s empirical/sunflowers/tests -p 'test_stage1_spatial_analysis.py' -v
```
