# Stage 13 unified boundary detection

Run from the repository root:

```sh
MPLBACKEND=Agg python empirical/cross_system_boundary_detection/scripts/13_unified_spatial_boundary_detection.py
PYTHONPATH=empirical/shared python -m pytest -q empirical/shared/test_spatial_boundary_detector.py empirical/cross_system_boundary_detection/tests/test_stage13.py
```

The detector preserves observed coordinates and gaps, and keeps SNAPP
population-tree units separate from substitution-distance units. Circular
shift and change-point nulls test different hypotheses.
