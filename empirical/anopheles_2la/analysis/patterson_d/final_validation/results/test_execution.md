# Test execution

The installed interpreters do not have `pytest` (`python -m pytest` and the repository `.venv-malaria310/bin/python -m pytest` both report `No module named pytest`). The current default Python 3.7 cannot import the existing test module because the analysis runner imports `importlib.metadata`, introduced in Python 3.8. Using the installed Python 3.10 runtime, all 18 existing and new test functions were invoked directly, including the two existing fixture-based tests with equivalent temporary-directory and monkeypatch fixtures: 18 passed, 0 failed. This exercises the assertions, but does not replace a conventional pytest invocation.

The saved-table derivation was run twice; SHA-256 outputs for `spatial_robustness.tsv` and `genealogy_correspondence.tsv` were identical across runs. The full 988-window raw-genotype scan was not rerun.
