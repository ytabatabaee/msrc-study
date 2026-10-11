import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from analyze_mai_rf import MAIN_CONDITIONS, get_msrc_collections  # noqa: E402


def test_msrc_matching_preserves_leading_zero_replicate_ids():
    collections = get_msrc_collections()
    key = (MAIN_CONDITIONS[1], "007", "msrc_no_suppression")
    assert key in collections
    assert (MAIN_CONDITIONS[1], "7", "msrc_no_suppression") not in collections


def test_missing_conditions_are_not_substituted_by_other_conditions():
    collections = get_msrc_collections()
    assert not any(key[0] == MAIN_CONDITIONS[0] for key in collections)
