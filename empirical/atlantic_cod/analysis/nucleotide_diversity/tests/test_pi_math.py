import importlib.util
from pathlib import Path
import unittest
SCRIPT=Path(__file__).resolve().parents[1]/'scripts/estimate_tracks.py'
spec=importlib.util.spec_from_file_location('track_estimator',SCRIPT)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
class PiMathTests(unittest.TestCase):
    def test_finite_sample_correction_manual_example(self):
        # Three ref copies and one alt: raw H=.375; correction 4/3 gives .5.
        self.assertAlmostEqual(mod.corrected_pi(1,4),.5)
    def test_single_copy_is_missing(self):
        self.assertIsNone(mod.corrected_pi(1,1))
if __name__=='__main__':unittest.main()
