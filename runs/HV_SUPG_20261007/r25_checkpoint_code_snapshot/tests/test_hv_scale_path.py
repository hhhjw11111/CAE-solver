"""A failed coefficient solve must roll back without promoting its field."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).absolute().parents[1]/'scripts'))
from run_hv_scale_path import ScalePath


class CoefficientRollback(unittest.TestCase):
    def test_failed_trial_keeps_last_accepted_coefficient(self):
        path = ScalePath(1., 2.)
        self.assertEqual(path.trial(), 1.25)
        self.assertTrue(path.reject())
        self.assertEqual(path.current, 1.)
        self.assertEqual(path.trial(), 1.125)
        path.accept(1.125)
        self.assertEqual(path.current, 1.125)
        self.assertEqual(path.trial(), 1.3125)

    def test_descending_path_and_target_are_actual_steps(self):
        path = ScalePath(1., .5)
        path.accept(.75)
        self.assertFalse(path.complete)
        path.accept(.5)
        self.assertTrue(path.complete)

    def test_short_final_failure_does_not_repeat_same_trial(self):
        path = ScalePath(1.98, 2.)
        self.assertFalse(path.reject())
        self.assertEqual(path.current, 1.98)
        self.assertFalse(path.complete)

    def test_unrequested_or_nonfinite_states_are_rejected(self):
        path = ScalePath(1., 2.)
        with self.assertRaises(ValueError):
            path.accept(2.)
        self.assertEqual(path.current, 1.)
        for target in [0, -1, float('nan'), float('inf')]:
            with self.assertRaises(ValueError):
                ScalePath(1., target)


if __name__ == '__main__':
    unittest.main()
