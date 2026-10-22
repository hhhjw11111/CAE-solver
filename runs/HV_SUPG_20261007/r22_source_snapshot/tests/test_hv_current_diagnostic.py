"""Check two analytic limits of the read-only current reconstruction."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).absolute().parents[1]/'scripts'))
from diagnose_hv_drift_diffusion_cancellation import indicators, KB, Q


class CurrentLimits(unittest.TestCase):
    def test_uniform_density_recovers_drift_power(self):
        e = np.array([[1e5, -2e5, 3e5]])
        zero = np.zeros_like(e)
        heat, _, dissipation, work, ratio = indicators(e, zero, zero, np.array([320.]), np.array([6e22]), np.array([1e-30]))
        np.testing.assert_allclose(dissipation, heat, rtol=1e-14)
        np.testing.assert_allclose(work, heat, rtol=1e-14)
        np.testing.assert_allclose(ratio, 1., rtol=1e-14)

    def test_local_boltzmann_balance_has_zero_current_but_nonzero_frozen_heat(self):
        e = np.array([[0., 0., 1e6]])
        t = np.array([300.]); vt = KB*t/Q
        heat, _, dissipation, work, ratio = indicators(e, -e/vt[:,None], e/vt[:,None], t, np.array([6e22]), np.array([1e-30]))
        self.assertGreater(float(heat[0]), 0.)
        self.assertLess(float(dissipation[0]/heat[0]), 1e-28)
        self.assertLess(abs(float(work[0]/heat[0])), 1e-14)
        self.assertLess(float(ratio[0]), 1e-14)


if __name__ == '__main__':
    unittest.main()
