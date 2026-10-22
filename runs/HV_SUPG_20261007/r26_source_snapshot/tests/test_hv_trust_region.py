"""The optional globalization preserves load targets, equations and old inputs."""
from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).absolute().parents[1]/'scripts'))
from run_hv_point import apply_trust_region,apply_load_path
from run_hv_sweep import make_input

class TrustRegion(unittest.TestCase):
    def test_default_is_byte_identical(self):
        text=make_input(Path('mesh'),8.,Path('seed'),6.,1.,False)
        self.assertEqual(apply_trust_region(text),text)

    def test_both_paths_retain_final_physics_and_gates(self):
        for cold in (False,True):
            old=apply_load_path(make_input(Path('mesh'),8.,Path('seed'),6.,1.,False),cold)
            new=apply_trust_region(old,True,.125)
            self.assertEqual(new.replace('\n newton_trust_region=.true.,trust_initial_radius=0.125,',''),old)
            self.assertIn('target_voltage_factor=1.0,target_joule_factor=1.0',new)
            self.assertIn('residual_tolerance=1e-9,update_tolerance=1e-9',new)

    def test_invalid_radius_rejected(self):
        for radius in (-1.,0.,float('nan'),float('inf'),5.):
            with self.assertRaises(ValueError):apply_trust_region('&solver_control',True,radius)

if __name__=='__main__':unittest.main()
