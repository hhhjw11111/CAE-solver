"""Numerical homotopy must not weaken final physics or default controls."""
from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).absolute().parents[1]/'scripts'))
from run_hv_point import apply_load_path
from run_hv_sweep import make_input

class LoadPath(unittest.TestCase):
    def test_default_input_is_identical(self):
        original=make_input(Path('mesh'),8.,Path('seed'),6.,1.,False)
        self.assertEqual(apply_load_path(original),original)

    def test_seed_is_reequilibrated_without_weakening_target(self):
        text=apply_load_path(make_input(Path('mesh'),8.,Path('seed'),6.,1.,False),True)
        for token in ['restart_voltage_factor=0.75','restart_joule_factor=0.0,',
                      'target_voltage=-8,','target_voltage_factor=1.0,target_joule_factor=1.0',
                      'residual_tolerance=1e-9,update_tolerance=1e-9',
                      'continuation_full_joule_during_voltage=.false.,']:
            self.assertIn(token,text)
        self.assertNotIn('restart_joule_factor=1.0,',text)

    def test_zero_start_uses_existing_native_zero_heat_default(self):
        text=apply_load_path(make_input(Path('mesh'),8.,None,0.,1.,False),True)
        self.assertNotIn('restart_nodes_file',text)
        self.assertIn('target_joule_factor=1.0',text)
        self.assertIn('continuation_full_joule_during_voltage=.false.,',text)

if __name__=='__main__':unittest.main()
