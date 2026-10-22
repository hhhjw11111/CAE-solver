from pathlib import Path
import re,sys,unittest
sys.path.insert(0,str(Path(__file__).absolute().parents[1]/'scripts'))
from run_hv_point import apply_load_path,apply_heat_restart
from run_hv_sweep import make_input

def fixture():return apply_load_path(make_input(Path('mesh'),8.,Path('seed'),8.,1.,False),True)

class HeatRestart(unittest.TestCase):
    def test_defaults_leave_input_identical(self):
        text=fixture();self.assertEqual(apply_heat_restart(text),text)

    def test_partial_heat_restart_preserves_target_and_tolerances(self):
        old=fixture();new=apply_heat_restart(old,.03125,.125,.25)
        self.assertIn('restart_joule_factor=0.03125,',new)
        self.assertIn('continuation_initial_step=0.125,',new)
        self.assertIn('continuation_maximum_step=0.25,',new)
        for token in ['target_voltage_factor=1.0,target_joule_factor=1.0',
                      'residual_tolerance=1e-9,update_tolerance=1e-9']:
            self.assertIn(token,new)
        for key in ['continuation_minimum_step','continuation_maximum_attempts','max_newton']:
            self.assertEqual(re.search(key+r'=([^,]+)',old)[1],re.search(key+r'=([^,]+)',new)[1])

    def test_invalid_controls_rejected(self):
        for controls in [(2.,None,None),(float('nan'),None,None),(.1,0.,.25),(.1,.5,.25),
                         (.1,.125,None),(.1,1e-5,.25),(.1,.125,float('inf'))]:
            with self.assertRaises(ValueError):apply_heat_restart(fixture(),*controls)
        with self.assertRaises(ValueError):
            apply_heat_restart(make_input(Path('mesh'),8.,Path('seed'),8.,1.,False),.1,.125,.25)

if __name__=='__main__':unittest.main()
