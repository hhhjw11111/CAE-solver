import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).absolute().parents[1]/'scripts'))
from prepare_hv_zero_heat_voltage_guess import voltage_guess

def fixture():
    row=dict(node='1',x_m='0',y_m='0',z_m='0',phi_V='-6',ux_m='1e-10',uy_m='0',uz_m='0',
             theta_K='0',eta_n='0',eta_p='0',n_m3='2',p_m3='3')
    return [row,dict(row,node='2',phi_V='0')],dict(applied_voltage_V=-6.,joule_factor=0.,residual=1e-16,update=1e-11)

class VoltageGuess(unittest.TestCase):
    def test_scale_is_only_a_guess_preserving_geometry_background_and_source(self):
        rows,load=fixture();before=[dict(r) for r in rows]
        result,backgrounds,ratio=voltage_guess(rows,load,8.)
        self.assertEqual(rows,before);self.assertEqual(ratio,8/6)
        self.assertEqual(float(result[0]['phi_V']),-8.)
        self.assertAlmostEqual(float(result[0]['ux_m']),1e-10*ratio,places=20)
        self.assertEqual(backgrounds,[2.,3.])
        self.assertEqual(result[0]['x_m'],rows[0]['x_m'])

    def test_nonzero_heat_or_unconverged_checkpoint_is_rejected(self):
        rows,load=fixture()
        for change in [dict(joule_factor=1.),dict(residual=2e-9),dict(update=float('nan'))]:
            with self.assertRaises(ValueError):voltage_guess(rows,dict(load,**change),8.)

    def test_nonuniform_heated_or_wrong_voltage_source_rejected(self):
        for change in [dict(theta_K='1'),dict(eta_n='.1'),dict(phi_V='-5')]:
            rows,load=fixture();rows[0].update(change)
            with self.assertRaises(ValueError):voltage_guess(rows,load,8.)

if __name__=='__main__':unittest.main()
