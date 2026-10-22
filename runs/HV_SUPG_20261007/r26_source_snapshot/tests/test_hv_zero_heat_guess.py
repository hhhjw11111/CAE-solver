import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).absolute().parents[1]/'scripts'))
from prepare_hv_zero_heat_guess import zero_heat_rows,RESET

def fixture():
    return dict(node='1',x_m='0',y_m='0',z_m='0',phi_V='-6',ux_m='1e-10',uy_m='0',uz_m='0',
                theta_K='20',eta_n='0.2',eta_p='-0.3',n_m3=str(2*math.exp(.2)),p_m3=str(3*math.exp(-.3)))

class ZeroHeatGuess(unittest.TestCase):
    def test_reset_is_only_an_initial_guess_without_mutating_source(self):
        source=fixture();before=dict(source);rows,backgrounds=zero_heat_rows([source])
        self.assertEqual(source,before)
        self.assertEqual(rows[0]['phi_V'],'-6')
        for key in RESET:self.assertEqual(float(rows[0][key]),0.)
        self.assertAlmostEqual(backgrounds[0],2.);self.assertAlmostEqual(backgrounds[1],3.)
        self.assertEqual(float(rows[0]['n_m3']),backgrounds[0])

    def test_nonfinite_or_nonpositive_density_is_rejected(self):
        for change in [{'theta_K':'nan'},{'n_m3':'0'},{'p_m3':'-1'},{'node':'2'}]:
            with self.assertRaises(ValueError):zero_heat_rows([dict(fixture(),**change)])

    def test_nonuniform_background_is_rejected(self):
        a=fixture();b=dict(a,node='2',n_m3='100')
        with self.assertRaises(ValueError):zero_heat_rows([a,b])

if __name__=='__main__':unittest.main()
