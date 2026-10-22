"""Checkpoint replay must not invent missing voltages or accept incomplete states."""
import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).absolute().parents[1]/'scripts'))
from run_hv_sweep import checkpoint_seed

class CheckpointReplay(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)

    def write(self,index=1,voltage=-1.0000000000000002,joule=1,residual=1e-13,nodes=True):
        load=self.root/f'accepted_{index:05d}_load.txt'
        node=self.root/f'accepted_{index:05d}_nodes.csv'
        load.write_text(f'status=ALGEBRAIC_CHECKPOINT_NOT_FLUX_OR_MESH_ACCEPTANCE\napplied_voltage_V={voltage}\njoule_factor={joule}\nresidual={residual}\nupdate=1e-12\n')
        if nodes:node.write_text('synthetic initial-field fixture')
        return node,load

    def test_roundoff_equivalent_voltage_can_be_replayed(self):
        expected=self.write()
        self.assertEqual(checkpoint_seed(self.root,1.),expected)

    def test_missing_voltage_is_not_interpolated(self):
        self.write(voltage=-1.)
        self.write(index=2,voltage=-2.)
        with self.assertRaisesRegex(ValueError,'interpolation is not allowed'):
            checkpoint_seed(self.root,1.5)

    def test_nonfinite_or_missing_state_is_rejected(self):
        self.write(residual=float('nan'))
        with self.assertRaises(ValueError):checkpoint_seed(self.root,1.)
        self.write(residual=1e-13)
        (self.root/'accepted_00001_nodes.csv').unlink()
        with self.assertRaises(ValueError):checkpoint_seed(self.root,1.)

    def test_partial_heat_is_not_the_full_heat_checkpoint(self):
        self.write(joule=0.)
        expected=self.write(index=2,joule=1.)
        self.assertEqual(checkpoint_seed(self.root,1.),expected)

if __name__=='__main__':unittest.main()
