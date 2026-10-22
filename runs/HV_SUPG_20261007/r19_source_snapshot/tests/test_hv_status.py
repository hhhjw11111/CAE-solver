"""Ensure incomplete/partial/NaN runs cannot become accepted NT points."""
import sys, tempfile, unittest,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_hv_sweep import classify,make_input,sha
from audit_hv_mesh import audit
from run_hv_trajectory import checkpoint_voltage

class StatusGate(unittest.TestCase):
    def test_missing_mesh_evidence_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            r=audit([Path(d)/f'L{i}' for i in range(3)],[0,8])
            self.assertFalse(r['all_mesh_pass'])
            self.assertTrue(all(v['status']=='MISSING_EVIDENCE' for v in r['voltages']))
    def test_complete_means_mesh_pending(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'nodes.csv').write_text('placeholder')
            valid=dict(final_voltage_factor=1,final_joule_factor=1,maximum_accepted_residual=1e-12,
                maximum_accepted_update=1e-12,theta_max_K=20,joule_W=.01,case_acceptance='PASS_NOT_MESH_ACCEPTANCE',
                maximum_stage_balance=1e-12,heat_balance_relative=1e-4,electron_balance_relative=1e-4,
                hole_balance_diagnostic=1e-4,joule_quadrature_4_5_relative=1e-6)
            def write(m): (p/'metrics.txt').write_text(''.join(f'{k}={v}\n' for k,v in m.items()))
            write(valid)
            self.assertEqual(classify(p,0)[0],'SOLVE_FLUX_PASS_MESH_PENDING')
            for k,v in [('final_joule_factor',.99),('final_voltage_factor',.99),
                        ('maximum_accepted_residual',1e-7),('maximum_accepted_update',1e-7),('theta_max_K','NaN')]:
                write(dict(valid,**{k:v}))
                self.assertEqual(classify(p,0)[0],'FAILED_NO_COMPLETE_STATE')
            write(dict(valid,case_acceptance='FAIL'))
            self.assertEqual(classify(p,50)[0],'SOLVED_FLUX_GATE_FAILED')
            for k in ['heat_balance_relative','electron_balance_relative','hole_balance_diagnostic',
                      'joule_quadrature_4_5_relative','maximum_stage_balance']:
                for value in ['NaN',float('inf'),-1.,.02]:
                    write(dict(valid,**{k:value}))
                    self.assertEqual(classify(p,0)[0],'SOLVED_FLUX_GATE_FAILED')
                missing=dict(valid);del missing[k];write(missing)
                self.assertEqual(classify(p,0)[0],'SOLVED_FLUX_GATE_FAILED')
    def test_full_heat_restart_path(self):
        t=make_input(Path('mesh.mesh'),8,Path('previous'),7.5,1,True)
        self.assertIn('restart_voltage_factor=0.9375',t)
        self.assertIn('restart_joule_factor=1.0',t)
        self.assertIn('continuation_full_joule_during_voltage=.true.',t)

    def test_copied_mesh_is_not_three_level_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            roots=[Path(d)/f'L{i}' for i in range(3)]
            mesh=Path(d)/'same.mesh';mesh.write_text('synthetic identical mesh')
            for root in roots:
                root.mkdir()
                (root/'manifest.json').write_text(json.dumps(dict(mesh=str(mesh),mesh_sha256=sha(mesh),
                    executable_sha256='same-executable',supg_scale=1,model='same-model')))
            r=audit(roots,[8])
            self.assertIn('THREE_DISTINCT_MESHES_REQUIRED',r['provenance_errors'])
            self.assertFalse(r['all_mesh_pass'])

    def test_partial_or_invalid_checkpoint_is_never_a_grid_point(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'load.txt'
            good=dict(applied_voltage_V=-6.5,joule_factor=1,residual=1e-12,update=1e-12)
            def write(m):p.write_text(''.join(f'{k}={v}\n' for k,v in m.items()))
            write(good);self.assertEqual(checkpoint_voltage(p),6.5)
            for change in [dict(joule_factor=.9),dict(residual='NaN'),dict(update=1e-5),dict(applied_voltage_V=-6.375)]:
                write(dict(good,**change));self.assertIsNone(checkpoint_voltage(p))
            p.write_text('applied_voltage_V=-6.5\n');self.assertIsNone(checkpoint_voltage(p))

if __name__=='__main__':unittest.main()
