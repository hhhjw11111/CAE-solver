"""Guard against promoting mixed, failed, or sensitive runs to an accepted study."""
import json, sys, tempfile, unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / 'scripts'))
from audit_hv_quadrature import audit, METRICS
from run_hv_sweep import sha


class QuadratureAuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name)
        mesh = root / 'mesh'; mesh.write_text('synthetic test fixture')
        exe = root / 'exe'; exe.write_text('synthetic executable fixture')
        self.runs = []
        for i, scale in enumerate([3,5]):
            d = root / f'run{i}'; d.mkdir(); self.runs.append(d)
            (d / 'nodes.csv').write_text('synthetic node fixture')
            m = dict(executable=str(exe), executable_sha256=sha(exe), mesh=str(mesh),
                     mesh_sha256=sha(mesh), model='FROZEN_DRIFT_HEAT_ZNO_NT', supg_scale=1.,log_weighted_quadrature_order=scale)
            (d / 'manifest.json').write_text(json.dumps(m))
            (d/'input.nml').write_text('synthetic input fixture')
            (d / 'execution.json').write_text(json.dumps(dict(exit_code=0,input_sha256=sha(d/'input.nml'),executable_sha256=sha(exe))))
            (d / 'field_audit.json').write_text(json.dumps(dict(passed=True)))
            data = dict(final_voltage_factor=1, final_joule_factor=1,
                        maximum_accepted_residual=1e-13, maximum_accepted_update=1e-12,
                        case_acceptance='PASS_NOT_MESH_ACCEPTANCE', voltage_V=6,
                        maximum_stage_balance=1e-12, heat_balance_relative=1e-4,
                        electron_balance_relative=1e-4, hole_balance_diagnostic=1e-4,
                        joule_quadrature_4_5_relative=1e-6,
                        model='NT', heat_source='DRIFT_CONDUCTIVITY_NOT_TOTAL_CURRENT_WORK',
                        source='FROZEN_CANDIDATE', supg_directional_tau='T', supg_scale=1.,log_weighted_quadrature_order=scale,
                        transport_discretization='HEX20_LOG_WEIGHTED_EXPERIMENTAL')
            data.update({k: 1 + i * .002 for k in METRICS})
            (d / 'metrics.txt').write_text('\n'.join(f'{k}={v}' for k, v in data.items()))

    def test_small_variation_is_not_full_validation(self):
        result=audit(self.runs)
        self.assertTrue(result['covered_metrics_passed'])
        self.assertFalse(result['fully_validated'])

    def test_hole_sensitivity_and_same_order_fail(self):
        path=self.runs[1]/'metrics.txt'
        original=path.read_text()
        path.write_text(original.replace('ring_abs_p_relative_mean=1.002','ring_abs_p_relative_mean=1.04'))
        self.assertFalse(audit(self.runs)['covered_metrics_passed'])
        path.write_text(original.replace('log_weighted_quadrature_order=5','log_weighted_quadrature_order=3'))
        manifest=self.runs[1]/'manifest.json'
        data=json.loads(manifest.read_text());data['log_weighted_quadrature_order']=3
        manifest.write_text(json.dumps(data))
        result=audit(self.runs)
        self.assertIn('ORDERS_MUST_BE_3_5',result['errors'])

    def test_changed_input_or_wrong_method_fail(self):
        (self.runs[0]/'input.nml').write_text('modified after execution')
        path=self.runs[1]/'metrics.txt'
        path.write_text(path.read_text().replace('HEX20_LOG_WEIGHTED_EXPERIMENTAL','HEX20_CONSISTENT_SUPG_EXPERIMENTAL'))
        result=audit(self.runs)
        self.assertFalse(result['covered_metrics_passed'])
        self.assertTrue(any('INPUT_HASH_MISMATCH' in e for e in result['errors']))
        self.assertTrue(any('QUADRATURE_STUDY_REQUIRES_LOG_WEIGHTED' in e for e in result['errors']))

    def test_local_carrier_difference_is_not_hidden_by_integrals(self):
        header='node,x_m,y_m,z_m,ux_m,uy_m,uz_m,phi_V,eta_n,eta_p,theta_K\n'
        for i,d in enumerate(self.runs):
            eta=1e-4*(1+.04*i)
            (d/'nodes.csv').write_text(header+f'1,0,0,0,0,0,0,1,{eta},0,1\n')
        result=audit(self.runs,include_local_fields=True)
        self.assertTrue(result['covered_metrics_passed'])
        self.assertFalse(result['all_requested_checks_passed'])
        self.assertFalse(result['local_field_comparisons']['eta_n_density_perturbation']['passed'])

if __name__=='__main__':unittest.main()
