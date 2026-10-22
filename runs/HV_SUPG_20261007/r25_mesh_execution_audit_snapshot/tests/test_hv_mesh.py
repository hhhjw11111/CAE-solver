"""A converged temperature must not mask unconverged holes or mismatched runs."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).absolute().parents[1]/'scripts'))
from audit_hv_mesh import audit, LIMITS
from run_hv_sweep import sha


class MeshGate(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name)
        exe = root/'exe'; exe.write_text('synthetic executable fixture')
        self.levels = []
        for i in range(3):
            d = root/f'L{i}'; d.mkdir(); self.levels.append(d)
            mesh = d/'mesh'; mesh.write_text(f'synthetic mesh fixture level {i}')
            manifest = dict(mesh=str(mesh), mesh_sha256=sha(mesh), executable=str(exe),
                            executable_sha256=sha(exe), supg_scale=1., model='FROZEN_DRIFT_HEAT_ZNO_NT')
            (d/'manifest.json').write_text(json.dumps(manifest))
            (d/'input.nml').write_text('synthetic input fixture')
            (d/'execution.json').write_text(json.dumps(dict(exit_code=0,
                executable_sha256=sha(exe),input_sha256=sha(d/'input.nml'))))
            (d/'field_audit.json').write_text(json.dumps(dict(passed=True)))
            (d/'nodes.csv').write_text('synthetic test fixture')
            metrics = dict(final_voltage_factor=1, final_joule_factor=1,
                           maximum_accepted_residual=1e-12, maximum_accepted_update=1e-12,
                           maximum_stage_balance=1e-12, heat_balance_relative=1e-4,
                           electron_balance_relative=1e-4, hole_balance_diagnostic=1e-4,
                           joule_quadrature_4_5_relative=1e-6,
                           case_acceptance='PASS_NOT_MESH_ACCEPTANCE', voltage_V=6., supg_scale=1.,
                           model='NT', heat_source='DRIFT_CONDUCTIVITY_NOT_TOTAL_CURRENT_WORK',
                           source='FROZEN_CANDIDATE', supg_directional_tau='T',
                           transport_discretization='HEX20_CONSISTENT_SUPG_EXPERIMENTAL')
            metrics.update({key: 1. + i*.002 for key in LIMITS})
            (d/'metrics.txt').write_text('\n'.join(f'{key}={value}' for key,value in metrics.items()))

    def test_mixed_quadrature_cannot_be_a_mesh_study(self):
        path=self.levels[2]/'manifest.json'
        manifest=json.loads(path.read_text());manifest['log_weighted_quadrature_order']=5
        path.write_text(json.dumps(manifest))
        with (self.levels[2]/'metrics.txt').open('a') as f:
            f.write('\nlog_weighted_quadrature_order=5')
        result=audit(self.levels,[6.])
        self.assertFalse(result['all_mesh_pass'])
        self.assertIn('INCONSISTENT_QUADRATURE_ORDER',result['provenance_errors'])

    def test_changed_run_input_cannot_pass(self):
        (self.levels[1]/'input.nml').write_text('changed after execution')
        result=audit(self.levels,[6.])
        self.assertFalse(result['all_mesh_pass'])
        self.assertEqual(result['voltages'][0]['status'],'EXECUTION_INPUT_MISMATCH')

    def test_run_executable_must_match_study_manifest(self):
        path=self.levels[1]/'execution.json'
        record=json.loads(path.read_text());record['executable_sha256']='different-build'
        path.write_text(json.dumps(record))
        result=audit(self.levels,[6.])
        self.assertFalse(result['all_mesh_pass'])
        self.assertEqual(result['voltages'][0]['status'],'EXECUTION_EXECUTABLE_MISMATCH')

    def test_holes_must_converge_at_frozen_two_percent_limit(self):
        self.assertTrue(audit(self.levels, [6.])['all_mesh_pass'])
        path = self.levels[2]/'metrics.txt'
        path.write_text(path.read_text().replace('ring_abs_p_relative_mean=1.004',
                                                'ring_abs_p_relative_mean=1.04'))
        result = audit(self.levels, [6.])
        self.assertFalse(result['all_mesh_pass'])
        self.assertEqual(result['voltages'][0]['status'], 'MESH_FAIL')
        self.assertGreater(result['voltages'][0]['changes']['ring_abs_p_relative_mean']['L1_L2'], .02)

    def test_missing_holes_and_wrong_voltage_cannot_pass(self):
        path = self.levels[2]/'metrics.txt'
        original = path.read_text()
        for altered in [original.replace('ring_abs_p_relative_mean=1.004', ''),
                        original.replace('voltage_V=6.0', 'voltage_V=5.0')]:
            path.write_text(altered)
            self.assertFalse(audit(self.levels, [6.])['all_mesh_pass'])

    def test_changed_executable_and_empty_study_cannot_pass(self):
        self.assertFalse(audit(self.levels, [])['all_mesh_pass'])
        manifest = json.loads((self.levels[0]/'manifest.json').read_text())
        Path(manifest['executable']).write_text('changed executable fixture')
        result = audit(self.levels, [6.])
        self.assertFalse(result['all_mesh_pass'])
        self.assertIn('EXECUTABLE_HASH_MISMATCH', result['provenance_errors'])

    def test_mixed_domain_integrals_cannot_pass(self):
        path=self.levels[2]/'manifest.json';m=json.loads(path.read_text())
        m['domain_symmetry']='Y_HALF';path.write_text(json.dumps(m))
        self.assertFalse(audit(self.levels,[6.])['all_mesh_pass'])

    def test_mixed_transport_weak_forms_cannot_pass(self):
        path=self.levels[2]/'metrics.txt'
        path.write_text(path.read_text().replace('HEX20_CONSISTENT_SUPG_EXPERIMENTAL',
                                                 'HEX20_LOG_WEIGHTED_EXPERIMENTAL'))
        result=audit(self.levels,[6.])
        self.assertFalse(result['all_mesh_pass'])
        self.assertIn('INCONSISTENT_TRANSPORT_DISCRETIZATION',result['provenance_errors'])


if __name__ == '__main__':
    unittest.main()
