"""Guard against promoting mixed, failed, or sensitive runs to an accepted study."""
import json, sys, tempfile, unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / 'scripts'))
from audit_hv_stabilization import audit, METRICS
from run_hv_sweep import sha


class StabilizationAuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name)
        mesh = root / 'mesh'; mesh.write_text('synthetic test fixture')
        exe = root / 'exe'; exe.write_text('synthetic executable fixture')
        self.runs = []
        for i, scale in enumerate([.5, 1., 2.]):
            d = root / f'run{i}'; d.mkdir(); self.runs.append(d)
            (d / 'nodes.csv').write_text('synthetic node fixture')
            m = dict(executable=str(exe), executable_sha256=sha(exe), mesh=str(mesh),
                     mesh_sha256=sha(mesh), model='FROZEN_DRIFT_HEAT_ZNO_NT', supg_scale=scale)
            (d / 'manifest.json').write_text(json.dumps(m))
            (d / 'execution.json').write_text(json.dumps(dict(exit_code=0)))
            (d / 'field_audit.json').write_text(json.dumps(dict(passed=True)))
            data = dict(final_voltage_factor=1, final_joule_factor=1,
                        maximum_accepted_residual=1e-13, maximum_accepted_update=1e-12,
                        case_acceptance='PASS_NOT_MESH_ACCEPTANCE', voltage_V=6,
                        maximum_stage_balance=1e-12, heat_balance_relative=1e-4,
                        electron_balance_relative=1e-4, hole_balance_diagnostic=1e-4,
                        joule_quadrature_4_5_relative=1e-6,
                        model='NT', heat_source='DRIFT_CONDUCTIVITY_NOT_TOTAL_CURRENT_WORK',
                        source='FROZEN_CANDIDATE', supg_directional_tau='T', supg_scale=scale)
            data.update({k: 1 + i * .002 for k in METRICS})
            (d / 'metrics.txt').write_text('\n'.join(f'{k}={v}' for k, v in data.items()))

    def test_small_variation_is_not_full_validation(self):
        result = audit(self.runs)
        self.assertTrue(result['covered_metrics_passed'])
        self.assertFalse(result['fully_validated'])

    def test_large_hole_metric_variation_fails(self):
        path = self.runs[2] / 'metrics.txt'
        path.write_text(path.read_text().replace('ring_abs_p_relative_mean=1.004',
                                                'ring_abs_p_relative_mean=1.04'))
        self.assertFalse(audit(self.runs)['covered_metrics_passed'])

    def test_failed_flux_and_missing_audit_do_not_pass(self):
        path = self.runs[0] / 'metrics.txt'
        path.write_text(path.read_text().replace('PASS_NOT_MESH_ACCEPTANCE', 'FAIL'))
        (self.runs[1] / 'field_audit.json').unlink()
        self.assertFalse(audit(self.runs)['covered_metrics_passed'])

    def test_wrong_voltage_and_stale_mesh_hash_fail(self):
        path = self.runs[0] / 'metrics.txt'
        path.write_text(path.read_text().replace('voltage_V=6', 'voltage_V=5'))
        manifest = json.loads((self.runs[0] / 'manifest.json').read_text())
        Path(manifest['mesh']).write_text('changed fixture')
        result = audit(self.runs)
        self.assertFalse(result['covered_metrics_passed'])
        self.assertIn('INCONSISTENT_VOLTAGE_V', result['errors'])

    def test_null_coefficient_is_rejected(self):
        path = self.runs[0] / 'manifest.json'
        data = json.loads(path.read_text()); data['supg_scale'] = None
        path.write_text(json.dumps(data))
        result = audit(self.runs)
        self.assertFalse(result['covered_metrics_passed'])
        self.assertIn('SCALES_MUST_BE_0.5_1_2', result['errors'])


if __name__ == '__main__':
    unittest.main()
