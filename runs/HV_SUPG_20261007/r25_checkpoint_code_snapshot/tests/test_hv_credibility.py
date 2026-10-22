"""An evidence ledger must not promote missing, stale or mislabelled data."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).absolute().parents[1] / 'scripts'))
from run_hv_sweep import sha
from summarize_hv_credibility import point_record, summarize, SOURCE


class CredibilityLedger(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.mesh = self.root / 'mesh'
        self.mesh.write_text('2 0\n1 -1e-5 0 0\n2 1e-5 1e-5 0\n')
        self.exe = self.root / 'exe'
        self.exe.write_text('test executable fixture, not a native solve')
        (self.root / 'input.nml').write_text('synthetic test input')
        self.manifest = dict(mesh=str(self.mesh), mesh_sha256=sha(self.mesh),
                             executable=str(self.exe), executable_sha256=sha(self.exe),
                             model='FROZEN_DRIFT_HEAT_ZNO_NT', supg_scale=1.)
        self.metrics = dict(model='NT', source=SOURCE, heat_source='DRIFT_CONDUCTIVITY_NOT_TOTAL_CURRENT_WORK',
                            voltage_V=0., supg_scale=1., final_voltage_factor=1., final_joule_factor=1.,
                            maximum_accepted_residual=0., maximum_accepted_update=0.,
                            maximum_stage_balance=0., heat_balance_relative=0.,
                            electron_balance_relative=0., hole_balance_diagnostic=0.,
                            joule_quadrature_4_5_relative=0., theta_max_K=0., joule_W=0.,
                            electron_right_A=0., hole_right_A=0., ring_abs_n_relative_mean=0.,
                            ring_abs_p_relative_mean=0., case_acceptance='PASS_NOT_MESH_ACCEPTANCE')
        self.write_records()

    def write_records(self, voltage=0.):
        (self.root / 'manifest.json').write_text(json.dumps(self.manifest))
        execution = dict(exit_code=0, timed_out=False, executable_sha256=sha(self.exe),
                         input_sha256=sha(self.root / 'input.nml'))
        (self.root / 'execution.json').write_text(json.dumps(execution))
        self.metrics['voltage_V'] = voltage
        (self.root / 'metrics.txt').write_text('\n'.join(f'{k}={v}' for k,v in self.metrics.items()))
        columns = ['node','x_m','y_m','z_m','ux_m','uy_m','uz_m','phi_V','eta_n','eta_p','theta_K','n_m3','p_m3']
        with (self.root / 'nodes.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(columns)
            writer.writerow([1,-1e-5,0,0,0,0,0,0,0,0,0,6e22,(2.9e-4)**2/6e22])
            writer.writerow([2,1e-5,1e-5,0,0,0,0,-voltage,0,0,0,6e22,(2.9e-4)**2/6e22])

    def record(self, voltage=0.):
        return point_record(voltage, self.root, self.root / 'manifest.json')

    def test_zero_needs_actual_exact_zero_field(self):
        self.assertTrue(self.record()['exact_zero_equilibrium_verified'])
        path = self.root / 'nodes.csv'
        with path.open() as stream:
            rows = list(csv.DictReader(stream))
        rows[0]['eta_p'] = '1e-15'
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
            writer.writeheader(); writer.writerows(rows)
        self.assertFalse(self.record()['fully_validated'])

    def test_transport_and_quadrature_identity_cannot_be_relabelled(self):
        self.manifest['transport_discretization'] = 'HEX20_LOG_WEIGHTED_EXPERIMENTAL'
        self.write_records()
        self.assertIn('TRANSPORT_DISCRETIZATION_MISMATCH', self.record()['evidence_errors'])
        self.metrics['transport_discretization'] = self.manifest['transport_discretization']
        self.manifest['log_weighted_quadrature_order'] = 5
        self.write_records()
        self.assertIn('QUADRATURE_ORDER_MISMATCH', self.record()['evidence_errors'])
        self.metrics['log_weighted_quadrature_order'] = 5
        self.write_records()
        self.assertTrue(self.record()['fully_validated'])

    def test_changed_input_or_execution_identity_is_invalid(self):
        (self.root / 'input.nml').write_text('changed after recorded execution')
        self.assertIn('INPUT_HASH_MISMATCH', self.record()['evidence_errors'])
        self.write_records()
        path = self.root / 'execution.json'
        data = json.loads(path.read_text()); data['executable_sha256'] = 'different'
        path.write_text(json.dumps(data))
        self.assertEqual(self.record()['state'], 'INVALID_POINT_EVIDENCE')

    def test_nonzero_pass_and_wrong_voltage_do_not_certify(self):
        self.write_records(6.)
        self.assertEqual(self.record(6.)['state'], 'SOLVE_FLUX_PASS_MESH_PENDING')
        self.assertFalse(self.record(6.)['fully_validated'])
        self.assertIn('VOLTAGE_MISMATCH', self.record(8.)['evidence_errors'])

    def test_half_integrals_are_explicit_and_temperature_is_not_doubled(self):
        self.manifest['domain_symmetry'] = 'Y_HALF'
        self.metrics['domain_symmetry'] = 'Y_HALF'
        for key in ['physical_exterior_heat_balance_relative', 'physical_exterior_electron_balance_relative',
                    'physical_exterior_hole_balance_diagnostic', 'symmetry_heat_flux_relative',
                    'symmetry_electron_flux_relative', 'symmetry_hole_flux_diagnostic']:
            self.metrics[key] = 0.
        self.metrics.update(joule_W=.004, electron_right_A=.0004, theta_max_K=20.)
        self.write_records(6.)
        record = self.record(6.)
        self.assertEqual(record['state'], 'SOLVE_FLUX_PASS_MESH_PENDING')
        self.assertEqual(record['metrics']['joule_W'], .004)
        self.assertEqual(record['reconstructed_full_integrals']['joule_W'], .008)
        self.assertEqual(record['metrics']['theta_max_K'], 20.)

    def test_missing_grid_points_stay_missing(self):
        plan = self.root / 'plan.json'
        plan.write_text(json.dumps(dict(points=[dict(voltage_V=0., directory='.')], studies=[])))
        result = summarize(plan)
        self.assertEqual(len(result['points']), 17)
        self.assertEqual(result['verified_exact_voltages_V'], [0.])
        self.assertEqual(result['points'][-1]['state'], 'NO_REGISTERED_POINT')
        self.assertFalse(result['complete_quantitative_0_8V'])
        self.assertIsNone(result['certified_nonzero_interval'])

    def test_off_grid_branch_sample_does_not_fill_next_grid_voltage(self):
        voltage = 6.459300885488345
        self.write_records(voltage)
        plan = self.root / 'plan.json'
        plan.write_text(json.dumps(dict(points=[], studies=[],
                                        off_grid_points=[dict(voltage_V=voltage, directory='.')])) )
        result = summarize(plan)
        self.assertEqual(result['points'][13]['voltage_V'], 6.5)
        self.assertEqual(result['points'][13]['state'], 'NO_REGISTERED_POINT')
        self.assertEqual(result['off_grid_points'][0]['state'], 'SOLVE_FLUX_PASS_MESH_PENDING')
        self.assertFalse(result['off_grid_points'][0]['fully_validated'])


if __name__ == '__main__':
    unittest.main()
