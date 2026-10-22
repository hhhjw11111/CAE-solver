"""Rebuild a voltage evidence ledger from native runs without certifying an interval.

The plan lists actual point directories and manifests, not interpolated values.
This tool rechecks input/executable/mesh provenance, native numeric gates and
nodal fields. Study results are recomputed from their source runs. Their scope
stays separate: failure on one mesh is not a result for every mesh at that voltage.
Only the exact zero equilibrium can be certified here; nonzero field-shape and
model-scope review are not automated by this inventory.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path

from audit_hv_mesh import audit as mesh_audit
from audit_hv_stabilization import audit as stabilization_audit
from run_hv_sweep import audit_nodes, classify, sha

POINT_METRICS = [
    'theta_max_K', 'joule_W', 'electron_right_A', 'hole_right_A',
    'heat_balance_relative', 'physical_exterior_heat_balance_relative',
    'electron_balance_relative', 'hole_balance_diagnostic',
    'ring_abs_n_relative_mean', 'ring_abs_p_relative_mean',
    'maximum_accepted_residual', 'maximum_accepted_update',
]
ZERO_FIELDS = ['ux_m', 'uy_m', 'uz_m', 'phi_V', 'eta_n', 'eta_p', 'theta_K']
ZERO_METRICS = ['theta_max_K', 'joule_W', 'electron_right_A', 'hole_right_A',
                'ring_abs_n_relative_mean', 'ring_abs_p_relative_mean']
SOURCE = 'ZNO_300K_CANDIDATE_V2; plan=F1_20260906'


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def absolute(path, base):
    # Preserve the ASCII junction path used by the native Fortran executable.
    return Path(os.path.abspath(base / path))


def exact_zero(directory, metrics):
    if metrics.get('voltage_V') != 0 or any(metrics.get(k) != 0 for k in ZERO_METRICS):
        return False
    count = 0
    with (directory / 'nodes.csv').open(encoding='utf-8-sig', newline='') as stream:
        for row in csv.DictReader(stream):
            count += 1
            if any(float(row[k]) != 0 for k in ZERO_FIELDS):
                return False
            # Frozen candidate values from physics/zno_candidate_m.f90; only this
            # material source is eligible for the explicit zero-equilibrium check.
            for key, reference in [('n_m3', 6e22), ('p_m3', (2.9e-4)**2 / 6e22)]:
                if not math.isclose(float(row[key]), reference, rel_tol=1e-14, abs_tol=0.):
                    return False
    return count > 0


def point_record(voltage, directory, manifest_path):
    record = dict(voltage_V=voltage, directory=str(directory),
                  manifest=str(manifest_path), state='MISSING_POINT_EVIDENCE',
                  fully_validated=False, evidence_errors=[], file_hashes={})
    if not (directory / 'execution.json').is_file():
        return record
    try:
        manifest = read_json(manifest_path)
        execution = read_json(directory / 'execution.json')
        for key in ['mesh', 'executable']:
            if sha(Path(manifest[key])) != manifest[key + '_sha256']:
                record['evidence_errors'].append(key.upper() + '_HASH_MISMATCH')
        if execution.get('executable_sha256') != manifest['executable_sha256']:
            record['evidence_errors'].append('EXECUTION_EXECUTABLE_MISMATCH')
        if sha(directory / 'input.nml') != execution.get('input_sha256'):
            record['evidence_errors'].append('INPUT_HASH_MISMATCH')
        if execution.get('timed_out') is not False:
            record['evidence_errors'].append('EXECUTION_TIMEOUT_OR_MISSING_STATUS')
        state, metrics = classify(directory, execution['exit_code'])
        domain = manifest.get('domain_symmetry', 'FULL')
        if metrics.get('voltage_V') != voltage:
            record['evidence_errors'].append('VOLTAGE_MISMATCH')
        if domain not in ('FULL', 'Y_HALF') or metrics.get('domain_symmetry', 'FULL') != domain:
            record['evidence_errors'].append('DOMAIN_MISMATCH')
        if (manifest.get('model') != 'FROZEN_DRIFT_HEAT_ZNO_NT' or metrics.get('model') != 'NT' or
                metrics.get('heat_source') != 'DRIFT_CONDUCTIVITY_NOT_TOTAL_CURRENT_WORK' or
                metrics.get('source') != SOURCE):
            record['evidence_errors'].append('MODEL_MISMATCH')
        if metrics.get('supg_scale') != manifest.get('supg_scale'):
            record['evidence_errors'].append('SUPG_SCALE_MISMATCH')
        if state != 'FAILED_NO_COMPLETE_STATE':
            fields = audit_nodes(directory / 'nodes.csv', Path(manifest['mesh']), voltage,
                                 half_y=domain == 'Y_HALF')
            record['field_audit_recomputed'] = fields
            if not fields['passed']:
                record['evidence_errors'].append('FIELD_AUDIT_FAILED')
        record.update(state=state, domain_symmetry=domain,
                      mesh_sha256=manifest['mesh_sha256'],
                      executable_sha256=manifest['executable_sha256'],
                      supg_scale=manifest['supg_scale'],
                      metrics={k: metrics[k] for k in POINT_METRICS if k in metrics})
        # Keep actual-domain numbers and explicitly named reconstructed integrals.
        factor = 2 if domain == 'Y_HALF' else 1
        record['full_domain_integral_factor'] = factor
        record['reconstructed_full_integrals'] = {
            k: factor * metrics[k] for k in ['joule_W', 'electron_right_A', 'hole_right_A']
            if isinstance(metrics.get(k), (float, int)) and math.isfinite(metrics[k])}
        record['exact_zero_equilibrium_verified'] = (
            not record['evidence_errors'] and state == 'SOLVE_FLUX_PASS_MESH_PENDING' and
            voltage == 0 and exact_zero(directory, metrics))
        if record['exact_zero_equilibrium_verified']:
            record['state'] = 'EXACT_ZERO_EQUILIBRIUM_VERIFIED'
        record['fully_validated'] = record['exact_zero_equilibrium_verified']
        if record['evidence_errors']:
            record['state'] = 'INVALID_POINT_EVIDENCE'
        record['file_hashes'] = {str(path): sha(path) for path in [manifest_path] +
                                [directory / name for name in
                                 ['input.nml', 'execution.json', 'metrics.txt', 'nodes.csv']] if path.is_file()}
    except (OSError, ValueError, TypeError, KeyError, OverflowError) as exc:
        record['state'] = 'INVALID_POINT_EVIDENCE'
        record['fully_validated'] = False
        record['evidence_errors'].append(type(exc).__name__ + ': ' + str(exc))
    # Invalid/nonfinite data remains a failed record, never a nonstandard JSON NaN.
    for key, value in list(record.get('metrics', {}).items()):
        if not isinstance(value, (float, int)) or not math.isfinite(value):
            record['metrics'][key] = None
    return record


def summarize(plan_path):
    plan = read_json(plan_path)
    base = plan_path.parent
    supplied = {}
    for point in plan['points']:
        voltage = point['voltage_V']
        if voltage in supplied or voltage not in [i / 2 for i in range(17)]:
            raise ValueError('Points must be unique members of the 0:0.5:8 voltage grid')
        supplied[voltage] = point
    points = []
    for i in range(17):
        voltage = i / 2
        if voltage not in supplied:
            points.append(dict(voltage_V=voltage, state='NO_REGISTERED_POINT', fully_validated=False))
            continue
        spec = supplied[voltage]
        directory = absolute(spec['directory'], base)
        manifest = absolute(spec.get('manifest', str(directory / 'manifest.json')), base)
        points.append(point_record(voltage, directory, manifest))
    studies = []
    for spec in plan.get('studies', []):
        directories = [absolute(p, base) for p in spec['runs']]
        if spec['kind'] == 'mesh':
            result = mesh_audit(directories, [spec['voltage_V']])
        elif spec['kind'] == 'stabilization':
            result = stabilization_audit(directories)
        else:
            raise ValueError('Unknown study kind')
        studies.append(dict(name=spec['name'], kind=spec['kind'], voltage_V=spec['voltage_V'],
                            scope=spec['scope'], runs=[str(d) for d in directories], result=result))
    supplementary = []
    for spec in plan.get('off_grid_points', []):
        voltage = spec['voltage_V']
        if not math.isfinite(voltage) or not 0 <= voltage <= 8 or voltage in [i / 2 for i in range(17)]:
            raise ValueError('Off-grid entries must be actual non-grid voltages within [0,8]')
        directory = absolute(spec['directory'], base)
        supplementary.append(point_record(voltage, directory, directory / 'manifest.json'))
    return dict(scope='ORIGINAL_MODEL_REGISTERED_EVIDENCE_INVENTORY',
                plan=str(plan_path), plan_sha256=sha(plan_path), points=points, studies=studies,
                off_grid_points=supplementary,
                verified_exact_voltages_V=[p['voltage_V'] for p in points if p['fully_validated']],
                certified_nonzero_interval=None, complete_quantitative_0_8V=False,
                model='FROZEN_DRIFT_HEAT_ZNO_NT',
                limitation='Studies keep their own voltage, mesh, executable and domain scope. '
                           'Nonzero field-shape/model review is not automated; no interpolation, '
                           'interval certificate, dynamical stability or experimental calibration implied.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new output directory')
    result = summarize(Path(os.path.abspath(args.plan)))
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'evidence.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    columns = ['voltage_V', 'state', 'domain_symmetry', 'fully_validated', 'directory',
               'full_domain_integral_factor', 'full_joule_W', 'full_electron_right_A'] + POINT_METRICS
    with (args.output / 'status_0_8V.csv').open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for point in result['points']:
            flat = dict(point, **point.get('metrics', {}))
            flat.update({'full_' + key: value for key, value in point.get('reconstructed_full_integrals', {}).items()})
            writer.writerow({key: flat.get(key, '') for key in columns})
    print(json.dumps(dict(report_written=True, complete_quantitative_0_8V=False,
                          verified_exact_voltages_V=result['verified_exact_voltages_V'])))


if __name__ == '__main__':
    main()
