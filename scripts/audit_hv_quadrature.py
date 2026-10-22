"""Audit log-weighted integration orders 3/5 at one voltage on the same mesh and executable.

The predeclared tolerance is 1% across both runs for the listed global and
annular metrics. Passing this comparison does not certify the mesh, the full
voltage interval, dynamical stability, or experimental material validity.
"""
import argparse, json, math
from pathlib import Path
from run_hv_sweep import classify, read_metrics, sha

METRICS = ['joule_W', 'theta_max_K', 'electron_right_A', 'ring_theta_mean_K',
           'ring_E_mean_V_m', 'ring_Q_mean_W_m3', 'ring_abs_n_relative_mean',
           'ring_abs_p_relative_mean', 'ring_Q_p95_W_m3', 'E_p95_V_m']
TOLERANCE = .01


def audit(runs, include_local_fields=False):
    records, errors = [], []
    if len(runs) != 2:
        errors.append('TWO_RUNS_REQUIRED')
    for directory in runs:
        try:
            manifest = json.loads((directory / 'manifest.json').read_text(encoding='utf-8-sig'))
            execution = json.loads((directory / 'execution.json').read_text(encoding='utf-8-sig'))
            if sha(directory/'input.nml')!=execution.get('input_sha256'):
                errors.append(f'{directory.name}:INPUT_HASH_MISMATCH')
            if execution.get('executable_sha256')!=manifest.get('executable_sha256'):
                errors.append(f'{directory.name}:EXECUTION_HASH_MISMATCH')
            state, m = classify(directory, execution['exit_code'])
            if m.get('transport_discretization')!='HEX20_LOG_WEIGHTED_EXPERIMENTAL':
                errors.append(f'{directory.name}:QUADRATURE_STUDY_REQUIRES_LOG_WEIGHTED')
            if state != 'SOLVE_FLUX_PASS_MESH_PENDING':
                errors.append(f'{directory.name}:SOLVE_OR_FLUX_FAIL')
            fields = json.loads((directory / 'field_audit.json').read_text(encoding='utf-8-sig'))
            if fields.get('passed') is not True:
                errors.append(f'{directory.name}:FIELD_AUDIT_FAILED')
            for path_key, hash_key in [('mesh', 'mesh_sha256'), ('executable', 'executable_sha256')]:
                if sha(Path(manifest[path_key])) != manifest[hash_key]:
                    errors.append(f'{directory.name}:{hash_key.upper()}_MISMATCH')
            if manifest.get('log_weighted_quadrature_order') != m.get('log_weighted_quadrature_order'):
                errors.append(f'{directory.name}:QUADRATURE_ORDER_MISMATCH')
            if manifest.get('domain_symmetry','FULL') != m.get('domain_symmetry','FULL'):
                errors.append(f'{directory.name}:DOMAIN_SYMMETRY_MISMATCH')
            if not all(isinstance(m.get(k), (float, int)) and math.isfinite(m[k]) for k in METRICS):
                errors.append(f'{directory.name}:MISSING_OR_NONFINITE_METRIC')
            records.append(dict(directory=str(directory.absolute()), manifest=manifest, metrics=m))
        except (OSError, ValueError, TypeError, KeyError) as exc:
            errors.append(f'{directory.name}:INCOMPLETE_EVIDENCE:{type(exc).__name__}')
    if len(records) == 2:
        if len({r['manifest'].get('domain_symmetry','FULL') for r in records})!=1:
            errors.append('INCONSISTENT_DOMAIN_SYMMETRY')
        scales = [r['manifest'].get('log_weighted_quadrature_order') for r in records]
        if (not all(isinstance(v, (int, float)) and math.isfinite(v) for v in scales) or
            sorted(scales) != [3,5]):
            errors.append('ORDERS_MUST_BE_3_5')
        for key in ['mesh_sha256', 'executable_sha256', 'model']:
            values = [r['manifest'].get(key) for r in records]
            if None in values or len(set(values)) != 1:
                errors.append('INCONSISTENT_' + key.upper())
        for key in ['voltage_V', 'model', 'heat_source', 'source', 'supg_directional_tau']:
            values = [r['metrics'].get(key) for r in records]
            if None in values or len(set(values)) != 1:
                errors.append('INCONSISTENT_' + key.upper())
    comparisons = {}
    if len(records) == 2 and not errors:
        records.sort(key=lambda r: r['manifest']['log_weighted_quadrature_order'])
        for key in METRICS:
            values = [r['metrics'][key] for r in records]
            span = (max(values) - min(values)) / max(max(abs(v) for v in values), 1e-300)
            comparisons[key] = dict(values=values, relative_span=span, limit=TOLERANCE,
                                    passed=span < TOLERANCE)
    passed = bool(comparisons) and not errors and all(c['passed'] for c in comparisons.values())
    local_fields={};local_error=None
    if include_local_fields and passed:
        try:
            import numpy as np
            from audit_hv_field_shape import difference
            a,b=[np.atleast_1d(np.genfromtxt(p/'nodes.csv',delimiter=',',names=True)) for p in runs]
            if len(a)!=len(b) or any(not np.array_equal(a[k],b[k]) for k in ['node','x_m','y_m','z_m']):
                raise ValueError('Different node coordinates or order')
            for key in ['phi_V','theta_K']:
                local_fields[key]=difference(a[key],b[key])
            scale=max(max(float(np.max(abs(v[k]))) for k in ['ux_m','uy_m','uz_m']) for v in [a,b])
            for key in ['ux_m','uy_m','uz_m']:
                local_fields[key]=difference(a[key],b[key],max(scale,1e-300))
            for key in ['eta_n','eta_p']:
                local_fields[key+'_density_perturbation']=difference(np.expm1(a[key]),np.expm1(b[key]))
            for item in local_fields.values():
                item.update(limit=TOLERANCE,passed=item['maximum_relative']<TOLERANCE)
        except (OSError,ValueError,TypeError,KeyError,OverflowError) as exc:
            local_error=type(exc).__name__+': '+str(exc)
    all_requested=passed and (not include_local_fields or
        (bool(local_fields) and local_error is None and all(c['passed'] for c in local_fields.values())))
    return dict(scope='SINGLE_VOLTAGE_SAME_MESH_QUADRATURE_SENSITIVITY',
                orders=[3,5], runs=[str(p.absolute()) for p in runs],
                errors=sorted(set(errors)), comparisons=comparisons,
                local_field_checks_requested=include_local_fields,local_field_comparisons=local_fields,
                local_field_error=local_error,all_requested_checks_passed=all_requested,
                node_hashes={str(p.absolute()):sha(p/'nodes.csv') for p in runs if (p/'nodes.csv').is_file()},
                covered_metrics_passed=passed, fully_validated=False,
                limitation='Still requires mesh convergence, field-shape review and model-scope validation; no full-interval certificate.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs', nargs=2, type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--include-local-fields',action='store_true',help='Additionally require 1%% maximum nodal-field agreement')
    a = p.parse_args(); result = audit(a.runs,a.include_local_fields)
    with a.output.open('x', encoding='utf-8') as f:
        json.dump(result, f, indent=2, allow_nan=False)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result['all_requested_checks_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
