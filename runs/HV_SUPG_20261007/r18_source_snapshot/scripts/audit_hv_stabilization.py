"""Audit SUPG scales 0.5/1/2 at one voltage on the same mesh and executable.

The predeclared tolerance is 1% across all three runs for the listed global and
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


def audit(runs):
    records, errors = [], []
    if len(runs) != 3:
        errors.append('THREE_RUNS_REQUIRED')
    for directory in runs:
        try:
            manifest = json.loads((directory / 'manifest.json').read_text(encoding='utf-8-sig'))
            execution = json.loads((directory / 'execution.json').read_text(encoding='utf-8-sig'))
            state, m = classify(directory, execution['exit_code'])
            if state != 'SOLVE_FLUX_PASS_MESH_PENDING':
                errors.append(f'{directory.name}:SOLVE_OR_FLUX_FAIL')
            fields = json.loads((directory / 'field_audit.json').read_text(encoding='utf-8-sig'))
            if fields.get('passed') is not True:
                errors.append(f'{directory.name}:FIELD_AUDIT_FAILED')
            for path_key, hash_key in [('mesh', 'mesh_sha256'), ('executable', 'executable_sha256')]:
                if sha(Path(manifest[path_key])) != manifest[hash_key]:
                    errors.append(f'{directory.name}:{hash_key.upper()}_MISMATCH')
            if manifest.get('supg_scale') != m.get('supg_scale'):
                errors.append(f'{directory.name}:SUPG_SCALE_MISMATCH')
            if not all(isinstance(m.get(k), (float, int)) and math.isfinite(m[k]) for k in METRICS):
                errors.append(f'{directory.name}:MISSING_OR_NONFINITE_METRIC')
            records.append(dict(directory=str(directory.absolute()), manifest=manifest, metrics=m))
        except (OSError, ValueError, TypeError, KeyError) as exc:
            errors.append(f'{directory.name}:INCOMPLETE_EVIDENCE:{type(exc).__name__}')
    if len(records) == 3:
        scales = [r['manifest'].get('supg_scale') for r in records]
        if (not all(isinstance(v, (int, float)) and math.isfinite(v) for v in scales) or
            sorted(scales) != [.5, 1., 2.]):
            errors.append('SCALES_MUST_BE_0.5_1_2')
        for key in ['mesh_sha256', 'executable_sha256', 'model']:
            values = [r['manifest'].get(key) for r in records]
            if None in values or len(set(values)) != 1:
                errors.append('INCONSISTENT_' + key.upper())
        for key in ['voltage_V', 'model', 'heat_source', 'source', 'supg_directional_tau']:
            values = [r['metrics'].get(key) for r in records]
            if None in values or len(set(values)) != 1:
                errors.append('INCONSISTENT_' + key.upper())
    comparisons = {}
    if len(records) == 3 and not errors:
        records.sort(key=lambda r: r['manifest']['supg_scale'])
        for key in METRICS:
            values = [r['metrics'][key] for r in records]
            span = (max(values) - min(values)) / max(max(abs(v) for v in values), 1e-300)
            comparisons[key] = dict(values=values, relative_span=span, limit=TOLERANCE,
                                    passed=span < TOLERANCE)
    passed = bool(comparisons) and not errors and all(c['passed'] for c in comparisons.values())
    return dict(scope='SINGLE_VOLTAGE_SAME_MESH_STABILIZATION_SENSITIVITY',
                scales=[.5, 1., 2.], runs=[str(p.absolute()) for p in runs],
                errors=sorted(set(errors)), comparisons=comparisons,
                covered_metrics_passed=passed, fully_validated=False,
                limitation='Still requires mesh convergence, field-shape review and model-scope validation; no full-interval certificate.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs', nargs=3, type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); result = audit(a.runs)
    with a.output.open('x', encoding='utf-8') as f:
        json.dump(result, f, indent=2, allow_nan=False)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result['covered_metrics_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
