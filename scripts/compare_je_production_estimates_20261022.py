"""Compare solved production states with the estimates pre-registered on 2026-10-21.

Review aid only: a deviation above the registered 2% flags the case for review; it
never accepts, rejects or replaces a solved state. Only states the runner accepted are
compared; zero-voltage states are listed with their absolute values. Y-half power is
doubled for the full plate exactly as in the production table; theta_max is not scaled.
"""
from pathlib import Path
import argparse
import csv
import time
import run_je_production_20261021 as base

r = base.r
ESTIMATES = r.STUDY/'paper_scope/production_release_scoped_20261021/preregistered_estimates.csv'
RELEASE_20261021 = ESTIMATES.parent/'release.json'


def estimates():
    with ESTIMATES.open(encoding='utf-8-sig') as f:
        return {row['case']: row for row in csv.DictReader(f)}


def estimate_key(job):
    # Estimates exist for the finest half-domain production states only.
    if job['domain'] != 'Y_HALF' or tuple(job['counts']) != base.FINEST: return None
    return f"R{job['radius']:02d}_{job['model']}_V{round(job['voltage']*10):03d}"


def relative(value, reference):
    return (value-reference)/reference if reference else None


def compare_state(job, metrics, estimate, limit):
    factor = 2 if job['domain'] == 'Y_HALF' else 1
    p_mw = factor*metrics['joule_W']*1e3
    theta = metrics['theta_max_K']
    row = dict(name=job['name'], case=estimate['case'], model=job['model'], radius_um=job['radius'], voltage_V=job['voltage'],
               joule_mW=p_mw, joule_mW_estimate=float(estimate['joule_mW_estimate']),
               theta_max_K=theta, theta_max_K_estimate=float(estimate['theta_max_K_estimate']))
    if job['voltage'] == 0:
        row.update(review='ZERO_VOLTAGE_ABSOLUTE_ONLY')
        return row
    rth = theta/(p_mw*1e-3)
    row.update(Rth_K_per_W=rth, Rth_K_per_W_estimate=float(estimate['Rth_K_per_W_estimate']),
               joule_deviation=relative(p_mw, row['joule_mW_estimate']),
               theta_deviation=relative(theta, row['theta_max_K_estimate']),
               Rth_deviation=relative(rth, float(estimate['Rth_K_per_W_estimate'])))
    worst = max(abs(row[k]) for k in ('joule_deviation', 'theta_deviation', 'Rth_deviation'))
    row.update(max_abs_deviation=worst, review='REVIEW_DEVIATION_ABOVE_REGISTERED_LIMIT' if worst > limit else 'WITHIN_REGISTERED_LIMIT')
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test', action='store_true', help='compare the certified trace 6/8 V states and write nothing')
    a = p.parse_args()
    limit = r.read(RELEASE_20261021)['estimate_model']['deviation_review_threshold_relative']
    table = estimates()
    rows, missing = [], []
    if a.self_test:
        trace = r.STUDY/'stage5_controlled_numerics_20261021/trace_joint_supg'
        for v in (6, 8):
            job = dict(name=f'A48_R22_Z18_S1_V0{v}0', model='NT', radius=3, voltage=float(v), domain='Y_HALF', counts=list(base.FINEST))
            rows.append(compare_state(job, base.sweep.read_metrics(trace/job['name']/'metrics.txt'), table[estimate_key(job)], limit))
        for row in rows: print({k: row[k] for k in ('case', 'joule_mW', 'joule_mW_estimate', 'theta_max_K', 'max_abs_deviation', 'review')})
        return 0
    for job in base.jobs():
        key = estimate_key(job)
        if key is None: continue
        directory = base.OUTPUT/job['name']
        state, metrics = base.classify_run(directory)
        if not base.accepted(directory):
            missing.append(dict(name=job['name'], state=state)); continue
        rows.append(compare_state(job, metrics, table[key], limit))
    out = base.OUTPUT/'audits'; out.mkdir(exist_ok=True)
    stamp = time.strftime('%Y%m%dT%H%M%S')
    if rows:
        fields = list(dict.fromkeys(k for row in rows for k in row))
        with (out/f'estimate_comparison_{stamp}.csv').open('x', encoding='utf-8-sig', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    r.save(out/f'estimate_comparison_{stamp}.json', dict(
        estimates=str(ESTIMATES), estimates_sha256=r.sha(ESTIMATES), review_limit_relative=limit,
        compared=len(rows), flagged=[x['name'] for x in rows if x['review'].startswith('REVIEW')],
        not_compared_not_accepted=missing, role='Review aid only; never a pass/fail gate or a substitute for solved states',
        script_sha256=r.sha(Path(__file__).resolve())))
    print(dict(compared=len(rows), flagged=sum(x['review'].startswith('REVIEW') for x in rows), not_accepted=len(missing)))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
