"""User-approved observable-scoped production release and pre-registered estimates.

Thresholds, physics and the 2026-10-21 FAIL records are unchanged. Each audited
observable is classified from the final joint mesh/SUPG report; failing hole
metrics are excluded from paper claims, not converted to PASS. The estimates are
registered before any production solve and are never used as results.
"""
from pathlib import Path
import csv
import hashlib
import json
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
SCOPE = ROOT/'runs/JE_RECONSTRUCTION_20261018/paper_scope'
REPORT = SCOPE/'final_joint_supg_review_20261021/report.json'
OUTPUT = SCOPE/'production_release_scoped_20261021'
LQ3_SOURCE = ROOT.parent/'output/paper_draft_20261008/data/radius_source.csv'
AUTHORIZATION = '同意按分观测量限定范围放行，从阶段1开始执行，并补上 LQ r=2/4 的对照。'
T0 = 300.0


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)


def classify(report):
    mesh = {r['voltages'][0]['voltage_V']: r['voltages'][0]['changes'] for r in report['mesh']}
    supg = {}
    for item in report['supg']:
        voltage = int(round(float(Path(item['runs'][0]).name.split('_V')[-1])/10))
        supg[voltage] = item['comparisons']
    rows = []
    for metric in mesh[6]:
        mesh_values = {v: mesh[v][metric] for v in (6, 8)}
        mesh_pass = all(x['L1_L2'] <= x['limit'] for x in mesh_values.values())
        supg_values = {v: supg[v].get(metric) for v in (6, 8)}
        supg_pass = all(x is None or x['passed'] for x in supg_values.values())
        if mesh_pass and supg_pass: use = 'PRIMARY_RESULT'
        elif mesh_pass: use = 'REPORT_WITH_SUPG_SENSITIVITY_BAND'
        else: use = 'EXCLUDED_FROM_PAPER_CLAIMS_FAIL_RETAINED'
        rows.append(dict(metric=metric,
            mesh_L1_L2_6V=mesh_values[6]['L1_L2'], mesh_L1_L2_8V=mesh_values[8]['L1_L2'],
            mesh_limit=mesh_values[6]['limit'], mesh_pass=mesh_pass,
            supg_span_6V=supg_values[6]['relative_span'] if supg_values[6] else None,
            supg_span_8V=supg_values[8]['relative_span'] if supg_values[8] else None,
            supg_limit=supg_values[6]['limit'] if supg_values[6] else None, supg_pass=supg_pass,
            paper_use=use,
            uncertainty_band_relative=max(mesh_values[6]['L1_L2'], mesh_values[8]['L1_L2'],
                *(x['relative_span'] for x in supg_values.values() if x))))
    rows.append(dict(metric='local_hole_field_and_hole_profiles', mesh_L1_L2_6V=None, mesh_L1_L2_8V=None,
        mesh_limit=None, mesh_pass=False, supg_span_6V=None, supg_span_8V=None, supg_limit=None,
        supg_pass=False, paper_use='EXCLUDED_FROM_PAPER_CLAIMS_FAIL_RETAINED', uncertainty_band_relative=None))
    return rows


def lq_reference():
    # 3 V LQ L0 states; LQ is linear in voltage, so P and theta scale with V^2.
    table = {}
    with LQ3_SOURCE.open(encoding='utf-8-sig') as stream:
        for row in csv.DictReader(stream):
            if row['model'] == 'LQ':
                table[int(row['radius_um'])] = (float(row['joule_W']), float(row['theta_max_K']))
    return table


def estimates(report):
    lq3 = lq_reference()
    power = {p['voltage_V']: p['full_reconstructed_JE_mW']*1e-3 for p in report['power'] if p['supg_scale'] == 1.0}
    theta = {}
    for item in report['supg']:
        v = int(round(float(Path(item['runs'][0]).name.split('_V')[-1])/10))
        theta[v] = item['comparisons']['theta_max_K']['values'][1]
    lq = lambda a, v: (lq3[a][0]*v*v/9, lq3[a][1]*v*v/9)
    # Two-parameter mobility/conductivity feedback calibrated on r=3 NT 6/8 V.
    calib = []
    for v in (6, 8):
        pl, tl = lq(3, v)
        calib.append((((pl/power[v])**0.5-1)*T0/theta[v], ((theta[v]/power[v])/(tl/pl)-1)/theta[v]))
    s = sum(x[0] for x in calib)/2; c = sum(x[1] for x in calib)/2

    def nt(a, v):
        pl, tl = lq(a, v)
        if v == 0: return 0., 0.
        rth, t = tl/pl, tl
        for _ in range(200):
            p = pl*(1+s*t/T0)**-2
            t = rth*(1+c*t)*p
        return p, t
    rows = []
    for i in range(17):
        v = 0.5*i
        for model, a in (('LQ', 3), ('NT', 3)):
            p, t = lq(a, v) if model == 'LQ' else nt(a, v)
            rows.append(dict(case=f'R{a:02d}_{model}_V{int(v*10):03d}', radius_um=a, model=model, voltage_V=v,
                             joule_mW_estimate=p*1e3, theta_max_K_estimate=t,
                             Rth_K_per_W_estimate=t/p if p else None))
    for a in (2, 4):
        for v in (6, 8):
            for model in ('NT', 'LQ'):
                p, t = lq(a, v) if model == 'LQ' else nt(a, v)
                rows.append(dict(case=f'R{a:02d}_{model}_V{v*10:03d}', radius_um=a, model=model, voltage_V=v,
                                 joule_mW_estimate=p*1e3, theta_max_K_estimate=t, Rth_K_per_W_estimate=t/p))
    model = dict(LQ='P,theta from 3 V LQ L0 states scaled by (V/3)^2 (linear transport)',
                 NT='P=P_LQ*(1+s*theta_max/T0)^-2; Rth=Rth_LQ*(1+c*theta_max); self-consistent fixed point',
                 s=s, c_per_K=c, T0_K=T0, calibration='r=3 finest half NT 6/8 V, SUPG=1',
                 status='ESTIMATE_ONLY_NOT_A_RESULT', deviation_review_threshold_relative=0.02)
    return rows, model


def main():
    report = json.loads(REPORT.read_text(encoding='utf-8-sig'))
    if report['production_release'] is not False or report['original_thresholds_unchanged'] is not True:
        raise RuntimeError('Unexpected final review state')
    OUTPUT.mkdir(exist_ok=False)
    rows = classify(report)
    with (OUTPUT/'observable_scope.csv').open('x', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    predicted, model = estimates(report)
    with (OUTPUT/'preregistered_estimates.csv').open('x', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(predicted[0]))
        writer.writeheader(); writer.writerows(predicted)
    save(OUTPUT/'release.json', dict(
        released_at=datetime.now(timezone.utc).isoformat(),
        authorization=dict(source='user chat 2026-10-21', text=AUTHORIZATION),
        release_type='OBSERVABLE_SCOPED_PRODUCTION', thresholds_changed=False, physics_changed=False,
        original_fail_records_retained=True, original_production_release_field_unchanged=False,
        primary=[r['metric'] for r in rows if r['paper_use'] == 'PRIMARY_RESULT'],
        with_band=[r['metric'] for r in rows if r['paper_use'] == 'REPORT_WITH_SUPG_SENSITIVITY_BAND'],
        excluded=[r['metric'] for r in rows if r['paper_use'].startswith('EXCLUDED')],
        production_matrix=dict(r3_LQ_V=[0.5*i for i in range(17)], r3_NT_V=[0.5*i for i in range(17)],
                               radius_NT=dict(radii_um=[2, 4], voltages_V=[6, 8]),
                               radius_LQ=dict(radii_um=[2, 4], voltages_V=[6, 8], added_by_user=True)),
        mesh='Y-half A48/R22/Z18 (full/half equivalence passed six times on A32/A40 at 6/8 V)',
        limitations=['Scope classification uses r=3 6/8 V evidence; other radii get their own 8 V two-level mesh check',
                     'LQ half-domain equivalence must be checked separately before LQ production',
                     'J.E heat source is a diagnostic current-work source, not a thermodynamic closure',
                     'Estimates are pre-registered for review only and never substitute for solved states'],
        estimate_model=model,
        evidence_hashes={str(REPORT): sha(REPORT), str(LQ3_SOURCE): sha(LQ3_SOURCE), str(Path(__file__)): sha(__file__),
                         'observable_scope.csv': sha(OUTPUT/'observable_scope.csv'),
                         'preregistered_estimates.csv': sha(OUTPUT/'preregistered_estimates.csv')}))
    print(json.dumps([(r['metric'], r['paper_use']) for r in rows], ensure_ascii=False, indent=1))
    print('s=%.4f c=%.3e' % (model['s'], model['c_per_K']))


if __name__ == '__main__':
    main()
