"""Acceptance summary for the 30-case reduced J.E production (S4 of the 2026-10-22 task package).

Reads only machine evidence (run directories, audits, estimate comparison) and writes
final_review/report.json. Original gates are applied as recorded by the audits; nothing is
reclassified. Derived paper quantities use accepted finest-mesh half-domain states only.
"""
from pathlib import Path
import json
import time
import run_je_production_reduced_20261022 as reduced

base, r = reduced.base, reduced.r
OUT = reduced.OUTPUT/'final_review'


def latest(pattern):
    hits = sorted((reduced.OUTPUT/'audits').glob(pattern))
    if not hits: raise RuntimeError('Missing '+pattern)
    return hits[-1]


def flatten_failures(audit):
    fails = []
    def walk(node, path):
        if isinstance(node, dict):
            if node.get('passed') is False and not any(isinstance(v, dict) and v.get('passed') is False for v in node.values()):
                fails.append(path)
            for k, v in node.items(): walk(v, f'{path}/{k}' if path else k)
    walk(audit, '')
    return sorted(set(fails))


def main():
    jobs = base.jobs()
    states, paper = [], {}
    for job in jobs:
        d = base.OUTPUT/job['name']
        state, m = base.classify_run(d)
        accepted = base.accepted(d)
        entry = dict(name=job['name'], group=job['group'], model=job['model'], radius_um=job['radius'],
                     voltage_V=job['voltage'], domain=job['domain'], counts=job['counts'], state=state, accepted=accepted,
                     seconds=r.read(d/'execution.json')['seconds'] if (d/'execution.json').exists() else None,
                     metrics_sha256=r.sha(d/'metrics.txt') if (d/'metrics.txt').exists() else None,
                     postprocess_recovery=(d/'postprocess_recovery.json').exists())
        if (d/'physical_power.json').exists():
            p = r.read(d/'physical_power.json')
            entry.update(physical_exterior_gap_relative=p['physical_exterior_gap_relative'],
                         independent_vs_native_relative=p['independent_vs_native_relative'])
        states.append(entry)
        if accepted and job['domain'] == 'Y_HALF' and tuple(job['counts']) == base.FINEST and job['voltage'] > 0:
            factor = 2
            paper[(job['model'], job['radius'], job['voltage'])] = dict(
                P_mW=factor*m['joule_W']*1e3, theta_max_K=m['theta_max_K'],
                Rth_K_per_W=m['theta_max_K']/(factor*m['joule_W']), I_n_mA=factor*m['electron_right_A']*1e3,
                ring_theta_mean_K=m['ring_theta_mean_K'], ring_E_mean_V_m=m['ring_E_mean_V_m'],
                ring_Q_mean_W_m3=m['ring_Q_mean_W_m3'], ring_Q_p95_W_m3=m['ring_Q_p95_W_m3'], E_p95_V_m=m['E_p95_V_m'],
                ring_abs_n_relative_mean=m['ring_abs_n_relative_mean'], minimum_mu_n_ratio=m['minimum_mu_n_ratio'],
                E_above_1e6_volume_fraction=m['E_above_1e6_volume_fraction'],
                stress_von_mises_max_Pa=m['stress_von_mises_max_Pa'])
    rel = lambda a, b: a/b-1
    nt_vs_lq = {}
    for v in sorted({k[2] for k in paper if k[0] == 'LQ' and k[1] == 3}):
        lq, nt = paper.get(('LQ', 3, v)), paper.get(('NT', 3, v))
        if lq and nt:
            nt_vs_lq[f'{v:g}'] = {k: rel(nt[k], lq[k]) for k in ('P_mW', 'theta_max_K', 'Rth_K_per_W', 'I_n_mA')}
    radius = {}
    for model in ('LQ', 'NT'):
        for v in (6., 8.):
            pts = {a: paper.get((model, a, v)) for a in (2, 3, 4)}
            if all(pts.values()):
                radius[f'{model}_{v:g}V'] = {
                    f'{k}_2to4': rel(pts[4][k], pts[2][k]) for k in
                    ('P_mW', 'theta_max_K', 'Rth_K_per_W', 'E_p95_V_m', 'ring_Q_mean_W_m3', 'ring_E_mean_V_m', 'ring_theta_mean_K')}
                radius[f'{model}_{v:g}V']['identity_check'] = (1+radius[f'{model}_{v:g}V']['P_mW_2to4'])*(1+radius[f'{model}_{v:g}V']['Rth_K_per_W_2to4'])-1
    lq3 = {v: paper[('LQ', 3, v)] for v in (2., 4., 6., 8.) if ('LQ', 3, v) in paper}
    quadratic = {f'{v:g}': lq3[v]['theta_max_K']/lq3[2.]['theta_max_K']/(v/2)**2-1 for v in lq3}
    # Kohlrausch phi-theta relation: insulated sides/hole, electrodes at T0, sigma(T) and k const give
    # int_0^theta k/sigma dT = phi(V-phi)/2, so theta_max sits on phi=V/2 and does not depend on geometry.
    # LQ (sigma const): theta_max = sigma0 V^2/(8k); NT with sigma ~ (T/T0)^-beta, beta=2:
    # theta_max = T0[(1 + 3 sigma0 V^2/(8 k T0))^(1/3) - 1]. Material values from zno_candidate (Table 1).
    q, n0, mu0, k, T0 = 1.602176634e-19, 6e22, 0.0205, 60.0, 300.0
    s0 = q*n0*mu0
    kohlrausch = {}
    for (model, a, v), st in sorted(paper.items()):
        ideal = s0*v*v/(8*k) if model == 'LQ' else T0*((1+3*s0*v*v/(8*k*T0))**(1/3)-1)
        kohlrausch[f'{model}_r{a}_V{v:g}'] = dict(theta_max_K=st['theta_max_K'], ideal_K=ideal,
                                                  relative=st['theta_max_K']/ideal-1)
    audit_path, compare_path = latest('audits_*.json'), latest('estimate_comparison_*.json')
    audit, compare = r.read(audit_path), r.read(compare_path)
    OUT.mkdir(exist_ok=True)
    report = dict(
        created=time.strftime('%Y-%m-%dT%H:%M:%S%z'), scope='REDUCED_30_CASE_TRACE_JE_PRODUCTION',
        release=str(base.RELEASE), release_sha256=r.sha(base.RELEASE),
        case_matrix_sha256=r.sha(reduced.INPUTS/'case_matrix.csv'),
        accepted=sum(s['accepted'] for s in states), total=len(states), states=states,
        audit_file=str(audit_path), audit_sha256=r.sha(audit_path), audit_failures=flatten_failures(audit),
        estimate_comparison_file=str(compare_path), estimate_comparison_sha256=r.sha(compare_path),
        estimate_flagged=compare['flagged'],
        paper_states={f'{k[0]}_r{k[1]}_V{k[2]:g}': v for k, v in sorted(paper.items())},
        derived=dict(nt_vs_lq_r3=nt_vs_lq, radius_2_to_4um=radius, lq_theta_quadratic_deviation=quadratic,
                     kohlrausch=dict(sigma0_S_per_m=s0, k_W_per_mK=k, T0_K=T0, beta=2,
                                     note='Ideal relation omits diffusion current, space charge, piezo/pyro coupling and D(T)',
                                     states=kohlrausch,
                                     max_abs_relative_LQ=max(abs(x['relative']) for kk, x in kohlrausch.items() if kk.startswith('LQ')),
                                     max_abs_relative_NT=max(abs(x['relative']) for kk, x in kohlrausch.items() if kk.startswith('NT')))),
        thresholds_changed=False, physics_changed=False, fail_records_retained=True,
        limitations=['Observable classification from r=3 6/8 V endpoints; intermediate voltages inherit it',
                     'r=2/4 SUPG sensitivity not measured; inherited from r=3',
                     'Hole quantities excluded from paper claims; LQ ring carrier perturbations are numerical zeros',
                     'J.E current-work heat source (R=0) is not a thermodynamic closure; no r=0 or N0 cases',
                     'NT V005 post-processing completed after an encoding crash without re-solving (postprocess_recovery.json)'])
    r.save(OUT/'report.json', report)
    print(json.dumps(dict(accepted=report['accepted'], total=report['total'], audit_failures=report['audit_failures'],
                          estimate_flagged=report['estimate_flagged'], derived=report['derived']), indent=1, ensure_ascii=False))


if __name__ == '__main__':
    main()
