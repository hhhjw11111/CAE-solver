"""Observable-scoped production release for the trace-metric SUPG variant; never launches a solver.

Supersedes production_release_scoped_20261021 only in the NT numerical variant: the
classification is recomputed from the stage5 controlled trace review, whose overall
FAIL and production_allowed_by_endpoint_prerequisite=false stay on record unchanged.
The 2026-10-21 release, its pre-registered estimates and all thresholds are kept.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import run_je_reconstruction_mesh_20261018 as r

PAPER = r.STUDY/'paper_scope'
TRACE_REVIEW = r.STUDY/'stage5_controlled_numerics_20261021/final_review/report.json'
PREVIOUS = PAPER/'production_release_scoped_20261021'
OUTPUT = PAPER/'production_release_scoped_20261022'
CERTIFIED_INPUT = r.STUDY/'stage5_controlled_numerics_20261021/trace_joint_supg/A48_R22_Z18_S1_V060/input.nml'
CERTIFIED_MANIFEST = CERTIFIED_INPUT.parent/'manifest.json'
AUTHORIZATION = ('Go with plan A: first write the scope-restricted release decision (it does not change the '
                 'existing FAIL), then rebuild the production inputs with the certified settings. Show me the '
                 'case matrix before starting any runs.')
EXCLUDED = ['ring_abs_p_relative_mean']


def classify(report):
    mesh = {v['voltage_V']: v for entry in report['mesh'] for v in entry['voltages']}
    supg = {}
    for entry in report['supg']:
        voltage = 6 if entry['runs'][0].endswith('V060') else 8 if entry['runs'][0].endswith('V080') else None
        if voltage is None or entry['scales'] != [0.5, 1.0, 2.0]: raise RuntimeError('Unexpected SUPG entry')
        supg[voltage] = entry['comparisons']
    if sorted(mesh) != [6, 8] or sorted(supg) != [6, 8]: raise RuntimeError('Trace review lacks 6/8 V')
    rows = []
    for metric in supg[6]:
        m6, m8 = mesh[6]['changes'][metric], mesh[8]['changes'][metric]
        s6, s8 = supg[6][metric], supg[8][metric]
        mesh_pass = max(m6['L1_L2'], m8['L1_L2']) < m6['limit']
        supg_pass = s6['passed'] and s8['passed']
        if metric in EXCLUDED: use = 'EXCLUDED_FROM_PAPER_CLAIMS_FAIL_RETAINED'
        elif mesh_pass and supg_pass: use = 'PRIMARY_RESULT'
        else: raise RuntimeError('Unexpected failing paper metric: '+metric)
        if metric in EXCLUDED and mesh_pass and supg_pass: raise RuntimeError('Excluded metric unexpectedly passed')
        rows.append(dict(metric=metric, mesh_L1_L2_6V=m6['L1_L2'], mesh_L1_L2_8V=m8['L1_L2'], mesh_limit=m6['limit'],
                         mesh_pass=mesh_pass, supg_span_6V=s6['relative_span'], supg_span_8V=s8['relative_span'],
                         supg_limit=s6['limit'], supg_pass=supg_pass, paper_use=use,
                         uncertainty_band_relative=max(m6['L1_L2'], m8['L1_L2'], s6['relative_span'], s8['relative_span'])))
    local_ok = all(entry['passed'] for entry in report['local'])
    rows.append(dict(metric='local_hole_field_and_hole_profiles', mesh_pass=local_ok, supg_pass=False,
                     paper_use='EXCLUDED_FROM_PAPER_CLAIMS_FAIL_RETAINED'))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    report = r.read(TRACE_REVIEW)
    if (report.get('production_allowed_by_endpoint_prerequisite') is not False or report.get('production_started') is not False
            or report.get('original_thresholds_unchanged') is not True):
        raise RuntimeError('Trace review state differs from the reviewed FAIL record')
    previous = r.read(PREVIOUS/'release.json')
    if previous.get('release_type') != 'OBSERVABLE_SCOPED_PRODUCTION': raise RuntimeError('Previous release missing')
    certified = r.read(CERTIFIED_MANIFEST)
    if certified['supg_directional_tau'] is not False or certified['runtime_env']['TPS_PARDISO_SPLIT_HOLES'] != '1':
        raise RuntimeError('Certified run is not the trace variant')
    rows = classify(report)
    OUTPUT.mkdir(exist_ok=True)
    with (OUTPUT/'observable_scope.csv').open('x', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    primary = [x['metric'] for x in rows if x['paper_use'] == 'PRIMARY_RESULT']
    excluded = [x['metric'] for x in rows if x['paper_use'].startswith('EXCLUDED')]
    r.save(OUTPUT/'release.json', dict(
        released_at=datetime.now(timezone.utc).isoformat(),
        authorization=dict(source='user chat 2026-10-22', text=AUTHORIZATION),
        release_type='OBSERVABLE_SCOPED_PRODUCTION', nt_variant='trace',
        thresholds_changed=False, physics_changed=False, original_fail_records_retained=True,
        trace_review_overall='FAIL', trace_review_production_allowed_field_unchanged=False,
        supersedes=dict(path=str(PREVIOUS/'release.json'), sha256=r.sha(PREVIOUS/'release.json'),
                        scope='NT numerical variant and observable classification only; record retained unchanged'),
        primary=primary, with_band=[], excluded=excluded,
        production_matrix=previous['production_matrix'],
        nt_numerics=dict(transport_supg=True, supg_scale=1.0, supg_directional_tau=False, TPS_PARDISO_SPLIT_HOLES='1',
                         domain='Y_HALF', r3_finest_half_mesh_sha256=certified['mesh_sha256'],
                         executable_sha256=certified['executable_sha256'], log_weighted_quadrature_order=3),
        lq_numerics='Unchanged: transport_supg=.false.; the nonlinear SUPG gate does not apply',
        mesh='Y-half A48/R22/Z18, byte-identical to the certified trace runs',
        preregistered_estimates=dict(path=str(PREVIOUS/'preregistered_estimates.csv'),
                                     sha256=r.sha(PREVIOUS/'preregistered_estimates.csv'),
                                     note='Registered 2026-10-21 before any production solve; kept, not recalibrated'),
        run_gate='No native solve until the user approves case_matrix.csv; the runner requires its SHA-256',
        limitations=[
            'Classification uses r=3 6/8 V trace evidence only; intermediate voltages inherit it',
            'r=2/4 inherit the r=3 SUPG classification; each gets its own 8 V A40/A48 mesh pair and NT local-field audit',
            'r=0 (no hole) is outside this matrix: the radius mesh generator supports 2/3/4 um only',
            'LQ half-domain equivalence is checked in the lq_check group before LQ results are used',
            'J.E heat source is a diagnostic current-work source, not a thermodynamic closure',
            'Hole quantities stay reported where produced but carry no paper claim'],
        evidence_hashes={str(p): r.sha(p) for p in [TRACE_REVIEW, CERTIFIED_INPUT, CERTIFIED_MANIFEST, Path(__file__).resolve(),
                                                     OUTPUT/'observable_scope.csv']}))
    print(dict(primary=primary, excluded=excluded, output=str(OUTPUT)))


if __name__ == '__main__':
    main()
