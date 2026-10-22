"""Trace-metric JE production: wraps run_je_production_20261021 without editing it.

Released by paper_scope/production_release_scoped_20261022/release.json (nt_variant=trace).
`prepare` writes every production input for review and launches nothing. `execute` needs
the SHA-256 of the reviewed case_matrix.csv, and checks each input the base runner builds
against the reviewed preview (seeded runs differ only in the seed's short path).
"""
from pathlib import Path
import argparse
import csv
import re
import time
import run_je_production_20261021 as base

r = base.r
RELEASE = r.STUDY/'paper_scope/production_release_scoped_20261022/release.json'
INPUTS = r.STUDY/'paper_scope/production_inputs_trace_20261022'
CERTIFIED = r.STUDY/'stage5_controlled_numerics_20261021/trace_joint_supg/A48_R22_Z18_S1_V060'
SEED_PLACEHOLDER = 'SEED_RESOLVED_AT_RUN_TIME'
# Keys a production input may differ in from the certified trace input; every other key must match.
PER_CASE_KEYS = {'case_id', 'mesh_file', 'target_voltage', 'hole_radius', 'restart_nodes_file', 'restart_voltage_factor',
                 'restart_joule_factor', 'continuation_initial_step', 'continuation_maximum_step'}
LQ_KEYS = {'material_profile', 'transport_supg', 'supg_directional_tau'}
PURPOSE = dict(lq_check='LQ full/half equivalence and r=3 8 V three-level mesh', lq_sweep='LQ 0-8 V curve, r=3',
               nt_sweep='NT 0-8 V curve, r=3 (continuation chain)', radius='Radius comparison, 6/8 V',
               radius_mesh='Radius-specific 8 V A40/A48 mesh pair')


def release():
    data = r.read(RELEASE)
    if (data.get('release_type') != 'OBSERVABLE_SCOPED_PRODUCTION' or data.get('nt_variant') != 'trace'
            or data.get('thresholds_changed') is not False or data.get('physics_changed') is not False):
        raise RuntimeError('Trace release record missing or altered')
    return data


def namelist(text):
    values = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line[0] in '&/': continue
        for item in re.split(r',(?=\s*[A-Za-z_]\w*\s*=)', line.rstrip(',')):
            key, value = item.split('=', 1)
            values[key.strip()] = value.strip()
    return values


def same(a, b):
    try:
        return all(float(x) == float(y) for x, y in zip(a.split(','), b.split(','))) and a.count(',') == b.count(',')
    except ValueError:
        return a == b


def preview(job):
    text, changed = base.make_input(job, None)
    if job['seed']:
        text = r.replace(text, 'restart_nodes_file', "'"+SEED_PLACEHOLDER+"'")
        text = r.replace(text, 'restart_joule_factor', '1.0')
    return text, changed


def seed_location(job):
    if not job['seed']: return ''
    kind, source = job['seed']
    return str((base.OUTPUT/source/'nodes.csv') if kind == 'job' else (base.OUTPUT/'seeds'/job['name']/'nodes.csv'))


def check_against_certified(job, text, certified):
    values = namelist(text)
    allowed = PER_CASE_KEYS | (LQ_KEYS if job['model'] == 'LQ' else set())
    keys = set(values) | set(certified)
    if job['domain'] == 'FULL': keys.discard('half_y_symmetry'); allowed = allowed | {'half_y_symmetry'}
    differing = sorted(k for k in keys if not same(values.get(k, '<absent>'), certified.get(k, '<absent>')))
    unexpected = [k for k in differing if k not in allowed]
    mesh = Path(base.mesh_of(job))
    checks = dict(differing_keys=differing, unexpected_differences=unexpected,
                  mesh_sha256=r.sha(mesh), material_profile=values['material_profile'],
                  transport_supg=values['transport_supg'], supg_directional_tau=values['supg_directional_tau'],
                  half_y_symmetry=values.get('half_y_symmetry'))
    if job['model'] == 'NT' and (values['transport_supg'] != '.true.' or values['supg_directional_tau'] != '.false.'
                                 or values['supg_scale'] != '1.0'):
        unexpected.append('NT_SUPG_SETTING')
    if job['model'] == 'LQ' and values['transport_supg'] != '.false.': unexpected.append('LQ_TRANSPORT_SUPG')
    return checks


def prepare():
    data = release()
    base.NT_VARIANT = 'trace'
    inventory = r.read(base.OUTPUT/'mesh_inventory.json')
    for entry in inventory:
        a, counts = entry['radius_um'], tuple(entry['counts'])
        if r.sha(base.full_mesh(a, counts)) != entry['full_mesh_sha256'] or r.sha(base.half_mesh(a, counts)) != entry['half_mesh_sha256']:
            raise RuntimeError('Mesh changed since inventory: '+entry['tag'])
    certified_manifest = r.read(CERTIFIED/'manifest.json')
    if len({r.sha(base.half_mesh(3, base.FINEST)), certified_manifest['mesh_sha256'],
            data['nt_numerics']['r3_finest_half_mesh_sha256']}) != 1:
        raise RuntimeError('r=3 finest mesh differs from the certified trace run')
    certified = namelist((CERTIFIED/'input.nml').read_text(encoding='utf-8-sig'))
    INPUTS.mkdir(exist_ok=False)
    rows, problems = [], []
    for job in base.jobs():
        text, changed = preview(job)
        folder = INPUTS/job['name']; folder.mkdir()
        with (folder/'input.nml').open('x', encoding='ascii') as f: f.write(text)
        checks = check_against_certified(job, text, certified)
        r.save(folder/'preview_manifest.json', dict(job=job, changed_input_keys=changed, seed_nodes_expected=seed_location(job),
               input_sha256=r.sha(folder/'input.nml'), checks=checks, native_runs=0,
               note='Preview; the runner rebuilds the identical text at run time, with the seed short path if seeded.'))
        if checks['unexpected_differences']: problems.append((job['name'], checks['unexpected_differences']))
        rows.append(dict(name=job['name'], group=job['group'], purpose=PURPOSE[job['group']], model=job['model'],
                         radius_um=job['radius'], mesh=base.tag(job['radius'], job['counts']), domain=job['domain'],
                         voltage_V=job['voltage'], seed=(job['seed'][0]+':'+job['seed'][1]) if job['seed'] else 'COLD_START',
                         restart_voltage_factor=job['restart'], continuation_step=job['step'] or '',
                         supg='trace tau, scale 1' if job['model'] == 'NT' else 'off (LQ)',
                         mesh_sha256=checks['mesh_sha256'], input_sha256=r.sha(folder/'input.nml'),
                         differing_from_certified=' '.join(checks['differing_keys'])))
    with (INPUTS/'case_matrix.csv').open('x', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    r.save(INPUTS/'provenance.json', dict(
        release=str(RELEASE), release_sha256=r.sha(RELEASE), certified_input=str(CERTIFIED/'input.nml'),
        certified_input_sha256=r.sha(CERTIFIED/'input.nml'), case_count=len(rows), unexpected_differences=problems,
        all_inputs_match_certified_settings=not problems, native_runs=0, production_started=False,
        script_hashes={str(p): r.sha(p) for p in [Path(__file__).resolve(), Path(base.__file__).resolve()]},
        decision_files={p.name: r.sha(p) for p in RELEASE.parent.iterdir() if p.is_file()},
        case_matrix_sha256=r.sha(INPUTS/'case_matrix.csv')))
    print(dict(cases=len(rows), unexpected=problems, case_matrix_sha256=r.sha(INPUTS/'case_matrix.csv')))
    return 1 if problems else 0


def execute(groups, confirmed_sha):
    release()
    matrix = INPUTS/'case_matrix.csv'
    if r.sha(matrix) != confirmed_sha: raise RuntimeError('case_matrix.csv is not the reviewed version')
    if not r.read(INPUTS/'provenance.json')['all_inputs_match_certified_settings']: raise RuntimeError('Input checks failed')
    base.NT_VARIANT, base.RELEASE = 'trace', RELEASE
    built = base.make_input

    def checked(job, seed):
        text, changed = built(job, seed)
        reviewed = (INPUTS/job['name']/'input.nml').read_text(encoding='ascii')
        if job['seed']: text_cmp = r.replace(text, 'restart_nodes_file', "'"+SEED_PLACEHOLDER+"'")
        else: text_cmp = text
        if text_cmp != reviewed: raise RuntimeError('Runtime input differs from reviewed preview: '+job['name'])
        return text, changed
    base.make_input = checked
    r.save(base.OUTPUT/f'batch_{time.strftime("%Y%m%dT%H%M%S")}_trace_wrapper_identity.json', dict(
        wrapper=str(Path(__file__).resolve()), wrapper_sha256=r.sha(Path(__file__).resolve()),
        release_sha256=r.sha(RELEASE), case_matrix_sha256=confirmed_sha, groups=groups))
    return base.execute(groups)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=['prepare', 'execute'])
    p.add_argument('--groups', nargs='+', default=['lq_check', 'lq_sweep', 'nt_sweep', 'radius', 'radius_mesh'])
    p.add_argument('--confirm-matrix', help='SHA-256 of the reviewed case_matrix.csv')
    a = p.parse_args()
    if a.command == 'prepare': return prepare()
    if not a.confirm_matrix: p.error('execute requires --confirm-matrix')
    return execute(a.groups, a.confirm_matrix)


if __name__ == '__main__':
    raise SystemExit(main())
