"""Reduced 30-case trace-metric JE production (user request 2026-10-22: 算例减少).

Importing this module patches the shared base runner in place, so every script that imports
it (figures, estimate comparison, audits) sees the reduced job list and the new output root.
Nothing is solved at import. Physics, tolerances, gates, SUPG variant, meshes, seeds logic and
the observable classification are those of the 2026-10-22 release; only the voltage lists shrink:
  LQ r=3: 0/2/4/6/8 V (cold starts)
  NT r=3: 0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8 V, each seeded from the previous accepted voltage with
          the same absolute 0.5 V continuation step as the 49-case matrix
  lq_check, radius and radius_mesh groups unchanged.
The stopped 49-case batch in stage5_production_20261021 is kept and never mixed with this one.
"""
from pathlib import Path
import argparse
import json
import shutil
import run_je_production_20261021 as base
import run_je_production_20261022 as trace

r = base.r
OLD_OUTPUT = base.OUTPUT
OUTPUT = r.STUDY/'stage5_production_reduced_20261022'
INPUTS = r.STUDY/'paper_scope/production_inputs_trace_reduced_20261022'
LQ_VOLTAGES = [0., 2., 4., 6., 8.]
NT_VOLTAGES = [0., .5, 1., 2., 3., 4., 5., 6., 7., 8.]
AUTHORIZATION = 'user chat 2026-10-22: 算例减少 把所有要求打包md 我让codex执行 确保任务完整'


def jobs():
    result = []

    def add(model, a, counts, v, domain='Y_HALF', seed=None, restart=0., step=None, group=''):
        name = f'{"F" if domain == "FULL" else "H"}_{base.tag(a, counts)}_{model}_V{round(v*10):03d}'
        result.append(dict(name=name, model=model, radius=a, counts=list(counts), voltage=v, domain=domain,
                           seed=seed, restart=restart, step=step, group=group))
        return name
    add('LQ', 3, base.COARSE, 8, 'FULL', group='lq_check')
    add('LQ', 3, base.COARSE, 8, group='lq_check')
    add('LQ', 3, base.MID, 8, group='lq_check')
    for v in LQ_VOLTAGES: add('LQ', 3, base.FINEST, v, group='lq_sweep')
    previous, previous_v = None, None
    for v in NT_VOLTAGES:
        if v <= .5: previous = add('NT', 3, base.FINEST, v, step=1.0, group='nt_sweep')
        else: previous = add('NT', 3, base.FINEST, v, seed=('job', previous), restart=previous_v/v, step=.5/v, group='nt_sweep')
        previous_v = v
    source6 = f'H_{base.tag(3, base.FINEST)}_NT_V060'
    for a in (2, 4):
        add('LQ', a, base.FINEST, 6, group='radius'); add('LQ', a, base.FINEST, 8, group='radius')
        add('LQ', a, base.MID, 8, group='radius_mesh')
        six = add('NT', a, base.FINEST, 6, seed=('radial', source6), restart=1.0, group='radius')
        eight = add('NT', a, base.FINEST, 8, seed=('job', six), restart=.75, step=.25, group='radius')
        add('NT', a, base.MID, 8, seed=('transfer', eight), restart=1.0, group='radius_mesh')
    return result


# Patch the shared modules: base.execute/by_name/audits/table and trace.prepare/execute look these up at call time.
base.jobs = jobs
base.OUTPUT = OUTPUT
base.NT_VARIANT = 'trace'
base.RELEASE = trace.RELEASE
trace.INPUTS = INPUTS


def stage_meshes():
    """Copy the already generated meshes and inventory; every file is re-hashed against the inventory."""
    OUTPUT.mkdir(exist_ok=True)
    if not (OUTPUT/'mesh_inventory.json').exists():
        shutil.copytree(OLD_OUTPUT/'meshes', OUTPUT/'meshes')
        shutil.copy2(OLD_OUTPUT/'mesh_inventory.json', OUTPUT/'mesh_inventory.json')
    for entry in r.read(OUTPUT/'mesh_inventory.json'):
        a, counts = entry['radius_um'], tuple(entry['counts'])
        if (r.sha(base.full_mesh(a, counts)) != entry['full_mesh_sha256']
                or r.sha(base.half_mesh(a, counts)) != entry['half_mesh_sha256']):
            raise RuntimeError('Copied mesh differs from inventory: '+entry['tag'])


def prepare():
    stage_meshes()
    status = trace.prepare()
    reduction = dict(authorization=AUTHORIZATION, previous_matrix=str(trace.r.STUDY/'paper_scope/production_inputs_trace_20261022/case_matrix.csv'),
                     previous_batch_stop=str(OLD_OUTPUT/'batch_20261022T000632_user_stop.json'),
                     lq_voltages=LQ_VOLTAGES, nt_voltages=NT_VOLTAGES, case_count=len(jobs()),
                     unchanged=['physics', 'tolerances', 'gates', 'SUPG trace variant', 'meshes', 'radius groups',
                                'lq_check', 'observable classification', 'hole exclusion'],
                     note='Only solved states are plotted; dropped voltages are never interpolated or scaled.')
    r.save(INPUTS/'matrix_reduction.json', reduction)
    return status


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=['prepare', 'execute', 'audit', 'table', 'compare'])
    p.add_argument('--groups', nargs='+', default=['lq_check', 'lq_sweep', 'nt_sweep', 'radius', 'radius_mesh'])
    p.add_argument('--confirm-matrix', help='SHA-256 of the user-approved reduced case_matrix.csv')
    a = p.parse_args()
    if a.command == 'prepare': return prepare()
    if a.command == 'execute':
        if not a.confirm_matrix: p.error('execute requires --confirm-matrix')
        stage_meshes()
        return trace.execute(a.groups, a.confirm_matrix)
    if a.command == 'audit': base.audits(); return 0
    if a.command == 'table': base.table(); return 0
    import compare_je_production_estimates_20261022 as compare
    import sys
    sys.argv = [sys.argv[0]]
    return compare.main()


if __name__ == '__main__':
    raise SystemExit(main())
