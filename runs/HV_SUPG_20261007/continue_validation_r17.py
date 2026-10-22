"""Orchestrate this run's remaining comparisons; child evidence is authoritative.

Run fine and coarse modes concurrently only with the existing small/large jobs
shown in this task. Each lane is sequential, keeping at most one large and one
small matrix factorization active. No old result directory is reused.
"""
from pathlib import Path
import argparse, csv, json, subprocess, sys, time

ROOT = Path(__file__).absolute().parents[3]
RUN = ROOT / 'solver_v2/runs/HV_SUPG_20261007'
EXE = ROOT / 'solver_v2/build/HV_SUPG_20261008_R17/TPSolverV2.exe'
MESH = {
    'fine': ROOT / 'solver_v2/cases/HV_SURFACE_20261008/XY1_Z3/H1_XY1_Z3.mesh',
    'coarse': ROOT / 'solver_v2/cases/HV_SURFACE_20261008/XY0_Z2/H1_XY0_Z2.mesh',
}
SCRIPT = ROOT / 'solver_v2/scripts'


def wait_record(path):
    start = time.monotonic()
    while time.monotonic() - start < 7400:
        try:
            return json.loads(path.read_text(encoding='utf-8-sig'))
        except (FileNotFoundError, json.JSONDecodeError):
            time.sleep(2)
    raise TimeoutError(str(path))


def point(mesh, output, nodes, voltage, scale, threads):
    return [sys.executable, str(SCRIPT/'run_hv_point.py'), '--exe', str(EXE),
            '--mesh', str(mesh), '--output', str(RUN/output), '--voltage', str(voltage),
            '--seed-nodes', str(nodes), '--seed-voltage', str(voltage),
            '--supg-scale', str(scale), '--omp-threads', str(threads), '--mkl-threads', str(threads)]


def main():
    p = argparse.ArgumentParser(); p.add_argument('mode', choices=['fine', 'coarse']); a = p.parse_args()
    records = []
    def execute(label, command):
        print('START', a.mode, label, flush=True)
        start = time.monotonic()
        with (RUN/f'queue_r17_{a.mode}_{label}.txt').open('xb') as f:
            result = subprocess.run(command, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
        records.append(dict(label=label, command=command, exit_code=result.returncode,
                            seconds=time.monotonic()-start))
        (RUN/f'queue_r17_{a.mode}_execution.json').write_text(json.dumps(dict(
            scope='ORCHESTRATION_ONLY_CHILD_RESULTS_ARE_AUTHORITATIVE', steps=records), indent=2))
        print('END', a.mode, label, 'exit', result.returncode, flush=True)

    if a.mode == 'fine':
        trace = RUN/'r17_inverse_XY1_Z3'
        wait_record(trace/'execution.json')
        with (trace/'inverse_branch.csv').open(newline='') as f:
            points = {int(r['point']): r for r in csv.DictReader(f)}
        for i in [3, 4]:
            if i not in points:
                print('MISSING_ACCEPTED_POINT', i, flush=True); continue
            v = float(points[i]['voltage_V'])
            execute(f'replay_{i}', point(MESH['fine'], f'r17_fine_replay_branch_{i:02d}',
                    trace/f'accepted_{i:05d}_nodes.csv', v, 1, 2))
        for scale, label in [(.5, '05'), (2, '2')]:
            execute(f'scale{label}_V060', point(MESH['fine'], f'r17_fine_scale{label}_V060',
                    RUN/'r17_XY1_Z3_V060/nodes.csv', 6, scale, 4))
        execute('sensitivity', [sys.executable, str(SCRIPT/'audit_hv_stabilization.py'), '--runs',
                str(RUN/'r17_fine_scale05_V060'), str(RUN/'r17_XY1_Z3_V060'),
                str(RUN/'r17_fine_scale2_V060'), '--output', str(RUN/'stabilization_fine_V060_r17.json')])
    else:
        wait_record(RUN/'r17_inverse_scale05/execution.json')
        execute('scale2_trace', [sys.executable, str(SCRIPT/'run_hv_temperature_trace.py'),
                '--exe', str(EXE), '--mesh', str(MESH['coarse']), '--output', str(RUN/'r17_inverse_scale2'),
                '--seed-nodes', str(RUN/'r13_XY0_Z2_V060/nodes.csv'), '--seed-voltage', '6',
                '--step-k', '5.4', '--points', '3', '--supg-scale', '2',
                '--omp-threads', '2', '--mkl-threads', '2', '--max-seconds', '5400'])
        for scale, label in [(.5, '05'), (2, '2')]:
            seed = RUN/f'r17_inverse_scale{label}/accepted_00001_nodes.csv'
            if not seed.is_file():
                print('MISSING_ACCEPTED_SEED', label, flush=True); continue
            execute(f'scale{label}_V060', point(MESH['coarse'], f'r17_scale{label}_V060', seed, 6, scale, 2))
        execute('sensitivity', [sys.executable, str(SCRIPT/'audit_hv_stabilization.py'), '--runs',
                str(RUN/'r17_scale05_V060'), str(RUN/'r17_scale1_V060'), str(RUN/'r17_scale2_V060'),
                '--output', str(RUN/'stabilization_coarse_V060_r17.json')])
    return 0 if records and all(r['exit_code'] == 0 for r in records) else 1


if __name__ == '__main__':
    raise SystemExit(main())
