"""One-off, sequential evidence completion for this task; not a scheduler."""
from pathlib import Path
import json
import subprocess
import sys
import time
import numpy as np

ROOT = Path(__file__).absolute().parents[3]
RUN = Path(__file__).absolute().parent
SCRIPTS = ROOT / 'solver_v2/scripts'
sys.path.insert(0, str(SCRIPTS))
from run_hv_sweep import read_metrics, sha


def main():
    start = time.monotonic()
    while not (RUN / 'r17_scale2_global_V060/execution.json').exists():
        if time.monotonic() - start > 2000:
            raise TimeoutError('Global-limiter comparison has not exited')
        time.sleep(2)
    commands = [
        ('coarse_sensitivity_global', [sys.executable, str(SCRIPTS / 'audit_hv_stabilization.py'),
            '--runs', str(RUN / 'r17_scale05_V060'), str(RUN / 'r17_scale1_V060'),
            str(RUN / 'r17_scale2_global_V060'), '--output',
            str(RUN / 'stabilization_coarse_global_V060_r17.json')]),
        ('production', [sys.executable, str(SCRIPTS / 'check_hv_production_regression.py'),
            '--baseline', str(ROOT / 'solver_v2/build/HV_BASELINE_20261007/TPSolverV2.exe'),
            '--candidate', str(ROOT / 'solver_v2/build/HV_SUPG_20261008_R17/TPSolverV2.exe'),
            '--case', str(ROOT / 'solver_v2/cases/F1_FAST_HOLE_L0/FAST_H1_LQ_L0.nml'),
            '--output', str(RUN / 'production_r17_serial')]),
    ]
    records = []
    for label, command in commands:
        print('START', label, flush=True)
        start = time.monotonic()
        with (RUN / f'final_queue_r17_{label}.txt').open('xb') as f:
            result = subprocess.run(command, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
        records.append(dict(label=label, command=command, exit_code=result.returncode,
                            seconds=time.monotonic()-start))
        (RUN / 'final_queue_r17_execution.json').write_text(json.dumps(records, indent=2))
        print('END', label, result.returncode, flush=True)

    prod = RUN / 'production_r17_serial'
    case = ROOT / 'solver_v2/cases/F1_FAST_HOLE_L0/FAST_H1_LQ_L0.nml'
    exe_paths = [ROOT / f'solver_v2/build/{b}/TPSolverV2.exe' for b in
                 ['HV_BASELINE_20261007', 'HV_SUPG_20261008_R17']]
    inputs = [case, case.parent / 'H1_L0.mesh', *exe_paths]
    fields = [np.genfromtxt(prod / k / 'nodes.csv', delimiter=',', names=True)
              for k in ['baseline', 'candidate']]
    metrics = [read_metrics(prod / k / 'metrics.txt') for k in ['baseline', 'candidate']]
    same_topology = all(np.array_equal(fields[0][k], fields[1][k])
                        for k in ['node', 'x_m', 'y_m', 'z_m'])
    scales = dict(ux_m=1e-10, uy_m=1e-10, uz_m=1e-10, phi_V=3., eta_n=.1, eta_p=.1, theta_K=10.)
    differences = {k: dict(max_absolute=float(np.max(np.abs(fields[0][k]-fields[1][k]))),
                          scaled_by_frozen_input=float(np.max(np.abs(fields[0][k]-fields[1][k]))/s))
                   for k, s in scales.items()}
    metric_differences = {k: dict(baseline=metrics[0][k], candidate=metrics[1][k],
                               absolute_difference=abs(metrics[0][k]-metrics[1][k]))
                          for k in ['theta_max_K', 'joule_W', 'electron_right_A',
                                    'heat_balance_relative', 'electron_balance_relative',
                                    'hole_balance_diagnostic']}
    report = json.loads((prod / 'report.json').read_text())
    result = dict(strict_byte_identical=report['passed'], same_topology=same_topology,
                  old_frozen_acceptance=[m['case_acceptance'] for m in metrics],
                  fields=differences, metrics=metric_differences,
                  maximum_scaled_node_difference=max(v['scaled_by_frozen_input'] for v in differences.values()),
                  provenance=[dict(path=str(p), sha256=sha(p)) for p in inputs],
                  interpretation='Numerical differences are separate from the unchanged strict-byte verdict; this is one representative LQ production case.')
    with (prod / 'numerical_comparison.json').open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
