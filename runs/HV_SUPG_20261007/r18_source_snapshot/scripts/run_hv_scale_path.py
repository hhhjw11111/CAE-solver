"""Continue the SUPG coefficient at fixed voltage, preserving every attempted solve.

This changes only a numerical coefficient. Each point is re-equilibrated with the
original nonlinear/flux gates; neither a failed trial nor an interpolated state
can become a checkpoint. A new executable first replays the starting state.
"""
import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import subprocess
import sys

from hv_runtime import add_storage_arguments
from run_hv_sweep import ROOT, classify, sha


@dataclass
class ScalePath:
    current: float
    target: float
    maximum_step: float = .25
    minimum_step: float = .015625

    def __post_init__(self):
        if (not all(math.isfinite(x) for x in (self.current, self.target, self.maximum_step, self.minimum_step))
                or min(self.current, self.target, self.minimum_step) <= 0
                or self.maximum_step < self.minimum_step):
            raise ValueError('Positive finite scales and consistent step bounds are required')
        self.step = self.maximum_step

    @property
    def complete(self):
        return self.current == self.target

    def trial(self):
        sign = 1 if self.target > self.current else -1
        return self.current + sign * min(self.step, abs(self.target - self.current))

    def accept(self, trial):
        if trial != self.trial():
            raise ValueError('Only the proposed, solved coefficient may be accepted')
        self.current = trial
        self.step = min(self.maximum_step, self.step * 1.5)

    def reject(self):
        # Cut back the actual failed increment, including a short final step.
        attempted = abs(self.trial() - self.current)
        if attempted / 2 < self.minimum_step:
            return False
        self.step = attempted / 2
        return True


def starting_evidence(directory, mesh):
    manifest = json.loads((directory/'manifest.json').read_text(encoding='utf-8-sig'))
    execution = json.loads((directory/'execution.json').read_text(encoding='utf-8-sig'))
    fields = json.loads((directory/'field_audit.json').read_text(encoding='utf-8-sig'))
    state, metrics = classify(directory, execution['exit_code'])
    if state != 'SOLVE_FLUX_PASS_MESH_PENDING' or fields.get('passed') is not True:
        raise ValueError('The starting state must pass the existing solve/flux/field gates')
    if manifest.get('mesh_sha256') != sha(mesh):
        raise ValueError('The starting state must use the same mesh')
    if (manifest.get('model') != 'FROZEN_DRIFT_HEAT_ZNO_NT'
            or manifest.get('supg_scale') != metrics.get('supg_scale')
            or metrics.get('supg_directional_tau') != 'T'):
        raise ValueError('Starting model/stabilization evidence is inconsistent')
    voltage = metrics['voltage_V']
    if not math.isfinite(voltage) or not 0 < voltage <= 8 or manifest.get('voltages') != [voltage]:
        raise ValueError('One fully solved voltage in (0,8] is required')
    return manifest, voltage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--seed-run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--target-scale', type=float, required=True)
    parser.add_argument('--maximum-step', type=float, default=.25)
    parser.add_argument('--minimum-step', type=float, default=.015625)
    parser.add_argument('--maximum-attempts', type=int, default=24)
    parser.add_argument('--max-seconds-per-point', type=float, default=3600)
    parser.add_argument('--omp-threads', type=int, default=4)
    parser.add_argument('--mkl-threads', type=int, default=4)
    add_storage_arguments(parser)
    args = parser.parse_args()
    for name in ['exe', 'mesh', 'seed_run', 'output']:
        setattr(args, name, getattr(args, name).absolute())
    if (args.maximum_attempts < 1 or min(args.omp_threads, args.mkl_threads, args.ooc_core_mb) < 1
            or not math.isfinite(args.max_seconds_per_point) or args.max_seconds_per_point <= 0):
        raise ValueError('Invalid runtime limits')
    seed_manifest, voltage = starting_evidence(args.seed_run, args.mesh)
    path = ScalePath(seed_manifest['supg_scale'], args.target_scale, args.maximum_step, args.minimum_step)
    executable_hash, mesh_hash = sha(args.exe), sha(args.mesh)
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = dict(executable=str(args.exe), executable_sha256=executable_hash,
                    mesh=str(args.mesh), mesh_sha256=mesh_hash, voltage=voltage,
                    seed_run=str(args.seed_run), seed_nodes_sha256=sha(args.seed_run/'nodes.csv'),
                    starting_scale=path.current, target_scale=path.target,
                    maximum_step=path.maximum_step, minimum_step=path.minimum_step,
                    status='NUMERICAL_COEFFICIENT_PATH_NOT_SENSITIVITY_OR_MESH_CERTIFICATE')
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2))
    previous = args.seed_run
    records = []
    reason = 'ATTEMPT_LIMIT'
    # Attempt zero is always an ordinary fixed-voltage replay in this executable.
    for index in range(args.maximum_attempts + 1):
        scale = path.current if index == 0 else path.trial()
        dest = args.output/f'attempt_{index:03d}_scale_{scale:.9g}'
        command = [sys.executable, str(ROOT/'scripts/run_hv_point.py'), '--exe', str(args.exe),
                   '--mesh', str(args.mesh), '--output', str(dest), '--voltage', str(voltage),
                   '--seed-nodes', str(previous/'nodes.csv'), '--seed-voltage', str(voltage),
                   '--supg-scale', str(scale), '--omp-threads', str(args.omp_threads),
                   '--mkl-threads', str(args.mkl_threads), '--storage', args.storage,
                   '--ooc-core-mb', str(args.ooc_core_mb), '--max-seconds', str(args.max_seconds_per_point)]
        if args.split_holes:
            command.append('--split-holes')
        if sha(args.exe) != executable_hash or sha(args.mesh) != mesh_hash:
            raise RuntimeError('Executable or mesh changed during the path')
        print(f'START attempt={index} scale={scale:.9g} seed={previous.name}', flush=True)
        with (args.output/f'attempt_{index:03d}_runner.txt').open('xb') as stream:
            completed = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
        summary = json.loads((dest/'summary.json').read_text()) if (dest/'summary.json').exists() else {}
        state = summary.get('state', 'RUNNER_FAILED')
        accepted = completed.returncode == 0 and state == 'SOLVE_FLUX_PASS_MESH_PENDING'
        record = dict(attempt=index, scale=scale, directory=str(dest), seed=str(previous),
                      exit_code=completed.returncode, state=state, accepted=accepted)
        records.append(record)
        with (args.output/'attempts.jsonl').open('a', encoding='utf-8') as stream:
            stream.write(json.dumps(record) + '\n')
        print(f'END attempt={index} scale={scale:.9g} state={state}', flush=True)
        if accepted:
            previous = dest
            if index > 0:
                path.accept(scale)
            if path.complete:
                reason = 'TARGET_SOLVE_FLUX_PASS'; break
        elif index == 0:
            reason = 'STARTING_REPLAY_FAILED'; break
        elif not path.reject():
            reason = 'MINIMUM_STEP_FAILED'; break
    result = dict(status=reason, target_reached=reason == 'TARGET_SOLVE_FLUX_PASS',
                  last_accepted_scale=path.current, last_accepted_directory=str(previous),
                  attempts=records, mesh_validated=False, stabilization_sensitivity_validated=False,
                  complete_quantitative_0_8V=False)
    (args.output/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result['target_reached'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
