"""Finish the runner's post-processing for a solved state whose post-processing crashed.

Used for H_r3_A48_R22_Z18_NT_V005 (batch 20261022T003858): the native solve and the flux/field
classification completed (summary state PASS), then the independent power audit raised
UnicodeDecodeError because audit_je_sensitivity_power_20261015.py reads manifest.json with the
locale codec (GBK). Run under PYTHONUTF8=1. Nothing is re-solved: this executes exactly the
remaining tail of run_je_production_20261021.run_native (heat-model check, physical power,
lossless sample compression, completion.json) and records the recovery.
"""
from pathlib import Path
import argparse
import sys
import time
import run_je_production_reduced_20261022 as reduced

base, r = reduced.base, reduced.r


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('name')
    p.add_argument('--reason', required=True)
    a = p.parse_args()
    if not sys.flags.utf8_mode: raise SystemExit('Run with PYTHONUTF8=1')
    job = base.by_name()[a.name]
    directory = base.OUTPUT/a.name
    if (directory/'completion.json').exists(): raise SystemExit('completion.json already exists')
    if (directory/'power_audit').exists(): raise SystemExit('power_audit exists; inspect manually')
    state, metrics = base.classify_run(directory)
    summary = r.read(directory/'summary.json')
    if state not in (base.PASS, base.ZERO_PASS) or summary['state'] != state:
        raise SystemExit(f'Not a solved accepted state: {state}')
    if not r.read(directory/'field_audit.json')['passed']: raise SystemExit('Field audit not passed')
    if metrics.get('heat_source') != 'E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE':
        raise RuntimeError('Wrong solved heat model')
    before = {f: r.sha(directory/f) for f in ('metrics.txt', 'nodes.csv', 'summary.json', 'input.nml', 'manifest.json')}
    completion = dict(state=state, fully_validated=False)
    if job['model'] == 'NT' and state == base.PASS:
        completion['physical_power'] = base.physical_power(directory, job)
    else:
        completion['power_audit'] = 'NOT_APPLIED' if state == base.PASS else 'ZERO_VOLTAGE'
    completion['field_samples'] = base.compress_samples(directory)
    r.save(directory/'completion.json', completion)
    after = {f: r.sha(directory/f) for f in before}
    if before != after: raise RuntimeError('Solved files changed during post-processing')
    r.save(directory/'postprocess_recovery.json', dict(
        recovered_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'), reason=a.reason, python_utf8_mode=True,
        native_solve_repeated=False, solved_file_hashes_unchanged=after, script=str(Path(__file__).resolve()),
        script_sha256=r.sha(Path(__file__).resolve()),
        note='Only the post-processing tail of run_native was executed; solve, metrics and classification are untouched.'))
    print(dict(name=a.name, state=state, completion=True))


if __name__ == '__main__':
    main()
