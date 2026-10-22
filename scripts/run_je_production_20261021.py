"""Observable-scoped JE production: 0--8 V LQ/NT at r=3 um and r=2/3/4 at 6/8 V.

Released by paper_scope/production_release_scoped_20261021/release.json. The
frozen executable, physics, tolerances and native flux gates are unchanged.
This runner only prepares Y-half meshes and inputs, sets continuation step
sizes (a path setting, not a tolerance) and executes one native solver at a
time. Every state is solved at its own voltage; nothing is interpolated or
scaled. LQ writes no continuation factors, so it is classified with the same
numeric gates minus those two fields. Zero-voltage states are checked as
equilibria with absolute limits because relative balances divide by zero.
"""
from pathlib import Path
import argparse
import csv
import gzip
import json
import math
import shutil
import subprocess
import sys
import time
import traceback
import run_je_reconstruction_mesh_20261018 as r
import run_hv_sweep as sweep
from run_je_paper_ooc8000_20261020 import preflight
from generate_je_paper_radius_mesh_20261018 import write as write_radius_mesh
from audit_hv_mesh import LIMITS, relative
from audit_hv_field_shape import audit as field_shape_audit
from audit_je_sensitivity_power_20261015 import run as power_audit
import audit_hv_half_y

OUTPUT = r.STUDY/'stage5_production_20261021'
RELEASE = r.STUDY/'paper_scope/production_release_scoped_20261021/release.json'
CANDIDATES = r.STUDY/'paper_scope/voltage_radius_scope_20261020_v2/candidate_inputs'
FINEST, MID, COARSE = (48, 22, 18), (40, 20, 16), (32, 18, 14)
OOC_MB = 8000
# Evidence the release relies on: regenerated r=3 meshes must be byte-identical.
KNOWN_FULL = {(3, FINEST): '0e24e5515b3b6942de571b14011841878e47cfd61efc64b006975b9a8742fcf0'}
KNOWN_HALF = {(3, FINEST): '3905ef1c2d9e333b880a3f347c088bf2d02f4e9c1eca595f0efbae762659f733'}
MESHES = [(3, FINEST), (3, MID), (3, COARSE), (2, FINEST), (2, MID), (4, FINEST), (4, MID)]
ZERO_LIMITS = dict(joule_W=1e-18, theta_max_K=1e-9, electron_right_A=1e-15)
GATES = {'maximum_stage_balance': 1e-8, 'heat_balance_relative': .01, 'electron_balance_relative': .01,
         'hole_balance_diagnostic': .01, 'joule_quadrature_4_5_relative': .001}
HALF_GATES = ['physical_exterior_heat_balance_relative', 'physical_exterior_electron_balance_relative',
              'physical_exterior_hole_balance_diagnostic', 'symmetry_heat_flux_relative',
              'symmetry_electron_flux_relative', 'symmetry_hole_flux_diagnostic']
PASS = r.PASS
ZERO_PASS = 'ZERO_VOLTAGE_EQUILIBRIUM_PASS'
TRACE = r.STUDY/'stage5_controlled_numerics_20261021/trace_joint_supg'
# 'trace': supg_directional_tau=.false. with the checked hole-block solve, exactly as
# the controlled trace study; 'directional': the frozen candidate template.
NT_VARIANT = 'directional'


def tag(a, counts):
    return f'r{a}_A{counts[0]}_R{counts[1]}_Z{counts[2]}'


def full_mesh(a, counts):
    return OUTPUT/'meshes/full'/tag(a, counts)/'H1.mesh'


def half_mesh(a, counts):
    return OUTPUT/'meshes/half'/tag(a, counts)/'H1_Y_HALF.mesh'


def prepare_meshes():
    inventory = []
    for a, counts in MESHES:
        folder = full_mesh(a, counts).parent
        if not folder.exists(): write_radius_mesh(folder, a, *counts)
        half = half_mesh(a, counts).parent
        if not half.exists():
            subprocess.run([sys.executable, str(r.ROOT/'scripts/generate_hv_half_y.py'), '--mesh',
                            str(full_mesh(a, counts)), '--output', str(half)], check=True, stdout=subprocess.DEVNULL)
        f, h = r.sha(full_mesh(a, counts)), r.sha(half_mesh(a, counts))
        if KNOWN_FULL.get((a, counts), f) != f or KNOWN_HALF.get((a, counts), h) != h:
            raise RuntimeError('Regenerated reference mesh differs: '+tag(a, counts))
        hm = r.read(half/'manifest.json')
        inventory.append(dict(tag=tag(a, counts), radius_um=a, counts=list(counts), full_mesh_sha256=f,
                              half_mesh_sha256=h, half_nodes=hm['nodes'], half_elements=hm['elements'],
                              minimum_sampled_jacobian=hm['minimum_sampled_jacobian'],
                              matches_released_reference=(a, counts) in KNOWN_HALF))
    return inventory


def jobs():
    result = []

    def add(model, a, counts, v, domain='Y_HALF', seed=None, restart=0., step=None, group=''):
        name = f'{"F" if domain == "FULL" else "H"}_{tag(a, counts)}_{model}_V{round(v*10):03d}'
        result.append(dict(name=name, model=model, radius=a, counts=list(counts), voltage=v, domain=domain,
                           seed=seed, restart=restart, step=step, group=group))
        return name
    add('LQ', 3, COARSE, 8, 'FULL', group='lq_check')
    add('LQ', 3, COARSE, 8, group='lq_check')
    add('LQ', 3, MID, 8, group='lq_check')
    for i in range(17): add('LQ', 3, FINEST, i/2, group='lq_sweep')
    previous = None
    for i in range(17):
        v = i/2
        if i <= 1: previous = add('NT', 3, FINEST, v, step=1.0, group='nt_sweep')
        else: previous = add('NT', 3, FINEST, v, seed=('job', previous), restart=(v-.5)/v, step=.5/v, group='nt_sweep')
    source6 = f'H_{tag(3, FINEST)}_NT_V060'
    for a in (2, 4):
        add('LQ', a, FINEST, 6, group='radius'); add('LQ', a, FINEST, 8, group='radius')
        add('LQ', a, MID, 8, group='radius_mesh')
        six = add('NT', a, FINEST, 6, seed=('radial', source6), restart=1.0, group='radius')
        eight = add('NT', a, FINEST, 8, seed=('job', six), restart=.75, step=.25, group='radius')
        add('NT', a, MID, 8, seed=('transfer', eight), restart=1.0, group='radius_mesh')
    return result


def by_name():
    return {j['name']: j for j in jobs()}


def mesh_of(job):
    a, counts = job['radius'], tuple(job['counts'])
    return full_mesh(a, counts) if job['domain'] == 'FULL' else half_mesh(a, counts)


def prepare_seed(job):
    if not job['seed']: return None
    kind, source = job['seed']
    folder = OUTPUT/source
    state = classify_run(folder)[0]
    if state != PASS: raise RuntimeError('Seed is not an accepted state: '+source)
    if kind == 'job': return folder/'nodes.csv'
    target = OUTPUT/'seeds'/job['name']
    if not target.exists():
        target.parent.mkdir(exist_ok=True)
        src = by_name()[source]
        if kind == 'radial':
            command = [sys.executable, str(r.ROOT/'scripts/map_radius_initial_guess_20261021.py'),
                       '--source-mesh', str(mesh_of(src)), '--source-nodes', str(folder/'nodes.csv'),
                       '--source-radius', f'{src["radius"]}e-6', '--target-mesh', str(mesh_of(job)),
                       '--target-radius', f'{job["radius"]}e-6', '--output', str(target)]
        else:
            command = [sys.executable, str(r.ROOT/'scripts/transfer_hv_initial_guess.py'), '--source-mesh',
                       str(mesh_of(src)), '--source-nodes', str(folder/'nodes.csv'), '--target-mesh',
                       str(mesh_of(job)), '--output', str(target)]
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
    report = r.read(target/'transfer_report.json')
    if report['target_mesh_sha256'] != r.sha(mesh_of(job)) or report['source_nodes_sha256'] != r.sha(folder/'nodes.csv'):
        raise RuntimeError('Seed mapping identity mismatch: '+job['name'])
    return target/'nodes.csv'


def template(job):
    # Only the zero-voltage NT candidate exists; every changed key is rewritten below.
    return CANDIDATES/('R03_LQ_V080' if job['model'] == 'LQ' else 'R03_NT_V000')/'input.nml'


def make_input(job, seed):
    text = (template(job)).read_text(encoding='utf-8-sig')
    changes = dict(case_id=f"'JE_PROD_{job['name']}'", mesh_file="'"+r.short(mesh_of(job))+"'",
                   target_voltage=repr(-float(job['voltage'])), hole_radius=repr(float(f'{job["radius"]}e-6')),
                   restart_nodes_file="'"+(r.short(seed) if seed else '')+"'",
                   restart_voltage_factor=repr(float(job['restart'])),
                   restart_joule_factor='1.0' if seed else '0.0')
    if job['step']:
        changes.update(continuation_initial_step=repr(job['step']), continuation_maximum_step=repr(job['step']))
    if job['model'] == 'NT' and NT_VARIANT == 'trace': changes['supg_directional_tau'] = '.false.'
    for key, value in changes.items(): text = r.replace(text, key, value)
    required = ['heat_source_current_work=.true.', 'residual_tolerance=1e-9,update_tolerance=1e-9',
                'target_joule_factor=1.0', "material_profile='ZNO_F1_"+job['model']+"'",
                'transport_supg=.true.,supg_scale=1.0,supg_directional_tau='+('.false.' if NT_VARIANT == 'trace' else '.true.')
                if job['model'] == 'NT' else 'transport_supg=.false.']
    for literal in required:
        if literal not in text: raise RuntimeError('Frozen input requirement missing: '+literal)
    if text.count('&problem_control') != 1 or 'half_y_symmetry' in text: raise RuntimeError('Unexpected template')
    if job['domain'] == 'Y_HALF': text = text.replace('&problem_control', '&problem_control\n half_y_symmetry=.true.,')
    return text, sorted(changes)


def classify_run(directory, rc=None):
    directory = Path(directory)
    if not (directory/'execution.json').exists(): return 'NOT_RUN', {}
    if rc is None: rc = r.read(directory/'execution.json')['exit_code']
    m = sweep.read_metrics(directory/'metrics.txt')
    if m.get('voltage_V') == 0:
        keys = ['maximum_accepted_residual', 'maximum_accepted_update', 'maximum_stage_balance', *ZERO_LIMITS]
        if not (directory/'nodes.csv').exists() or not all(isinstance(m.get(k), float) and math.isfinite(m[k]) for k in keys):
            return 'FAILED_NO_COMPLETE_STATE', m
        ok = (m['maximum_accepted_residual'] < 1e-9 and m['maximum_accepted_update'] < 1e-9 and
              m['maximum_stage_balance'] < 1e-8 and all(abs(m[k]) < lim for k, lim in ZERO_LIMITS.items()))
        return (ZERO_PASS if ok else 'ZERO_VOLTAGE_EQUILIBRIUM_FAIL'), m
    if m.get('model') != 'LQ': return sweep.classify(directory, rc)
    required = ['maximum_accepted_residual', 'maximum_accepted_update', 'theta_max_K', 'joule_W']
    complete = m.get('solution_mode') == 'THREE_STAGE_LQ_NO_RETURN_FEEDBACK' and (directory/'nodes.csv').exists()
    complete = complete and all(isinstance(m.get(k), float) and math.isfinite(m[k]) for k in required)
    complete = complete and 0 <= m['maximum_accepted_residual'] < 1e-9 and 0 <= m['maximum_accepted_update'] < 1e-9
    if not complete: return 'FAILED_NO_COMPLETE_STATE', m
    limits = dict(GATES)
    if m.get('domain_symmetry') == 'Y_HALF': limits.update({k: .01 for k in HALF_GATES})
    gates = all(isinstance(m.get(k), float) and math.isfinite(m[k]) and 0 <= m[k] < lim for k, lim in limits.items())
    if rc == 0 and m.get('case_acceptance') == 'PASS_NOT_MESH_ACCEPTANCE' and gates: return PASS, m
    return 'SOLVED_FLUX_GATE_FAILED', m


def accepted(directory):
    return classify_run(directory)[0] in (PASS, ZERO_PASS) and (Path(directory)/'completion.json').exists()


def run_native(job, gate):
    preflight()
    directory = OUTPUT/job['name']
    if shutil.disk_usage(OUTPUT).free < r.START_DISK: raise RuntimeError('RESOURCE_STOP: fewer than 16 GiB free')
    seed = prepare_seed(job)
    text, changed = make_input(job, seed)
    directory.mkdir(exist_ok=False)
    inp = directory/'input.nml'
    inp.write_text(text, encoding='ascii')
    env = r.runtime_env()
    env.update(OMP_NUM_THREADS=str(r.THREADS), MKL_NUM_THREADS=str(r.THREADS))
    trace = job['model'] == 'NT' and NT_VARIANT == 'trace'
    env = r.configure_storage(env, Path(r.short(directory)), 'OOC', OOC_MB, trace)
    mesh = mesh_of(job)
    r.save(directory/'manifest.json', dict(executable=str(r.EXE), executable_sha256=r.EXPECTED_EXE, mesh=str(mesh),
        mesh_sha256=r.sha(mesh), model='JE_CURRENT_WORK_DIAGNOSTIC_ZNO_'+job['model'], model_contract='JE_SIMPLIFIED_R0_V1',
        voltages=[job['voltage']], radius_um=job['radius'], counts=job['counts'], group=job['group'],
        supg_scale=1.0 if job['model'] == 'NT' else None, domain_symmetry=job['domain'],
        supg_directional_tau=(not trace) if job['model'] == 'NT' else None,
        nt_variant=NT_VARIANT if job['model'] == 'NT' else None, checked_hole_block_solve=trace,
        log_weighted_quadrature_order=3, runtime_env=r.recorded_runtime(env),
        seed=dict(kind=job['seed'][0], source=job['seed'][1], nodes=str(seed), sha256=r.sha(seed)) if seed else None,
        restart_voltage_factor=job['restart'], continuation_step_override=job['step'],
        path_setting_note='Continuation step size changes the load path only; tolerances and gates are frozen.',
        initial_guess_status='REQUIRES_REEQUILIBRATION' if seed else 'COLD_START', heat_source_unchanged=True,
        input_template=str(template(job)), changed_input_keys=changed,
        wall_time_limit_seconds=r.TIMEOUT, release=str(RELEASE), release_sha256=r.sha(RELEASE), verification_gate=gate,
        resource_limits=dict(start_disk_bytes=r.START_DISK, runtime_stop_disk_bytes=r.STOP_DISK,
                             stop_available_ram_bytes=r.MIN_AVAILABLE_RAM, low_ram_consecutive_samples=3)))
    started = time.monotonic()
    reason, process, monitor = None, None, None
    peak = dict(peak_rss_bytes=0, peak_private_bytes=0, peak_pagefile_bytes=0, peak_scratch_bytes=0,
                minimum_disk_free_bytes=shutil.disk_usage(directory).free, minimum_available_physical_bytes=None)
    try:
        with (directory/'stdout.txt').open('xb') as log, (directory/'resource_samples.jsonl').open('x', encoding='ascii') as samples:
            process = subprocess.Popen([r.short(r.EXE), 'run-f1', r.short(inp), r.short(directory)], cwd=r.short(r.ROOT),
                                       env=env, stdout=log, stderr=subprocess.STDOUT)
            r.save(directory/'process.json', dict(pid=process.pid, executable=str(r.EXE), owned_by_this_run=True))
            monitor = r.WindowsMonitor(process.pid)
            low_ram = 0
            while process.poll() is None:
                try:
                    sample = monitor.sample()
                except OSError:
                    if process.poll() is not None: break
                    raise
                sample.update(elapsed_seconds=time.monotonic()-started, disk_free_bytes=shutil.disk_usage(directory).free,
                              scratch_bytes=r.scratch_bytes(directory/'pardiso_ooc'))
                samples.write(json.dumps(sample)+'\n'); samples.flush()
                for key, source in [('peak_rss_bytes', 'peak_rss_bytes'), ('peak_private_bytes', 'private_bytes'),
                                    ('peak_pagefile_bytes', 'peak_pagefile_bytes'), ('peak_scratch_bytes', 'scratch_bytes')]:
                    peak[key] = max(peak[key], sample[source])
                peak['minimum_disk_free_bytes'] = min(peak['minimum_disk_free_bytes'], sample['disk_free_bytes'])
                previous = peak['minimum_available_physical_bytes']
                peak['minimum_available_physical_bytes'] = min(previous if previous is not None else sample['available_physical_bytes'],
                                                               sample['available_physical_bytes'])
                low_ram = low_ram+1 if sample['available_physical_bytes'] < r.MIN_AVAILABLE_RAM else 0
                if sample['elapsed_seconds'] >= r.TIMEOUT: reason = 'HARD_TIMEOUT_7200_SECONDS'
                elif sample['disk_free_bytes'] < r.STOP_DISK: reason = 'RESOURCE_STOP_DISK_BELOW_4_GIB'
                elif low_ram >= 3: reason = 'RESOURCE_STOP_AVAILABLE_RAM_BELOW_512_MIB'
                if reason:
                    print(reason+'; stopping only owned PID '+str(process.pid), flush=True)
                    process.kill(); break
                time.sleep(1)
            process.wait()
    except BaseException as exc:
        reason = reason or 'MONITOR_OR_RUNNER_EXCEPTION: '+repr(exc)
        if process is not None and process.poll() is None:
            process.kill(); process.wait()
        raise
    finally:
        if monitor is not None: monitor.close()
        if process is not None:
            peak['peak_scratch_bytes'] = max(peak['peak_scratch_bytes'], r.scratch_bytes(directory/'pardiso_ooc'))
            r.save(directory/'execution.json', dict(exit_code=process.returncode, seconds=time.monotonic()-started,
                   timed_out=reason == 'HARD_TIMEOUT_7200_SECONDS', stop_reason=reason, input_sha256=r.sha(inp),
                   executable_sha256=r.EXPECTED_EXE))
            r.save(directory/'resource_usage.json', dict(**peak, sampling_seconds=1, scope='Native process only'))
            r.cleanup_owned_factors(directory, process)
    state, metrics = classify_run(directory, process.returncode)
    if (directory/'nodes.csv').exists():
        fields = r.audit_nodes(directory/'nodes.csv', mesh, job['voltage'], half_y=job['domain'] == 'Y_HALF')
        r.save(directory/'field_audit.json', fields)
        if not fields['passed']: state = 'FAILED_FIELD_AUDIT'
    if reason: state = 'RESOURCE_OR_TIMEOUT_STOP'
    r.save(directory/'summary.json', dict(state=state, mesh_validated=False, power_gap_validated=False,
                                          elapsed_seconds=time.monotonic()-started, stop_reason=reason))
    if state not in (PASS, ZERO_PASS): return state
    if metrics.get('heat_source') != 'E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE':
        raise RuntimeError('Wrong solved heat model')
    completion = dict(state=state, fully_validated=False)
    if job['model'] == 'NT' and state == PASS:
        completion['physical_power'] = physical_power(directory, job)
    else:
        completion['power_audit'] = ('NOT_APPLIED: the independent power assembler implements exponential '
                                     'carriers with SUPG; native independent flux gates apply') if state == PASS else 'ZERO_VOLTAGE'
    completion['field_samples'] = compress_samples(directory)
    r.save(directory/'completion.json', completion)
    return state


def physical_power(directory, job):
    power_audit(directory)
    b = r.read(directory/'power_audit/order5.json')
    met = sweep.read_metrics(directory/'metrics.txt')
    if job['domain'] == 'Y_HALF':
        cut = b['boundary_sets']['y-']['power_W']
        physical = b['all_exterior_power_W']-cut
        reference = met['physical_exterior_terminal_power_W']
        factor = 2
    else:
        cut, physical, reference, factor = 0., b['all_exterior_power_W'], met['terminal_power_W'], 1
    difference = abs(physical-reference)/max(abs(physical), abs(reference), 1e-300)
    if difference >= 1e-8: raise RuntimeError('Independent physical-exterior power does not match native audit')
    result = dict(domain_symmetry=job['domain'], voltage_V=met['voltage_V'], raw_JE_W=b['JE_W'],
                  raw_symmetry_cut_power_W=cut, raw_physical_exterior_power_W=physical,
                  physical_exterior_gap_relative=(b['JE_W']-physical)/b['JE_W'],
                  reconstructed_full_JE_W=factor*b['JE_W'], independent_vs_native_relative=difference,
                  comparison_limit=1e-8, power_gap_acceptance_threshold=None,
                  note='Gap is diagnostic; never heat compensation.')
    r.save(directory/'physical_power.json', result)
    return result


def compress_samples(directory):
    source = directory/'field_samples.csv'
    if not source.exists(): return None
    digest, size = r.sha(source), source.stat().st_size
    target = directory/'field_samples.csv.gz'
    with source.open('rb') as fi, gzip.open(target, 'xb', compresslevel=6) as fo: shutil.copyfileobj(fi, fo, 1 << 22)
    with gzip.open(target, 'rb') as check:
        import hashlib
        h = hashlib.sha256()
        for chunk in iter(lambda: check.read(1 << 22), b''): h.update(chunk)
    if h.hexdigest() != digest: raise RuntimeError('Compressed field samples do not round-trip')
    source.unlink()
    return dict(file='field_samples.csv.gz', uncompressed_sha256=digest, uncompressed_bytes=size,
                note='Lossless gzip; decompressed content verified byte-identical before removal of the CSV.')


def wait_for_trace(poll=120):
    """Block until the controlled trace study has ended and no native solver is alive."""
    while True:
        if (TRACE/'final_audit_v2.json').exists() or (TRACE/'stopped_v2.json').exists():
            try:
                preflight(); break
            except RuntimeError:
                pass
        time.sleep(poll)
    final = TRACE/'final_audit_v2.json'
    passed = final.exists() and r.read(final).get('all_original_numerical_checks_passed') is True
    decision = dict(decided_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'), trace_final=str(final) if final.exists() else None,
                    trace_final_sha256=r.sha(final) if final.exists() else None,
                    trace_stopped=(TRACE/'stopped_v2.json').exists(), all_original_numerical_checks_passed=passed)
    if passed:
        decision.update(nt_variant='trace', models=['LQ', 'NT'],
                        note='Controlled trace study passed every original mesh, SUPG and local gate at r=3 6/8 V; '
                             'NT production uses that variant. The observable-scoped release stays on record.')
    else:
        decision.update(nt_variant=None, models=['LQ'],
                        note='Trace study did not pass all original gates; LQ only. NT production waits for the user.')
    OUTPUT.mkdir(exist_ok=True)
    r.save(OUTPUT/f'nt_variant_decision_{time.strftime("%Y%m%dT%H%M%S")}.json', decision)
    return decision


def execute(groups, models=('LQ', 'NT')):
    release = r.read(RELEASE)
    if release.get('release_type') != 'OBSERVABLE_SCOPED_PRODUCTION' or release.get('thresholds_changed') is not False:
        raise RuntimeError('Scoped release record missing or altered')
    gate = r.validate_protocol('paper')
    OUTPUT.mkdir(exist_ok=True)
    if not (OUTPUT/'mesh_inventory.json').exists(): r.save(OUTPUT/'mesh_inventory.json', prepare_meshes())
    identity = dict(release_sha256=r.sha(RELEASE), executable_sha256=r.EXPECTED_EXE, ooc_mb=OOC_MB, threads=r.THREADS,
                    nt_variant=NT_VARIANT, models=list(models),
                    script_hashes={str(p): r.sha(p) for p in [Path(__file__), r.ROOT/'scripts/map_radius_initial_guess_20261021.py',
                                   r.ROOT/'scripts/transfer_hv_initial_guess.py', r.ROOT/'scripts/run_hv_sweep.py',
                                   r.ROOT/'scripts/run_je_reconstruction_mesh_20261018.py',
                                   r.ROOT/'scripts/audit_je_sensitivity_power_20261015.py']})
    stamp = time.strftime('%Y%m%dT%H%M%S')
    r.save(OUTPUT/f'batch_{stamp}_identity.json', dict(**identity, groups=groups, verification_gate=gate))
    events, failed = [], set()
    for job in jobs():
        if job['group'] not in groups or job['model'] not in models: continue
        directory = OUTPUT/job['name']
        if directory.exists() and job['model'] == 'NT' and r.read(directory/'manifest.json').get('nt_variant') != NT_VARIANT:
            failed.add(job['name']); events.append(dict(name=job['name'], state='EXISTING_OTHER_NT_VARIANT_NOT_MIXED')); continue
        if directory.exists():
            events.append(dict(name=job['name'], state='REUSED_'+classify_run(directory)[0] if accepted(directory)
                               else 'EXISTING_INCOMPLETE_NOT_RETRIED'))
            if not accepted(directory): failed.add(job['name'])
            continue
        if job['seed'] and job['seed'][1] in failed:
            failed.add(job['name']); events.append(dict(name=job['name'], state='SKIPPED_SEED_FAILED')); continue
        print('START '+job['name'], flush=True)
        try:
            state = run_native(job, gate)
        except Exception as exc:
            state = 'RUNNER_ERROR: '+repr(exc)
            traceback.print_exc()
        events.append(dict(name=job['name'], state=state))
        print('DONE '+job['name']+' '+state, flush=True)
        if state not in (PASS, ZERO_PASS): failed.add(job['name'])
        r.save(OUTPUT/f'batch_{stamp}_progress_{len(events):03d}.json', dict(events=events))
    r.save(OUTPUT/f'batch_{stamp}_finished.json', dict(events=events, failed=sorted(failed), fully_validated=False))
    return 0 if not failed else 1


def compare(a, b):
    ma, mb = sweep.read_metrics(OUTPUT/a/'metrics.txt'), sweep.read_metrics(OUTPUT/b/'metrics.txt')
    return {k: dict(coarse=ma.get(k), fine=mb.get(k), relative=relative(ma[k], mb[k]) if k in ma and k in mb else None,
                    limit=v, passed=bool(k in ma and k in mb and relative(ma[k], mb[k]) is not None and relative(ma[k], mb[k]) < v))
            for k, v in LIMITS.items()}


def audits():
    out = OUTPUT/'audits'; out.mkdir(exist_ok=True)
    report = {}
    full, half = OUTPUT/f'F_{tag(3, COARSE)}_LQ_V080', OUTPUT/f'H_{tag(3, COARSE)}_LQ_V080'
    if accepted(full) and accepted(half):
        audit_hv_half_y.classify = classify_run   # LQ-aware classification; comparison logic unchanged
        report['lq_full_half_equivalence'] = audit_hv_half_y.audit(full, half)
    levels = [f'H_{tag(3, c)}_LQ_V080' for c in (COARSE, MID, FINEST)]
    if all(accepted(OUTPUT/x) for x in levels):
        report['lq_r3_three_level_8V'] = dict(L0_L1=compare(levels[0], levels[1]), L1_L2=compare(levels[1], levels[2]))
    for a in (2, 4):
        for model in ('LQ', 'NT'):
            pair = [f'H_{tag(a, c)}_{model}_V080' for c in (MID, FINEST)]
            if all(accepted(OUTPUT/x) for x in pair):
                entry = dict(L1_L2=compare(*pair))
                seed = OUTPUT/'seeds'/pair[0]
                if model == 'NT' and seed.exists():
                    entry['local_fields'] = field_shape_audit(seed, OUTPUT/pair[1], OUTPUT/pair[0])
                report[f'{model}_r{a}_two_level_8V'] = entry
    stamp = time.strftime('%Y%m%dT%H%M%S')
    r.save(out/f'audits_{stamp}.json', dict(report, excluded_metrics_note='ring_abs_p_relative_mean is reported but excluded from paper claims'))
    print(json.dumps({k: (v.get('passed') if isinstance(v, dict) and 'passed' in v else 'see file') for k, v in report.items()}, indent=1))


def table():
    rows = []
    for job in jobs():
        d = OUTPUT/job['name']; state, m = classify_run(d)
        row = dict(name=job['name'], group=job['group'], model=job['model'], radius_um=job['radius'],
                   mesh=tag(job['radius'], job['counts']), domain=job['domain'], voltage_V=job['voltage'], state=state)
        factor = 2 if job['domain'] == 'Y_HALF' else 1
        for k in ['joule_W', 'electron_right_A']:
            row['full_'+k] = factor*m[k] if isinstance(m.get(k), float) else None
        for k in ['theta_max_K', 'ring_theta_mean_K', 'ring_E_mean_V_m', 'ring_Q_mean_W_m3', 'ring_Q_p95_W_m3',
                  'E_p95_V_m', 'ring_abs_n_relative_mean', 'minimum_mu_n_ratio', 'stress_von_mises_max_Pa', 'displacement_max_m']:
            row[k] = m.get(k) if isinstance(m.get(k), float) else None
        if (d/'execution.json').exists(): row['seconds'] = r.read(d/'execution.json')['seconds']
        rows.append(row)
    with (OUTPUT/'production_table.csv').open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=['plan', 'meshes', 'execute', 'after-trace', 'audit', 'table'])
    p.add_argument('--groups', nargs='+', default=['lq_check', 'lq_sweep', 'nt_sweep', 'radius', 'radius_mesh'])
    p.add_argument('--nt-variant', choices=['directional', 'trace'], default='directional')
    a = p.parse_args()
    global NT_VARIANT
    NT_VARIANT = a.nt_variant
    if a.command == 'after-trace':
        decision = wait_for_trace()
        if decision['nt_variant']: NT_VARIANT = decision['nt_variant']
        return execute(a.groups, decision['models'])
    if a.command == 'plan': print(json.dumps(jobs(), indent=1)); return 0
    if a.command == 'meshes':
        OUTPUT.mkdir(exist_ok=True); inv = prepare_meshes(); print(json.dumps(inv, indent=1))
        if not (OUTPUT/'mesh_inventory.json').exists(): r.save(OUTPUT/'mesh_inventory.json', inv)
        return 0
    if a.command == 'execute': return execute(a.groups)
    if a.command == 'audit': audits(); return 0
    table(); return 0


if __name__ == '__main__':
    raise SystemExit(main())
