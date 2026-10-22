"""JE half-domain execution adapter; frozen native model and original gates.

The resource monitor is copied from run_je_reconstruction_mesh_20261018.run_native.
Changes: robust process preflight, explicit half_y_symmetry input/manifest, and
half-domain nodal auditing. OOC remains 6000 MB; native code is untouched.
"""
from pathlib import Path
from run_je_reconstruction_mesh_20261018 import (
    shutil, time, subprocess, json, BASE, ROOT, EXE, EXPECTED_EXE,
    THREADS, OOC_MB, TIMEOUT, START_DISK, STOP_DISK, MIN_AVAILABLE_RAM,
    short, replace, runtime_env, configure_storage, recorded_runtime,
    sha, save, WindowsMonitor, scratch_bytes, cleanup_owned_factors,
    classify, audit_nodes, PASS)
from run_je_paper_ooc8000_20261020 import preflight

def run_native(job,mesh,seed,directory,gate):
    preflight()
    if shutil.disk_usage(directory.parent).free < START_DISK:
        raise RuntimeError('RESOURCE_STOP: fewer than 16 GiB available before the next solve')
    directory.mkdir(exist_ok=False)
    voltage,scale = job['voltage'],job['supg_scale']
    text = (BASE/f'input_V{voltage*10:03d}.nml').read_text(encoding='utf-8-sig')
    changes = dict(case_id="'JE_RECON_"+job['name']+"'",mesh_file="'"+short(mesh)+"'",
                   restart_nodes_file="'"+short(seed)+"'",restart_voltage_factor='1.0',supg_scale=str(scale))
    for key,value in changes.items(): text = replace(text,key,value)
    for literal in ('heat_source_current_work=.true.','residual_tolerance=1e-9,update_tolerance=1e-9',
                    'target_joule_factor=1.0','restart_joule_factor=1.0'):
        if literal not in text: raise RuntimeError('Frozen input requirement missing: '+literal)
    if text.count('&problem_control') != 1 or 'half_y_symmetry' in text:
        raise RuntimeError('Unexpected frozen half-domain template')
    text = text.replace('&problem_control', '&problem_control\n half_y_symmetry=.true.,')
    changes['half_y_symmetry'] = '.true.'
    inp = directory/'input.nml'
    inp.write_text(text,encoding='ascii')
    env = runtime_env()
    env.update(OMP_NUM_THREADS=str(THREADS),MKL_NUM_THREADS=str(THREADS))
    env = configure_storage(env,Path(short(directory)),'OOC',OOC_MB,False)
    manifest = dict(executable=str(EXE),executable_sha256=EXPECTED_EXE,mesh=str(mesh),mesh_sha256=sha(mesh),
        model='JE_CURRENT_WORK_DIAGNOSTIC_ZNO_NT',model_contract='JE_SIMPLIFIED_R0_V1',voltages=[voltage],
        supg_scale=scale,domain_symmetry='Y_HALF',log_weighted_quadrature_order=3,counts=job['counts'],
        runtime_env=recorded_runtime(env),seed_nodes=str(seed),seed_sha256=sha(seed),
        initial_guess_status='REQUIRES_REEQUILIBRATION',heat_source_unchanged=True,input_template=str(BASE/f'input_V{voltage*10:03d}.nml'),
        changed_input_keys=list(changes),wall_time_limit_seconds=TIMEOUT,stage3_gate=gate,
        resource_limits=dict(start_disk_bytes=START_DISK,runtime_stop_disk_bytes=STOP_DISK,
                             stop_available_ram_bytes=MIN_AVAILABLE_RAM,low_ram_consecutive_samples=3))
    save(directory/'manifest.json',manifest)
    started = time.monotonic()
    reason,process,monitor = None,None,None
    peak = dict(peak_rss_bytes=0,peak_private_bytes=0,peak_pagefile_bytes=0,peak_scratch_bytes=0,
                minimum_disk_free_bytes=shutil.disk_usage(directory).free,minimum_available_physical_bytes=None)
    try:
        with (directory/'stdout.txt').open('xb') as log, (directory/'resource_samples.jsonl').open('x',encoding='ascii') as samples:
            process = subprocess.Popen([short(EXE),'run-f1',short(inp),short(directory)],cwd=short(ROOT),env=env,
                                       stdout=log,stderr=subprocess.STDOUT)
            save(directory/'process.json',dict(pid=process.pid,executable=str(EXE),owned_by_this_run=True))
            monitor = WindowsMonitor(process.pid)
            low_ram = 0
            while process.poll() is None:
                try:
                    sample = monitor.sample()
                except OSError:
                    if process.poll() is not None: break
                    raise
                sample.update(elapsed_seconds=time.monotonic()-started,
                              disk_free_bytes=shutil.disk_usage(directory).free,
                              scratch_bytes=scratch_bytes(directory/'pardiso_ooc'))
                samples.write(json.dumps(sample)+'\n');samples.flush()
                for key,source in [('peak_rss_bytes','peak_rss_bytes'),('peak_private_bytes','private_bytes'),
                                   ('peak_pagefile_bytes','peak_pagefile_bytes'),('peak_scratch_bytes','scratch_bytes')]:
                    peak[key] = max(peak[key],sample[source])
                peak['minimum_disk_free_bytes'] = min(peak['minimum_disk_free_bytes'],sample['disk_free_bytes'])
                previous = peak['minimum_available_physical_bytes']
                peak['minimum_available_physical_bytes'] = min(previous if previous is not None else sample['available_physical_bytes'],sample['available_physical_bytes'])
                low_ram = low_ram+1 if sample['available_physical_bytes'] < MIN_AVAILABLE_RAM else 0
                if sample['elapsed_seconds'] >= TIMEOUT: reason = 'HARD_TIMEOUT_7200_SECONDS'
                elif sample['disk_free_bytes'] < STOP_DISK: reason = 'RESOURCE_STOP_DISK_BELOW_4_GIB'
                elif low_ram >= 3: reason = 'RESOURCE_STOP_AVAILABLE_RAM_BELOW_512_MIB'
                if reason:
                    print(reason+'; stopping only owned PID '+str(process.pid),flush=True)
                    process.kill();break
                time.sleep(1)
            process.wait()
    except BaseException as exc:
        reason = reason or 'MONITOR_OR_RUNNER_EXCEPTION: '+repr(exc)
        if process is not None and process.poll() is None:
            process.kill();process.wait()
        raise
    finally:
        if monitor is not None: monitor.close()
        if process is not None:
            peak['peak_scratch_bytes'] = max(peak['peak_scratch_bytes'],scratch_bytes(directory/'pardiso_ooc'))
            save(directory/'execution.json',dict(exit_code=process.returncode,seconds=time.monotonic()-started,
                 timed_out=reason=='HARD_TIMEOUT_7200_SECONDS',stop_reason=reason,input_sha256=sha(inp),executable_sha256=EXPECTED_EXE))
            save(directory/'resource_usage.json',dict(**peak,sampling_seconds=1,scope='Native process only; private peak sampled; RSS/pagefile peaks include Windows counters'))
            cleanup_owned_factors(directory,process)
    state,metrics = classify(directory,process.returncode)
    if (directory/'nodes.csv').exists():
        fields = audit_nodes(directory/'nodes.csv',mesh,voltage,half_y=True)
        save(directory/'field_audit.json',fields)
        if not fields['passed']: state = 'FAILED_FIELD_AUDIT'
    if reason: state = 'RESOURCE_OR_TIMEOUT_STOP'
    save(directory/'summary.json',dict(state=state,mesh_validated=False,power_gap_validated=False,
                                       elapsed_seconds=time.monotonic()-started,stop_reason=reason))
    if state != PASS: raise RuntimeError('NATIVE_BATCH_STOP: '+job['name']+' '+state)
    if metrics.get('heat_source') != 'E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE':
        raise RuntimeError('Wrong solved heat model')

