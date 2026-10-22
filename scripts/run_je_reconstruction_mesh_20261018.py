"""Review or execute explicit, sequential batches of the frozen JE mesh study.

Default is read-only plan output. --execute requires the selected verification
scope's PASS gate. Original and paper scopes never share output roots.
The first executable batch is probe (FULL A40/R14/Z14, SUPG=1, 6 V).
No batch certifies the whole study; production remains separately gated.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import ctypes
import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
from run_je_sensitivity_resource_20261016 import EXE, BASE, short, replace
from run_hv_sweep import sha, classify, audit_nodes
from run_hv_regression import runtime_env
from hv_runtime import configure_storage, recorded_runtime
from generate_je_reconstruction_mesh_20261018 import write as write_mesh
from audit_hv_field_shape import audit as field_comparison
from audit_hv_mesh import audit as mesh_comparison
from audit_hv_stabilization import audit as supg_comparison

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT/'runs/JE_RECONSTRUCTION_20261018'
OLD = ROOT/'runs/JE_SENSITIVITY_20261015'
DIRECTIONAL = ROOT/'runs/JE_FINE_CROSS_20261016/GLOBAL_DIRECTIONAL'
EXPECTED_EXE = '405e953c8cb0b1eca2a61eccc011945892924b43fca6805e87fd8fec19006d92'
PASS = 'SOLVE_FLUX_PASS_MESH_PENDING'
THREADS, OOC_MB, TIMEOUT = 8, 6000, 7200
GIB = 1024**3
START_DISK = 16*GIB
STOP_DISK = 4*GIB
MIN_AVAILABLE_RAM = 512*1024**2
PAPER_SCOPE = 'JE_PAPER_MINIMAL_20261018'
PAPER_REQUIRED_IDS = {'A1','A2','A1JJE','JEGRAD','G02JE','F1MJE','G04JE','SUPG','SUPGJ','SUPGJT'}


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def save(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)


def name(counts, voltage, scale=1.):
    a, r, z = counts
    return f'A{a}_R{r}_Z{z}_S{scale:g}_V{voltage*10:03d}'


def baseline(voltage):
    return OLD/f'A32_S1_V{voltage*10:03d}'


def radial18(voltage):
    return DIRECTIONAL/name((32,18,14),voltage)


def jobs(batch, output, joint_level):
    result = []
    def add(counts, v, source, kind, scale=1.):
        result.append(dict(name=name(counts,v,scale), counts=list(counts), voltage=v, supg_scale=scale,
                           source=str(source), directory=str(output/name(counts,v,scale)), kind=kind))
    if batch == 'probe':
        add((40,14,14),6,baseline(6),'angular')
    elif batch == 'angular':
        add((40,14,14),8,baseline(8),'angular')
        for v in (6,8): add((48,14,14),v,output/name((40,14,14),v),'angular')
    elif batch == 'radial':
        for r in (20,22):
            for v in (6,8):
                source = radial18(v) if r == 20 else output/name((32,20,14),v)
                add((32,r,14),v,source,'radial')
    elif batch == 'joint':
        counts = (40,20,16) if joint_level == 1 else (48,22,18)
        for v in (6,8):
            source = radial18(v) if joint_level == 1 else output/name((40,20,16),v)
            add(counts,v,source,'joint')
    elif batch == 'supg':
        counts = (40,20,16) if joint_level == 1 else (48,22,18)
        for scale in (.5,2.):
            for v in (6,8): add(counts,v,output/name(counts,v),'SUPG',scale)
    return result


def verification_gate(verification_scope):
    if verification_scope == 'original':
        return STUDY/'stage3_verification/acceptance.json'
    if verification_scope == 'paper':
        return STUDY/'paper_scope/verification_reuse/acceptance.json'
    raise ValueError('Unknown verification scope: '+str(verification_scope))


def validate_protocol(verification_scope='original'):
    contract = read(STUDY/'contract.json')
    gate_path = verification_gate(verification_scope)
    if not gate_path.is_file(): raise RuntimeError('STAGE3_GATE_MISSING: '+str(gate_path))
    gate_sha = sha(gate_path)
    gate = read(gate_path)
    reuse_evidence = None
    if verification_scope == 'original':
        # Preserve the original full-suite gate. The excluded A2P failure is not
        # converted into PASS and the paper gate is never an implicit fallback.
        if (gate.get('all_passed') is not True or gate.get('core_sources_unchanged') is not True
                or gate.get('build_manifests_verified') is not True):
            raise RuntimeError('STAGE3_GATE_NOT_PASSED: no mesh solver may start')
    else:
        if (gate.get('scope') != PAPER_SCOPE or gate.get('all_required_passed') is not True
                or gate.get('core_sources_unchanged') is not True or gate.get('build_manifests_verified') is not True):
            raise RuntimeError('PAPER_VERIFICATION_GATE_NOT_PASSED')
        if (set(gate.get('required_ids',[])) != PAPER_REQUIRED_IDS or len(gate.get('required_ids',[])) != 10
                or gate.get('targeted_passed') != 10 or gate.get('targeted_total') != 10
                or gate.get('legacy_passed') != 16 or gate.get('legacy_total') != 16):
            raise RuntimeError('PAPER_REQUIRED_CHECKS_INCOMPLETE: require ten targeted checks and sixteen legacy CLI comparisons')
        # This rehashes original logs, sources and both build manifests and
        # reconstructs every required result. A stored boolean is insufficient.
        from audit_je_paper_reuse_20261018 import validate_acceptance
        revalidated = validate_acceptance(gate_path,expected_exe=EXE)
        if (revalidated.get('scope') != PAPER_SCOPE or revalidated.get('all_required_passed') is not True
                or set(revalidated.get('required_ids',[])) != PAPER_REQUIRED_IDS
                or revalidated.get('targeted_passed') != 10 or revalidated.get('legacy_passed') != 16):
            raise RuntimeError('PAPER_RECOMPUTED_VERIFICATION_FAILED')
        reuse_evidence = dict(auditor_sha256=sha(ROOT/'scripts/audit_je_paper_reuse_20261018.py'),
                              evidence_file_hashes=revalidated['evidence_file_hashes'],
                              targeted_passed=10,targeted_total=10,legacy_passed=16,legacy_total=16)
    if sha(gate_path) != gate_sha: raise RuntimeError('Verification acceptance changed during validation')
    if sha(EXE) != EXPECTED_EXE: raise RuntimeError('Historical JE executable identity changed')
    if sha(ROOT/contract['protocol']) != contract['protocol_sha256']:
        raise RuntimeError('Frozen migration protocol changed')
    if len(contract['core_sources_before']) != 40:
        raise RuntimeError('Frozen core inventory must contain all forty production source files')
    for item in contract['core_sources_before']:
        if sha(ROOT.parent/item['path']) != item['sha256']:
            raise RuntimeError('Production source changed: '+item['path'])
    for path, digest in contract['frozen_gate_scripts'].items():
        if sha(ROOT/path) != digest: raise RuntimeError('Frozen acceptance script changed: '+path)
    return dict(stage3_acceptance=str(gate_path),stage3_sha256=gate_sha,
                verification_scope=verification_scope,scope=PAPER_SCOPE if verification_scope=='paper' else 'JE_RECONSTRUCTION_ORIGINAL_20261018',
                recomputed_paper_evidence=reuse_evidence,
                contract_sha256=sha(STUDY/'contract.json'),executable_sha256=EXPECTED_EXE)


def accepted(directory, require_new_completion=False):
    directory = Path(directory)
    execution, manifest = read(directory/'execution.json'), read(directory/'manifest.json')
    state, metrics = classify(directory,execution['exit_code'])
    if state != PASS or read(directory/'field_audit.json').get('passed') is not True:
        raise RuntimeError('Seed/result lacks solve, flux or field acceptance: '+str(directory))
    if execution['input_sha256'] != sha(directory/'input.nml'):
        raise RuntimeError('Changed input: '+str(directory))
    if manifest['executable_sha256'] != EXPECTED_EXE or execution['executable_sha256'] != EXPECTED_EXE:
        raise RuntimeError('Seed/result executable differs from frozen study')
    if sha(Path(manifest['mesh'])) != manifest['mesh_sha256']:
        raise RuntimeError('Changed seed mesh')
    if manifest.get('domain_symmetry','FULL') != 'FULL' or metrics.get('domain_symmetry','FULL') != 'FULL':
        raise RuntimeError('This entry implements FULL only')
    if metrics.get('heat_source') != 'E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE':
        raise RuntimeError('Seed is not the JE branch')
    if require_new_completion and read(directory/'completion.json').get('power_audit_completed') is not True:
        raise RuntimeError('Previous attempt incomplete; preserve it and use a new output root')
    return manifest


def no_native_solver():
    # Read-only inventory; never stop another task's process.
    command = "@(Get-Process -Name TPSolverV2 -ErrorAction SilentlyContinue | Select-Object Id,ProcessName) | ConvertTo-Json -Compress"
    proc = subprocess.run(['powershell','-NoProfile','-NonInteractive','-Command',command],
                          capture_output=True,check=True)
    if proc.stdout.strip(): raise RuntimeError('ANOTHER_NATIVE_SOLVER_ACTIVE: '+proc.stdout.decode('utf-8-sig').strip())


class WindowsMonitor:
    def __init__(self,pid):
        from ctypes import wintypes
        class Counters(ctypes.Structure):
            _fields_ = [('cb',wintypes.DWORD),('PageFaultCount',wintypes.DWORD)] + [
                (key,ctypes.c_size_t) for key in ('PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage',
                'QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage',
                'PeakPagefileUsage','PrivateUsage')]
        class Memory(ctypes.Structure):
            _fields_ = [('dwLength',wintypes.DWORD),('dwMemoryLoad',wintypes.DWORD)] + [
                (key,ctypes.c_ulonglong) for key in ('ullTotalPhys','ullAvailPhys','ullTotalPageFile',
                'ullAvailPageFile','ullTotalVirtual','ullAvailVirtual','ullAvailExtendedVirtual')]
        self.Counters,self.Memory = Counters,Memory
        self.kernel = ctypes.WinDLL('kernel32',use_last_error=True)
        self.psapi = ctypes.WinDLL('psapi',use_last_error=True)
        self.kernel.OpenProcess.argtypes = [wintypes.DWORD,wintypes.BOOL,wintypes.DWORD]
        self.kernel.OpenProcess.restype = wintypes.HANDLE
        self.kernel.CloseHandle.argtypes = [wintypes.HANDLE]
        self.psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE,ctypes.POINTER(Counters),wintypes.DWORD]
        self.kernel.GlobalMemoryStatusEx.argtypes = [ctypes.POINTER(Memory)]
        self.handle = self.kernel.OpenProcess(0x0400|0x0010,False,pid)
        if not self.handle: raise ctypes.WinError(ctypes.get_last_error())

    def sample(self):
        counters,memory = self.Counters(),self.Memory()
        counters.cb,memory.dwLength = ctypes.sizeof(counters),ctypes.sizeof(memory)
        if not self.psapi.GetProcessMemoryInfo(self.handle,ctypes.byref(counters),counters.cb):
            raise ctypes.WinError(ctypes.get_last_error())
        if not self.kernel.GlobalMemoryStatusEx(ctypes.byref(memory)):
            raise ctypes.WinError(ctypes.get_last_error())
        return dict(rss_bytes=counters.WorkingSetSize,peak_rss_bytes=counters.PeakWorkingSetSize,
                    private_bytes=counters.PrivateUsage,peak_pagefile_bytes=counters.PeakPagefileUsage,
                    available_physical_bytes=memory.ullAvailPhys)

    def close(self):
        self.kernel.CloseHandle(self.handle)


def scratch_bytes(scratch):
    total = 0
    for path in scratch.iterdir():
        try:
            if path.is_file(): total += path.stat().st_size
        except FileNotFoundError:
            # PARDISO may replace a factor between directory enumeration and stat.
            continue
    return total


def cleanup_owned_factors(directory,process):
    if process.poll() is None: raise RuntimeError('Cannot clean factors while this native process is alive')
    directory = directory.resolve()
    scratch = (directory/'pardiso_ooc').resolve(strict=True)
    if scratch.parent != directory or not directory.is_relative_to(STUDY.resolve()):
        raise RuntimeError('Unsafe scratch path; no factors removed')
    selected = []
    for path in scratch.iterdir():
        if not re.fullmatch(r'factors(?:[0-9_.-].*)?',path.name): continue
        if path.is_symlink() or not path.is_file() or path.resolve().parent != scratch:
            raise RuntimeError('Unexpected factor target; no recursive cleanup is allowed')
        stat = path.stat()
        selected.append(dict(path=str(path.resolve()),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns))
    save(directory/'temporary_factor_cleanup_plan.json',dict(pid=process.pid,exit_code=process.returncode,
         scope='Only this fresh run directory; original scientific results and run_ooc.cfg retained',files=selected))
    removed = []
    for entry in selected:
        path = Path(entry['path'])
        stat = path.stat()
        if path.resolve().parent != scratch or stat.st_size != entry['bytes'] or stat.st_mtime_ns != entry['mtime_ns']:
            raise RuntimeError('Factor changed after cleanup inventory')
        path.unlink()
        removed.append(entry)
    save(directory/'temporary_factor_cleanup.json',dict(completed=True,files=removed,
         removed_bytes=sum(x['bytes'] for x in removed),native_exit_code=process.returncode))


def run_native(job,mesh,seed,directory,gate):
    no_native_solver()
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
    inp = directory/'input.nml'
    inp.write_text(text,encoding='ascii')
    env = runtime_env()
    env.update(OMP_NUM_THREADS=str(THREADS),MKL_NUM_THREADS=str(THREADS))
    env = configure_storage(env,Path(short(directory)),'OOC',OOC_MB,False)
    manifest = dict(executable=str(EXE),executable_sha256=EXPECTED_EXE,mesh=str(mesh),mesh_sha256=sha(mesh),
        model='JE_CURRENT_WORK_DIAGNOSTIC_ZNO_NT',model_contract='JE_SIMPLIFIED_R0_V1',voltages=[voltage],
        supg_scale=scale,domain_symmetry='FULL',log_weighted_quadrature_order=3,counts=job['counts'],
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
        fields = audit_nodes(directory/'nodes.csv',mesh,voltage)
        save(directory/'field_audit.json',fields)
        if not fields['passed']: state = 'FAILED_FIELD_AUDIT'
    if reason: state = 'RESOURCE_OR_TIMEOUT_STOP'
    save(directory/'summary.json',dict(state=state,mesh_validated=False,power_gap_validated=False,
                                       elapsed_seconds=time.monotonic()-started,stop_reason=reason))
    if state != PASS: raise RuntimeError('NATIVE_BATCH_STOP: '+job['name']+' '+state)
    if metrics.get('heat_source') != 'E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE':
        raise RuntimeError('Wrong solved heat model')


def ensure_mesh(output,counts):
    folder = output/'meshes'/name(counts,0).split('_S')[0]
    if folder.exists():
        manifest = read(folder/'manifest.json')
        mesh = folder/'H1.mesh'
        if list(manifest['counts'].values()) != counts or sha(mesh) != manifest['mesh_sha256']:
            raise RuntimeError('Existing new mesh identity mismatch')
        return mesh
    mesh,_ = write_mesh(folder,*counts)
    return mesh


def prepare_mapping(output,job,mesh):
    source = Path(job['source'])
    manifest = accepted(source,require_new_completion=source.is_relative_to(output))
    _,metrics = classify(source,read(source/'execution.json')['exit_code'])
    if metrics['voltage_V'] != job['voltage']: raise RuntimeError('Seed voltage mismatch')
    mapping = output/'transfers'/job['name']
    if not mapping.exists():
        subprocess.run([sys.executable,str(ROOT/'scripts/transfer_hv_initial_guess.py'),
                        '--source-mesh',manifest['mesh'],'--source-nodes',str(source/'nodes.csv'),
                        '--target-mesh',str(mesh),'--output',str(mapping)],check=True)
    report = read(mapping/'transfer_report.json')
    if report['source_nodes_sha256'] != sha(source/'nodes.csv') or report['target_mesh_sha256'] != sha(mesh):
        raise RuntimeError('Mapping identity mismatch')
    if report['boundary_extension_nodes'] != 0:
        raise RuntimeError('Mapping used boundary extrapolation; preserve and review before any solve')
    return mapping


def batch_audits(batch,output,joint_level):
    result = dict(batch=batch,fully_validated=False,studies=[])
    for v in (6,8):
        if batch == 'angular':
            levels = [baseline(v)]+[output/name((a,14,14),v) for a in (40,48)]
        elif batch == 'radial':
            levels = [radial18(v)]+[output/name((32,r,14),v) for r in (20,22)]
        elif batch == 'joint' and joint_level == 2:
            levels = [radial18(v),output/name((40,20,16),v),output/name((48,22,18),v)]
        elif batch == 'supg':
            counts = (40,20,16) if joint_level == 1 else (48,22,18)
            result['studies'].append(supg_comparison([output/name(counts,v,s) for s in (.5,1.,2.)]))
            continue
        else: continue
        result['studies'].append(mesh_comparison(levels,[v]))
    result['limitation'] = 'All original gates retained; local comparisons are saved per run. Independent directions and one joint candidate are not full spatial/SUPG certification.'
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=['probe','angular','radial','joint','supg'],required=True)
    parser.add_argument('--joint-level',type=int,choices=[1,2],default=1,help='1=A40/R20/Z16, 2=A48/R22/Z18; explicit joint or SUPG batch only')
    parser.add_argument('--verification-scope',choices=['original','paper'],default='original',
                        help='Explicit original full-suite gate or user-approved minimal paper verification; no fallback')
    parser.add_argument('--output-root',type=Path,help='Default: stage4_mesh for original; stage4_paper_mesh for paper')
    parser.add_argument('--execute',action='store_true',help='Run only after stage3 PASS; omission prints a read-only plan')
    args = parser.parse_args()
    default_name = 'stage4_paper_mesh' if args.verification_scope == 'paper' else 'stage4_mesh'
    output = (args.output_root or STUDY/default_name).resolve()
    if not output.is_relative_to(STUDY.resolve()) or output == STUDY.resolve():
        raise ValueError('Output must be a new study subdirectory beneath JE_RECONSTRUCTION_20261018')
    is_paper_root = output.parent == STUDY.resolve() and (output.name=='stage4_paper_mesh' or output.name.startswith('stage4_paper_mesh_'))
    if args.verification_scope == 'paper' and not is_paper_root:
        raise ValueError('Paper scope requires its separate stage4_paper_mesh[_suffix] output root')
    if args.verification_scope == 'original' and any(p.name.startswith('stage4_paper_mesh') for p in [output,*output.parents]):
        raise ValueError('Original scope cannot use a paper output root')
    queue = jobs(args.batch,output,args.joint_level)
    plan = dict(batch=args.batch,jobs=queue,domain='FULL',threads=THREADS,ooc_core_mb=OOC_MB,hard_timeout_seconds=TIMEOUT,
                verification_scope=args.verification_scope,
                scope=PAPER_SCOPE if args.verification_scope=='paper' else 'JE_RECONSTRUCTION_ORIGINAL_20261018',
                stage3_gate=str(verification_gate(args.verification_scope)),executable=str(EXE),
                requires_existing_probe=args.batch!='probe',fully_validated=False)
    print(json.dumps(plan,indent=2,ensure_ascii=False),flush=True)
    if not args.execute: return 0
    if os.name != 'nt': raise RuntimeError('This executor requires Windows Intel runtime')
    gate = validate_protocol(args.verification_scope)
    no_native_solver()
    if args.batch != 'probe': accepted(output/name((40,14,14),6),require_new_completion=True)
    if args.batch in ('joint','supg'):
        # Joint batches follow completed independent studies, without promoting
        # any failed mesh comparison to PASS. Scientific review remains explicit.
        for v in (6,8):
            for counts in ((48,14,14),(32,22,14)):
                accepted(output/name(counts,v),require_new_completion=True)
    output.mkdir(parents=True,exist_ok=True)
    root_manifest = output/'runner_manifest.json'
    scripts = [Path(__file__),ROOT/'scripts/generate_je_reconstruction_mesh_20261018.py',
               ROOT/'scripts/transfer_hv_initial_guess.py',ROOT/'scripts/run_je_sensitivity_resource_20261016.py',
               ROOT/'scripts/hv_runtime.py',ROOT/'scripts/run_hv_sweep.py',ROOT/'scripts/diagnose_je_power_gap_20261015.py']
    if args.verification_scope == 'paper': scripts.append(ROOT/'scripts/audit_je_paper_reuse_20261018.py')
    identity = dict(script_hashes={str(p):sha(p) for p in scripts},executable_sha256=EXPECTED_EXE,
                    contract_sha256=sha(STUDY/'contract.json'),domain_symmetry='FULL',
                    verification_scope=args.verification_scope,verification_evidence=gate)
    if root_manifest.exists():
        if read(root_manifest) != identity: raise RuntimeError('Runner identity changed; use a fresh output root')
    else: save(root_manifest,identity)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    save(output/f'batch_{stamp}_plan.json',dict(**plan,stage3_gate_evidence=gate))
    events = []
    try:
        for job in queue:
            if validate_protocol(args.verification_scope) != gate:
                raise RuntimeError('Verification evidence changed within the batch')
            directory = Path(job['directory'])
            if directory.exists():
                prior = accepted(directory,require_new_completion=True)
                if prior.get('counts') != job['counts'] or prior.get('supg_scale') != job['supg_scale'] or prior.get('voltages') != [job['voltage']]:
                    raise RuntimeError('Existing result does not match the requested batch point')
                events.append(dict(name=job['name'],state='EXPLICITLY_REUSED_COMPLETED_SAME_ROOT_RESULT'))
                continue
            if shutil.disk_usage(output).free < START_DISK: raise RuntimeError('RESOURCE_STOP_BEFORE_PREPARATION')
            mesh = ensure_mesh(output,job['counts'])
            mapping = prepare_mapping(output,job,mesh)
            run_native(job,mesh,mapping/'nodes.csv',directory,gate)
            subprocess.run([sys.executable,str(ROOT/'scripts/audit_je_sensitivity_power_20261015.py'),
                            '--runs',str(directory)],check=True)
            comparison = field_comparison(mapping,Path(job['source']),directory)
            save(directory/'mapped_field_comparison.json',comparison)
            save(directory/'completion.json',dict(power_audit_completed=True,field_comparison_passed=comparison['passed'],
                 fully_validated=False,source_run=job['source'],geometry_and_transfer_only_initial_guess=True))
            events.append(dict(name=job['name'],state=PASS,local_comparison_passed=comparison['passed']))
            print('COMPLETED '+job['name']+'; local comparison passed='+str(comparison['passed']),flush=True)
        report = batch_audits(args.batch,output,args.joint_level)
        save(output/f'batch_{stamp}_audits.json',report)
        save(output/f'batch_{stamp}_finished.json',dict(all_requested_attempts_completed=True,events=events,fully_validated=False))
    except BaseException as exc:
        save(output/f'batch_{stamp}_stopped.json',dict(reason=repr(exc),events=events,fully_validated=False,
             no_failed_run_is_skipped=True,next_action='Review preserved evidence; no automatic retry or next native job'))
        raise
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)
