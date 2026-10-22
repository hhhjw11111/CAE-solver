"""Explicit OOC-memory retry after the preserved FULL joint PARDISO -9 failure.

Only the PARDISO in-memory budget changes (6000 -> 8000 MB). Frozen production
code, native inputs, executable, gates, thread count and resource stops remain.
Each batch is explicit, sequential, and stops at the first failed point.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json
import os
import shutil
import subprocess
import sys
import run_je_reconstruction_mesh_20261018 as r

ORIGINAL = r.STUDY/'stage4_paper_mesh'
OUTPUT = r.STUDY/'stage4_paper_mesh_ooc8000_20261020'
FAILED = ORIGINAL/'A40_R20_Z16_S1_V060'
MB = 8000


def preflight():
    command = "@(Get-Process | Where-Object {$_.ProcessName -match '^(TPSolverV2|TPSOLV~[0-9]+)$'} | Select-Object Id,ProcessName) | ConvertTo-Json -Compress"
    result = subprocess.run(['powershell','-NoProfile','-NonInteractive','-Command',command],capture_output=True,check=True)
    if result.stdout.strip(): raise RuntimeError('ANOTHER_NATIVE_SOLVER_ACTIVE: '+result.stdout.decode('utf-8-sig'))


def validate():
    gate = r.validate_protocol('paper')
    ex = r.read(FAILED/'execution.json')
    if ex['exit_code']==0 or 'PARDISO phase 13 failed with error -9' not in (FAILED/'stdout.txt').read_text(encoding='utf-8'):
        raise RuntimeError('Required original OOC -9 evidence absent')
    if ex['input_sha256']!=r.sha(FAILED/'input.nml') or ex['executable_sha256']!=r.EXPECTED_EXE:
        raise RuntimeError('Original failed-attempt identity changed')
    original_identity = r.read(ORIGINAL/'runner_manifest.json')
    for path,digest in original_identity['script_hashes'].items():
        if r.sha(Path(path))!=digest: raise RuntimeError('Frozen dependency changed: '+path)
    for voltage in (6,8):
        for counts in ((48,14,14),(32,22,14)):
            r.accepted(ORIGINAL/r.name(counts,voltage),require_new_completion=True)
    r.accepted(ORIGINAL/r.name((40,14,14),6),require_new_completion=True)
    return gate


def mapping(job, mesh):
    prior = ORIGINAL/'transfers'/job['name']
    if prior.exists():
        report = r.read(prior/'transfer_report.json')
        source = Path(job['source'])
        manifest = r.accepted(source)
        if (report['source_nodes_sha256']!=r.sha(source/'nodes.csv') or
            report['source_mesh_sha256']!=manifest['mesh_sha256'] or
            report['target_mesh_sha256']!=r.sha(mesh) or report['boundary_extension_nodes']!=0):
            raise RuntimeError('Existing transfer identity mismatch')
        return prior
    return r.prepare_mapping(OUTPUT,job,mesh)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=['joint','supg'],required=True)
    parser.add_argument('--joint-level',type=int,choices=[1,2],default=1)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    jobs=r.jobs(args.batch,OUTPUT,args.joint_level)
    plan=dict(batch=args.batch,joint_level=args.joint_level,jobs=jobs,ooc_core_mb=MB,
        threads=r.THREADS,timeout_seconds=r.TIMEOUT,original_failed_run=str(FAILED),
        physics_changed=False,thresholds_changed=False,production_acceptance=False)
    print(json.dumps(plan,indent=2,ensure_ascii=False),flush=True)
    if not args.execute:return
    if os.name!='nt':raise RuntimeError('Windows runtime required')
    os.environ['PYTHONUTF8']='1'
    preflight()
    gate=validate()
    dependencies=[Path(__file__).resolve(),Path(r.__file__).resolve(),r.ROOT/'scripts/audit_je_sensitivity_power_20261015.py']
    identity=dict(scope='JE_PAPER_OOC8000_RESOURCE_RETRY_20261020',stage3_gate=gate,
        script_hashes={str(p):r.sha(p) for p in dependencies},
        original_runner_identity_sha256=r.sha(ORIGINAL/'runner_manifest.json'),
        original_failed_evidence={str(FAILED/p):r.sha(FAILED/p) for p in ('execution.json','input.nml','manifest.json','stdout.txt','resource_usage.json')},
        executable_sha256=r.EXPECTED_EXE,ooc_core_mb=MB,threads=r.THREADS,
        hard_timeout_seconds=r.TIMEOUT,start_disk_bytes=r.START_DISK,stop_disk_bytes=r.STOP_DISK,
        stop_available_ram_bytes=r.MIN_AVAILABLE_RAM,low_ram_consecutive_samples=3,
        rationale_url='https://www.intel.com/content/www/us/en/developer/articles/training/how-to-use-ooc-pardiso.html',
        changed_runtime_option='MKL_PARDISO_OOC_MAX_CORE_SIZE: 6000 -> 8000 MB',
        no_model_or_numerical_threshold_changes=True,original_failure_preserved=True)
    OUTPUT.mkdir(exist_ok=True)
    manifest=OUTPUT/'runner_manifest.json'
    if manifest.exists():
        if r.read(manifest)!=identity:raise RuntimeError('Retry identity changed; fresh output required')
    else:r.save(manifest,identity)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    r.save(OUTPUT/f'batch_{stamp}_plan.json',plan)
    events=[]
    # The imported monitor/solver function is unchanged on disk. This explicit
    # resource setting is recorded by its native manifest and OOC config.
    r.OOC_MB=MB
    try:
        for job in jobs:
            preflight()
            if validate()!=gate:raise RuntimeError('Verification identity changed')
            directory=Path(job['directory'])
            if directory.exists():
                old=r.accepted(directory,require_new_completion=True)
                if old.get('counts')!=job['counts'] or old.get('supg_scale')!=job['supg_scale'] or old.get('voltages')!=[job['voltage']]:
                    raise RuntimeError('Existing retry result mismatch')
                if old['runtime_env']['MKL_PARDISO_OOC_MAX_CORE_SIZE']!=str(MB):raise RuntimeError('OOC budget mismatch')
                events.append(dict(name=job['name'],state='REUSED_COMPLETED_RETRY'));continue
            if shutil.disk_usage(OUTPUT).free<r.START_DISK:raise RuntimeError('RESOURCE_STOP_BEFORE_PREPARATION')
            oldmesh=ORIGINAL/'meshes'/r.name(job['counts'],0).split('_S')[0]
            mesh=r.ensure_mesh(ORIGINAL if oldmesh.exists() else OUTPUT,job['counts'])
            transfer=mapping(job,mesh)
            r.run_native(job,mesh,transfer/'nodes.csv',directory,gate)
            subprocess.run([sys.executable,'-X','utf8',str(r.ROOT/'scripts/audit_je_sensitivity_power_20261015.py'),
                '--runs',str(directory)],check=True)
            comparison=r.field_comparison(transfer,Path(job['source']),directory)
            r.save(directory/'mapped_field_comparison.json',comparison)
            r.save(directory/'completion.json',dict(power_audit_completed=True,field_comparison_passed=comparison['passed'],
                fully_validated=False,source_run=job['source'],mapping_directory=str(transfer),
                geometry_and_transfer_only_initial_guess=True,resource_retry_only=True))
            if job['name']==FAILED.name and r.sha(directory/'input.nml')!=r.sha(FAILED/'input.nml'):
                raise RuntimeError('Resource-only retry did not preserve the exact native input')
            events.append(dict(name=job['name'],state=r.PASS,local_comparison_passed=comparison['passed']))
            print('COMPLETED '+job['name']+'; local comparison passed='+str(comparison['passed']),flush=True)
        r.save(OUTPUT/f'batch_{stamp}_audits.json',r.batch_audits(args.batch,OUTPUT,args.joint_level))
        r.save(OUTPUT/f'batch_{stamp}_finished.json',dict(all_requested_attempts_completed=True,events=events,fully_validated=False))
    except BaseException as exc:
        r.save(OUTPUT/f'batch_{stamp}_stopped.json',dict(reason=repr(exc),events=events,fully_validated=False,
            next_action='Review preserved failure; no automatic next job or budget increase'))
        raise


if __name__=='__main__':main()
