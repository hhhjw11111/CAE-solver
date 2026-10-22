"""Resume interrupted paper batches in a fresh root, reusing complete states by reference.

Run this as an independent hidden background process. No completed data are
copied, renamed, or rerun. Missing exit status never becomes a numerical FAIL.
The frozen 8000 MB resource policy and all original scientific gates remain.
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
import run_je_paper_ooc8000_20261020 as policy

PRIOR = policy.OUTPUT
OUTPUT = r.STUDY/'stage4_paper_mesh_resumed_20261020'
INTERRUPTED = PRIOR/'A40_R20_Z16_S1_V080'


def resolve_complete(name):
    for root in (OUTPUT,PRIOR,policy.ORIGINAL):
        directory=root/name
        if (directory/'completion.json').exists():
            r.accepted(directory,require_new_completion=True)
            return directory
    raise RuntimeError('No completed same-point source: '+name)


def evidence():
    gate=policy.validate()
    prior_identity=r.read(PRIOR/'runner_manifest.json')
    for path,digest in prior_identity['script_hashes'].items():
        if r.sha(Path(path))!=digest:raise RuntimeError('Prior retry code changed: '+path)
    for filename in ('execution.json','completion.json'):
        if (INTERRUPTED/filename).exists():
            raise RuntimeError('Interruption premise changed: '+filename)
    owner=r.read(INTERRUPTED/'process.json')
    if owner.get('owned_by_this_run') is not True:raise RuntimeError('Missing process ownership')
    frozen={str(INTERRUPTED/p):r.sha(INTERRUPTED/p) for p in
        ('input.nml','manifest.json','process.json','stdout.txt','metrics.txt','resource_samples.jsonl')}
    completed=resolve_complete('A40_R20_Z16_S1_V060')
    return dict(gate=gate,prior_runner_identity_sha256=r.sha(PRIOR/'runner_manifest.json'),
        interrupted_evidence=frozen,interrupted_pid=owner['pid'],
        interruption_status='PROCESS_ABSENT_NO_EXIT_OR_COMPLETION_RECORD_CAUSE_UNCONFIRMED',
        completed_joint1_6V=str(completed),
        completed_joint1_6V_completion_sha256=r.sha(completed/'completion.json'))


def queue(batch,level):
    jobs=r.jobs(batch,OUTPUT,level)
    for job in jobs:
        if batch=='joint' and level==2:
            job['source']=str(resolve_complete(r.name((40,20,16),job['voltage'])))
        elif batch=='supg':
            job['source']=str(resolve_complete(r.name(job['counts'],job['voltage'])))
    return jobs


def transfer(job,mesh):
    source=Path(job['source'])
    source_manifest=r.accepted(source,require_new_completion=source.parent in (PRIOR,OUTPUT))
    for root in (OUTPUT,PRIOR,policy.ORIGINAL):
        directory=root/'transfers'/job['name']
        if not directory.exists():continue
        report=r.read(directory/'transfer_report.json')
        if (report['source_nodes_sha256']==r.sha(source/'nodes.csv') and
            report['source_mesh_sha256']==source_manifest['mesh_sha256'] and
            report['target_mesh_sha256']==r.sha(mesh) and report['boundary_extension_nodes']==0):
            return directory
        if root==OUTPUT:raise RuntimeError('New-root transfer identity mismatch')
    return r.prepare_mapping(OUTPUT,job,mesh)


def batch_audit(batch,level):
    studies=[]
    for voltage in (6,8):
        if batch=='joint' and level==2:
            studies.append(r.mesh_comparison([r.radial18(voltage),
                resolve_complete(r.name((40,20,16),voltage)),
                resolve_complete(r.name((48,22,18),voltage))],[voltage]))
        elif batch=='supg':
            counts=(40,20,16) if level==1 else (48,22,18)
            studies.append(r.supg_comparison([resolve_complete(r.name(counts,voltage,s)) for s in (.5,1.,2.)]))
    return dict(batch=batch,joint_level=level,studies=studies,fully_validated=False,
        limitation='Explicit references to actual completed runs across immutable roots. No production acceptance.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=['joint','supg'],required=True)
    parser.add_argument('--joint-level',choices=[1,2],type=int,default=1)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    jobs=queue(args.batch,args.joint_level)
    plan=dict(batch=args.batch,joint_level=args.joint_level,jobs=jobs,ooc_core_mb=8000,threads=8,
        timeout_seconds=7200,completed_points_reused_by_reference=True,fully_validated=False)
    print(json.dumps(plan,indent=2,ensure_ascii=False),flush=True)
    if not args.execute:return
    os.environ['PYTHONUTF8']='1'
    policy.preflight()
    initial=evidence()
    identity=dict(scope='JE_PAPER_INTERRUPTION_RESUME_20261020',evidence=initial,
        scripts={str(Path(p).resolve()):r.sha(Path(p)) for p in (__file__,policy.__file__,r.__file__)},
        unchanged_resource_policy=dict(ooc_core_mb=8000,threads=r.THREADS,timeout=r.TIMEOUT,
            start_disk=r.START_DISK,stop_disk=r.STOP_DISK,stop_available_ram=r.MIN_AVAILABLE_RAM),
        original_results_preserved=True,physics_changed=False,thresholds_changed=False)
    OUTPUT.mkdir(exist_ok=True)
    path=OUTPUT/'runner_manifest.json'
    if path.exists():
        if r.read(path)!=identity:raise RuntimeError('Resume identity changed')
    else:r.save(path,identity)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    r.save(OUTPUT/f'batch_{stamp}_plan.json',plan)
    events=[]
    r.OOC_MB=8000
    try:
        for job in jobs:
            policy.preflight()
            if evidence()!=initial:raise RuntimeError('Gate or source evidence changed')
            try:prior=resolve_complete(job['name'])
            except RuntimeError as exc:
                if not str(exc).startswith('No completed same-point source:'):raise
                prior=None
            if prior is not None:
                manifest=r.read(prior/'manifest.json')
                if manifest.get('counts')!=job['counts'] or manifest.get('supg_scale')!=job['supg_scale'] or manifest.get('voltages')!=[job['voltage']]:
                    raise RuntimeError('Reuse does not match requested point')
                events.append(dict(name=job['name'],state='EXPLICIT_REFERENCE_REUSE',directory=str(prior)))
                print('REUSED '+job['name']+' from '+str(prior),flush=True);continue
            dest=Path(job['directory'])
            if dest.exists():raise RuntimeError('Incomplete attempt preserved; use a new resume root')
            if shutil.disk_usage(OUTPUT).free<r.START_DISK:raise RuntimeError('RESOURCE_STOP_BEFORE_PREPARATION')
            candidates=[root for root in (policy.ORIGINAL,PRIOR) if
                (root/'meshes'/r.name(job['counts'],0).split('_S')[0]).exists()]
            mesh=r.ensure_mesh(candidates[0] if candidates else OUTPUT,job['counts'])
            mapped=transfer(job,mesh)
            r.run_native(job,mesh,mapped/'nodes.csv',dest,initial['gate'])
            subprocess.run([sys.executable,'-X','utf8',str(r.ROOT/'scripts/audit_je_sensitivity_power_20261015.py'),
                '--runs',str(dest)],check=True)
            fields=r.field_comparison(mapped,Path(job['source']),dest)
            r.save(dest/'mapped_field_comparison.json',fields)
            r.save(dest/'completion.json',dict(power_audit_completed=True,field_comparison_passed=fields['passed'],
                fully_validated=False,source_run=job['source'],mapping_directory=str(mapped),
                geometry_and_transfer_only_initial_guess=True))
            if job['name']==INTERRUPTED.name and r.sha(dest/'input.nml')!=r.sha(INTERRUPTED/'input.nml'):
                raise RuntimeError('Resumed 8V input differs from interrupted input')
            events.append(dict(name=job['name'],state=r.PASS,directory=str(dest),local_comparison_passed=fields['passed']))
            print('COMPLETED '+job['name']+'; local comparison passed='+str(fields['passed']),flush=True)
        r.save(OUTPUT/f'batch_{stamp}_audits.json',batch_audit(args.batch,args.joint_level))
        r.save(OUTPUT/f'batch_{stamp}_finished.json',dict(all_requested_attempts_completed=True,events=events,fully_validated=False))
    except BaseException as exc:
        r.save(OUTPUT/f'batch_{stamp}_stopped.json',dict(reason=repr(exc),events=events,fully_validated=False))
        raise


if __name__=='__main__':main()
