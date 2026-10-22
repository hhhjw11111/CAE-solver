"""One adaptive 0--8 V path, exact half-volt checkpoints, independent point audits.

Preserves continuation predictor history between output voltages. Each saved grid
point is replayed with the same executable at its true voltage and full Joule heat
for independent flux/field checks. A checkpoint alone never becomes an accepted point.
"""
from pathlib import Path
import argparse,csv,json,os,subprocess,time,math
from run_hv_sweep import ROOT,FIELDS,make_input,sha,read_metrics,classify,audit_nodes
from run_hv_regression import runtime_env

def checkpoint_voltage(meta):
    m=read_metrics(meta)
    if not all(isinstance(m.get(k),(float,int)) and math.isfinite(m[k]) for k in
               ['applied_voltage_V','joule_factor','residual','update']):return None
    v=-m['applied_voltage_V'];grid=round(v*2)/2
    if (m['joule_factor']!=1 or m['residual']>=1e-9 or m['update']>=1e-9 or
        not 0<=v<=8 or abs(v-grid)>1e-10):return None
    return grid

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe',type=Path,required=True);p.add_argument('--mesh',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--supg-scale',type=float,default=1.)
    p.add_argument('--seed-nodes',type=Path);p.add_argument('--seed-voltage',type=float,default=0.)
    p.add_argument('--mkl-threads',type=int,default=4);p.add_argument('--max-seconds',type=float,default=28800)
    p.add_argument('--audit-during-trajectory',action='store_true',
        help='Allow simultaneous solver processes. Default defers flux audits to keep peak memory bounded.')
    a=p.parse_args()
    a.exe=Path(os.path.abspath(a.exe));a.mesh=Path(os.path.abspath(a.mesh));a.output=Path(os.path.abspath(a.output))
    if not 0<=a.seed_voltage<=.5 or (a.seed_voltage>0 and a.seed_nodes is None):
        raise ValueError('A full grid trajectory may only seed at 0--0.5 V.')
    if a.seed_nodes:a.seed_nodes=Path(os.path.abspath(a.seed_nodes))
    a.output.mkdir(parents=True,exist_ok=False);trajectory=a.output/'trajectory';trajectory.mkdir()
    env=runtime_env();env['MKL_NUM_THREADS']=str(a.mkl_threads)
    manifest=dict(executable=str(a.exe),executable_sha256=sha(a.exe),mesh=str(a.mesh),mesh_sha256=sha(a.mesh),
        voltages=[i/2 for i in range(17)],supg_scale=a.supg_scale,model='FROZEN_DRIFT_HEAT_ZNO_NT',
        acceptance='SOLVE_AND_FLUX_ONLY; MESH_AND_STABILIZATION_STUDIES_REQUIRED',
        voltage_checkpoint_spacing=.0625,MKL_NUM_THREADS=a.mkl_threads,OMP_NUM_THREADS=env['OMP_NUM_THREADS'],
        audit_during_trajectory=a.audit_during_trajectory)
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
    rows=[dict(voltage=v,state='WAITING_FOR_CHECKPOINT',exit_code='') for v in manifest['voltages']]
    def save():
        with (a.output/'status.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=FIELDS);w.writeheader();w.writerows(rows)
    def audit_point(v,seed,metadata=None):
        dest=a.output/f'V{round(v*10):03d}';dest.mkdir();inp=dest/'input.nml'
        text=make_input(a.mesh,v,seed.parent if seed and v>0 else None,v,a.supg_scale,False)
        if seed and v>0:text=text.replace((seed.parent/'nodes.csv').as_posix(),seed.as_posix())
        inp.write_text(text,encoding='ascii');i=round(v*2);rows[i]['state']='INDEPENDENT_AUDIT_RUNNING';save()
        if seed:
            (dest/'seed.json').write_text(json.dumps(dict(path=str(seed),sha256=sha(seed),
                status='REPLAY_OF_ALGEBRAIC_CHECKPOINT_REQUIRES_INDEPENDENT_AUDIT',
                original_load=read_metrics(metadata) if metadata else None),indent=2))
        start=time.monotonic();timed_out=False
        with (dest/'stdout.txt').open('wb') as f:
            proc=subprocess.Popen([str(a.exe),'run-f1',str(inp),str(dest)],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
            (dest/'process.json').write_text(json.dumps(dict(pid=proc.pid)))
            try:rc=proc.wait(timeout=1800)
            except subprocess.TimeoutExpired:proc.kill();rc=proc.wait();timed_out=True
        (dest/'execution.json').write_text(json.dumps(dict(exit_code=rc,seconds=time.monotonic()-start,
            timed_out=timed_out,input_sha256=sha(inp),executable_sha256=sha(a.exe)),indent=2))
        state,m=classify(dest,rc)
        if state!='FAILED_NO_COMPLETE_STATE':
            fields=audit_nodes(dest/'nodes.csv',a.mesh,v)
            (dest/'field_audit.json').write_text(json.dumps(fields,indent=2,allow_nan=False))
            if not fields['passed']:state='FAILED_FIELD_AUDIT'
        rows[i].update({k:m[k] for k in FIELDS if k in m});rows[i].update(state=state,exit_code=rc)
        save();print(f'GRID {v:g} V: {state}',flush=True)
    save()
    # Zero is solved and independently audited, even for a nonzero seed.
    audit_point(0.,None)
    if rows[0]['state']!='SOLVE_FLUX_PASS_MESH_PENDING':
        (a.output/'summary.json').write_text(json.dumps(dict(all_requested_solve_flux_pass=False,
            reason='ZERO_STATE_AUDIT_FAILED',mesh_validated=False,complete_quantitative_0_8V=False),indent=2))
        return 1
    text=make_input(a.mesh,8,a.seed_nodes.parent if a.seed_nodes else None,a.seed_voltage,a.supg_scale,False)
    if a.seed_nodes:text=text.replace((a.seed_nodes.parent/'nodes.csv').as_posix(),a.seed_nodes.as_posix())
    text=text.replace('continuation_initial_step=0.03125,','continuation_initial_step=0.0625,')
    text=text.replace('continuation_secant_predictor=.true.,',
        'continuation_secant_predictor=.true.,continuation_voltage_checkpoint_spacing=0.0625,')
    inp=trajectory/'input.nml';inp.write_text(text,encoding='ascii');start=time.monotonic();timed_out=False
    if a.seed_nodes:(trajectory/'seed.json').write_text(json.dumps(dict(path=str(a.seed_nodes),
        sha256=sha(a.seed_nodes),voltage=a.seed_voltage,status='INITIAL_GUESS_REQUIRES_REEQUILIBRATION'),indent=2))
    with (trajectory/'stdout.txt').open('wb') as f:
        process=subprocess.Popen([str(a.exe),'run-f1',str(inp),str(trajectory)],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        (trajectory/'process.json').write_text(json.dumps(dict(pid=process.pid)))
        seen=set();pending={}
        while True:
            returncode=process.poll()
            for meta in sorted(trajectory.glob('accepted_*_load.txt')):
                if meta.name in seen:continue
                v=checkpoint_voltage(meta)
                # A file may be visible before its writer closes it; retry it.
                if v is None:continue
                seed=meta.with_name(meta.name.replace('_load.txt','_nodes.csv'))
                if not seed.exists():continue
                seen.add(meta.name)
                if v>0 and rows[round(v*2)]['state']=='WAITING_FOR_CHECKPOINT':
                    if a.audit_during_trajectory:
                        audit_point(v,seed,meta)
                    else:
                        pending[v]=(seed,meta)
                        rows[round(v*2)]['state']='CHECKPOINT_AWAITING_FLUX_AUDIT'
                        save();print(f'GRID {v:g} V: algebraic checkpoint; flux audit deferred',flush=True)
            if returncode is not None:break
            if time.monotonic()-start>a.max_seconds:
                process.kill();process.wait();timed_out=True;continue
            time.sleep(3)
        rc=process.wait()
    (trajectory/'execution.json').write_text(json.dumps(dict(exit_code=rc,seconds=time.monotonic()-start,
        timed_out=timed_out,input_sha256=sha(inp),executable_sha256=sha(a.exe)),indent=2))
    for v,(seed,meta) in sorted(pending.items()):audit_point(v,seed,meta)
    for row in rows:
        if row['state']=='WAITING_FOR_CHECKPOINT':row['state']='NO_ACCEPTED_CHECKPOINT'
    save();passed=all(r['state']=='SOLVE_FLUX_PASS_MESH_PENDING' for r in rows)
    (a.output/'summary.json').write_text(json.dumps(dict(all_requested_solve_flux_pass=passed,
        trajectory_exit_code=rc,mesh_validated=False,complete_quantitative_0_8V=False),indent=2))
    return 0 if passed else 1

if __name__=='__main__':raise SystemExit(main())
