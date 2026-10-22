"""Finite follow-through to 8 V and representative studies; never certifies missing fields."""
from pathlib import Path
import json, subprocess, sys, time

BASE=Path(__file__).absolute().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_hv_sweep import checkpoint_seed, sha

OUT=BASE/'r25_high_voltage_followthrough'
EXE=ROOT/'build/HV_LOG_INVERSE_20261008_R25_CLEAN/TPSolverV2.exe'
A64=ROOT/'cases/HV_CONTROLLED_ANGLE_20261008/A64_half/H1_Y_HALF.mesh'
A48=ROOT/'cases/HV_CONTROLLED_ANGLE_20261008/A48_half/H1_Y_HALF.mesh'
A32=Path(json.loads((BASE/'r25_log_L2_V060_replay/manifest.json').read_text())['mesh'])
OUT.mkdir(exist_ok=False)
state=dict(status='WAITING_FOR_LOWER_VOLTAGE_QUEUE',complete_quantitative_0_8V=False,
           script_sha256=sha(Path(__file__)),executable_sha256=sha(EXE),operations=[])
def save():
    p=OUT/'progress.tmp';p.write_text(json.dumps(state,indent=2));p.replace(OUT/'progress.json')
def command(script,args):return [sys.executable,str(ROOT/'scripts'/script)]+list(map(str,args))
def run(script,args):
    cmd=command(script,args);i=len(state['operations'])
    item=dict(command=cmd,status='RUNNING');state['operations'].append(item);save()
    print('START',script,flush=True)
    with (OUT/f'operation_{i:03d}.txt').open('w') as stream:
        rc=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT).returncode
    item.update(status='FINISHED',exit_code=rc);save();print('END',script,rc,flush=True)
    return rc
def point(mesh,dest,v,seed,seed_v=None):
    return ['--exe',EXE,'--mesh',mesh,'--output',dest,'--voltage',repr(v),
            '--seed-nodes',seed,'--seed-voltage',repr(v if seed_v is None else seed_v),
            '--half-y','--log-weighted','--storage','OOC','--ooc-core-mb','5000',
            '--split-holes','--omp-threads','2','--mkl-threads','2']
def rel(path):return str(path.relative_to(BASE))

high=None
try:
    save();deadline=time.monotonic()+20000
    while True:
        try:previous=json.loads((BASE/'r25_followthrough/progress.json').read_text())
        except (FileNotFoundError,json.JSONDecodeError):previous={}
        if (previous.get('status') in ('FINITE_QUEUE_FINISHED_NOT_INTERVAL_CERTIFICATE','STOPPED_WITH_ERROR')
            and (BASE/'r25_log_A64_path_V005_V060/execution.json').exists()):break
        if time.monotonic()>deadline:raise TimeoutError('Lower-voltage work is not complete; no third large solve started')
        time.sleep(15)
    state['status']='RUNNING_HIGH_VOLTAGE_AND_REPRESENTATIVES';save()
    plan=json.loads((BASE/'r25_evidence_plan_stage3.json').read_text())
    for row in previous.get('points',[]):
        v=row['voltage_V']
        if v not in (0.5,6.0):
            plan['points']=[p for p in plan['points'] if p['voltage_V']!=v]
            plan['points'].append(dict(voltage_V=v,directory=rel(Path(row['directory']))))
    for mapping in [BASE/'r25_followthrough/radial_R1_R2_mapping']:
        if (mapping/'transfer_report.json').exists():
            plan['studies'].append(dict(name='6V fixed A64/z2 radial three-level study',kind='mesh',voltage_V=6,
                scope='Radial-only 14/16/18 comparison',runs=['r25_log_R0_A64_Z2_V060','r25_log_R1_A64_Z2_V060','r25_log_A64_Z2_V060']))
            plan['studies'].append(dict(name='6V radial R1-R2 local field',kind='field_shape',voltage_V=6,
                scope='2% maximum local field criterion',runs=['r25_log_R1_A64_Z2_V060','r25_log_A64_Z2_V060'],mapping=rel(mapping)))
    highdir=OUT/'target_V080'
    seed=BASE/'r25_followthrough/V060/nodes.csv'
    if not seed.exists():seed=BASE/'r25_log_A64_Z2_V060/nodes.csv'
    highcmd=command('run_hv_point.py',point(A64,highdir,8.,seed,6.)+['--checkpoint-volts','0.5','--max-seconds','14400'])
    state['high_voltage_command']=highcmd;save()
    highlog=(OUT/'high_voltage_driver.txt').open('w')
    high=subprocess.Popen(highcmd,stdout=highlog,stderr=subprocess.STDOUT)
    state['high_voltage_driver_pid']=high.pid;save()
    for v in (3.5,5.0):
        label=f'V{int(v*10):03d}';fine=BASE/'r25_followthrough'/label
        if not (fine/'nodes.csv').exists():continue
        coarse=[]
        for name,mesh in [('A32',A32),('A48',A48)]:
            mapping=OUT/(label+'_'+name+'_seed');dest=OUT/(label+'_'+name)
            rc=run('transfer_hv_initial_guess.py',['--source-mesh',A64,'--source-nodes',fine/'nodes.csv',
                '--target-mesh',mesh,'--output',mapping])
            if rc:break
            run('run_hv_point.py',point(mesh,dest,v,mapping/'nodes.csv'));coarse.append(dest)
        if len(coarse)==2:
            plan['studies'].append(dict(name=f'{v:g}V controlled angular32/48/64',kind='mesh',voltage_V=v,
                scope='Fixed radial and thickness layout',runs=[rel(p) for p in coarse+[fine]]))
            mapping=OUT/(label+'_A48_to_A64')
            if run('transfer_hv_initial_guess.py',['--source-mesh',A48,'--source-nodes',coarse[1]/'nodes.csv',
                    '--target-mesh',A64,'--output',mapping])==0:
                plan['studies'].append(dict(name=f'{v:g}V A48-A64 local field',kind='field_shape',voltage_V=v,
                    scope='2% maximum local field criterion',runs=[rel(coarse[1]),rel(fine)],mapping=rel(mapping)))
        perturbed=OUT/(label+'_Q5_seed');q5=OUT/(label+'_Q5')
        if run('perturb_hv_initial_guess.py',['--source-nodes',fine/'nodes.csv','--output',perturbed,
                '--temperature-amplitude-k','0.01'])==0:
            run('run_hv_point.py',point(A64,q5,v,perturbed/'nodes.csv')+['--quadrature-order','5'])
            plan['studies'].append(dict(name=f'{v:g}V A64 quadrature3/5',kind='quadrature',voltage_V=v,
                scope='Ten scalar/ring and seven local fields at 1%',include_local_fields=True,runs=[rel(fine),rel(q5)]))
    high.wait();highlog.close();state['high_voltage_driver_exit_code']=high.returncode;save()
    for v in (6.5,7.,7.5,8.):
        try:nodes,load=checkpoint_seed(highdir,v)
        except ValueError:continue
        dest=OUT/f'V{int(v*10):03d}'
        run('run_hv_point.py',point(A64,dest,v,nodes))
        plan['points']=[p for p in plan['points'] if p['voltage_V']!=v]
        plan['points'].append(dict(voltage_V=v,directory=rel(dest)))
    planpath=BASE/'r25_evidence_plan_after_8v_queue.json'
    with planpath.open('x') as stream:json.dump(plan,stream,indent=2)
    run('summarize_hv_credibility.py',['--plan',planpath,'--output',BASE/'delivery_r25_after_8v_queue'])
    state['status']='FINITE_QUEUE_FINISHED_REVIEW_REQUIRED'
    state['high_voltage_attempt']=json.loads((highdir/'summary.json').read_text())
except Exception as exc:
    state.update(status='STOPPED_WITH_ERROR',error=f'{type(exc).__name__}: {exc}')
    if high is not None:state['high_voltage_still_running']=high.poll() is None
    raise
finally:save()
