"""Finite fixed-mesh high-voltage replay queue, with no mesh-refinement tasks."""
from pathlib import Path
import csv,json,math,subprocess,sys,time

HERE=Path(__file__).absolute().parent
BASE=HERE.parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_hv_sweep import checkpoint_seed,read_metrics,sha
from summarize_hv_credibility import point_record

OUT=HERE/'results';OUT.mkdir(exist_ok=False)
EXE=ROOT/'build/HV_LOG_INVERSE_20261008_R25_CLEAN/TPSolverV2.exe'
MESH=ROOT/'cases/HV_CONTROLLED_ANGLE_20261008/A64_half/H1_Y_HALF.mesh'
SOURCES={'full_heat_voltage':BASE/'r25_high_voltage_followthrough/target_V080',
         'voltage_then_heat':BASE/'r25_log_A64_voltage_then_heat_V080'}
GRID=[6.,6.25,6.5,7.,7.5,8.]
records=[];processed=set()
state=dict(status='RUNNING',scope='FIXED_EXISTING_MESH_NO_SPATIAL_CERTIFICATE',
           script_sha256=sha(Path(__file__)),mesh_sha256=sha(MESH),
           executable_sha256=sha(EXE),records=records,attempts={})

def add_record(v,directory,source):
    record=point_record(v,directory,directory/'manifest.json')
    record['route']=source;records.append(record)

def save():
    state['all_requested_solve_flux_pass']=all(any(r['voltage_V']==v and
        r['state']=='SOLVE_FLUX_PASS_MESH_PENDING' for r in records) for v in [6.,6.5,7.,7.5,8.])
    state['complete_quantitative_0_8V']=False
    temporary=OUT/'progress.tmp';temporary.write_text(json.dumps(state,indent=2,allow_nan=False))
    temporary.replace(OUT/'progress.json')
    with (OUT/'status_6_8V.csv').open('w',newline='',encoding='utf-8-sig') as stream:
        keys=['voltage_V','state','theta_max_K','joule_half_W','electron_right_half_A',
              'physical_exterior_heat_balance_relative','directory']
        writer=csv.DictWriter(stream,fieldnames=keys);writer.writeheader()
        for v in GRID:
            available=[r for r in records if r['voltage_V']==v]
            preferred=[r for r in available if r['state']=='SOLVE_FLUX_PASS_MESH_PENDING']
            if not available:
                missing='NO_FULL_HEAT_RESULT_OBTAINED' if state['status']=='FINISHED_FIXED_MESH_REVIEW_REQUIRED' else 'NO_FULL_HEAT_RESULT_YET'
                writer.writerow(dict(voltage_V=v,state=missing));continue
            r=(preferred or available)[-1];m=r.get('metrics',{})
            writer.writerow(dict(voltage_V=v,state=r['state'],directory=r['directory'],
                theta_max_K=m.get('theta_max_K'),joule_half_W=m.get('joule_W'),
                electron_right_half_A=m.get('electron_right_A'),
                physical_exterior_heat_balance_relative=m.get('physical_exterior_heat_balance_relative')))

def replay(route,source,v):
    nodes,load=checkpoint_seed(source,v)
    dest=OUT/f'{route}_V{round(v*100):04d}'
    cmd=[sys.executable,str(ROOT/'scripts/run_hv_point.py'),'--exe',str(EXE),'--mesh',str(MESH),
        '--output',str(dest),'--voltage',repr(v),'--seed-nodes',str(nodes),'--seed-voltage',repr(v),
        '--half-y','--log-weighted','--storage','OOC','--ooc-core-mb','5000','--split-holes',
        '--omp-threads','2','--mkl-threads','2']
    operation=dict(route=route,voltage_V=v,command=cmd,checkpoint_load_sha256=sha(load),
                   checkpoint_nodes_sha256=sha(nodes))
    state['current_replay']=operation;save();print('REPLAY',route,v,flush=True)
    with (OUT/f'{route}_V{round(v*100):04d}_driver.txt').open('w') as stream:
        rc=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT).returncode
    operation['exit_code']=rc
    add_record(v,dest,route)
    records[-1]['checkpoint_provenance']=operation
    state.pop('current_replay');save()

try:
    add_record(6.,BASE/'r25_followthrough/V060','completed_full_heat_path')
    existing=BASE/'r25_log_A64_fixed_V0625'
    if (existing/'execution.json').exists():add_record(6.25,existing,'full_heat_voltage')
    deadline=time.monotonic()+32400
    while True:
        for route,source in SOURCES.items():
            info=dict(directory=str(source),finished=(source/'execution.json').exists(),
                      accepted_full_heat_voltages=[])
            for f in sorted(source.glob('accepted_*_load.txt')):
                try:
                    load=read_metrics(f);v=abs(load['applied_voltage_V'])
                    if not math.isfinite(v):continue
                    info['last_checkpoint_voltage_V']=v
                    info['last_checkpoint_heat_factor']=load.get('joule_factor')
                    checkpoint_seed(source,v)
                    if v not in info['accepted_full_heat_voltages']:info['accepted_full_heat_voltages'].append(v)
                except (OSError,ValueError,KeyError):continue
            if (source/'summary.json').exists():info['native_summary']=json.loads((source/'summary.json').read_text())
            state['attempts'][route]=info
        save()
        # A replay may allocate another large system. Wait until its source solve
        # has exited, so only the other continuation and this replay can be active.
        for route,source in SOURCES.items():
            if not state['attempts'][route]['finished']:continue
            for v in GRID[1:]:
                if (route,v) in processed:continue
                if v==6.25 and any(r['voltage_V']==v and r['route']==route for r in records):
                    processed.add((route,v));continue
                try:checkpoint_seed(source,v)
                except ValueError:processed.add((route,v));continue
                replay(route,source,v);processed.add((route,v))
        if all(i['finished'] for i in state['attempts'].values()):break
        if time.monotonic()>deadline:raise TimeoutError('High-voltage attempts did not finish; incomplete states remain unaccepted')
        time.sleep(20)
    state['status']='FINISHED_FIXED_MESH_REVIEW_REQUIRED'
except Exception as exc:
    state.update(status='STOPPED_WITH_ERROR',error=f'{type(exc).__name__}: {exc}')
    raise
finally:save()
