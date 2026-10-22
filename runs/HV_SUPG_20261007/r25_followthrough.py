"""Run-owned finite validation queue. No interpolation or interval certification."""
from pathlib import Path
import json, subprocess, sys, time

BASE=Path(__file__).absolute().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_hv_sweep import checkpoint_seed, sha

EXE=ROOT/'build/HV_LOG_INVERSE_20261008_R25_CLEAN/TPSolverV2.exe'
MESH=ROOT/'cases/HV_CONTROLLED_ANGLE_20261008/A64_half/H1_Y_HALF.mesh'
PATH_RUN=BASE/'r25_log_A64_path_V005_V060'
OUT=BASE/'r25_followthrough'
OUT.mkdir(exist_ok=False)
state=dict(status='RUNNING',complete_quantitative_0_8V=False,operations=[],points=[],
           script_sha256=sha(Path(__file__)),executable_sha256=sha(EXE))
def save():
    temporary=OUT/'progress.tmp'
    temporary.write_text(json.dumps(state,indent=2))
    temporary.replace(OUT/'progress.json')

def run(script,args):
    command=[sys.executable,str(ROOT/'scripts'/script)]+list(map(str,args))
    index=len(state['operations'])
    item=dict(command=command,status='RUNNING')
    state['operations'].append(item);save()
    print('START',script,flush=True)
    with (OUT/f'operation_{index:03d}.txt').open('w') as stream:
        result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
    item.update(status='FINISHED',exit_code=result.returncode);save()
    print('END',script,result.returncode,flush=True)
    return result.returncode

def point_args(mesh,dest,v,seed):
    return ['--exe',EXE,'--mesh',mesh,'--output',dest,'--voltage',repr(v),
            '--seed-nodes',seed,'--seed-voltage',repr(v),'--half-y','--log-weighted',
            '--storage','OOC','--ooc-core-mb','5000','--split-holes',
            '--omp-threads','2','--mkl-threads','2']

try:
    save()
    deadline=time.monotonic()+18000
    prerequisite=BASE/'r25_log_R1_A64_Z2_V060'
    while not (prerequisite/'summary.json').exists():
        if time.monotonic()>deadline:raise TimeoutError('R1 completion not recorded')
        time.sleep(15)
    r0=BASE/'r25_log_R0_A64_Z2_V060'
    run('run_hv_point.py',point_args(ROOT/'cases/HV_CONTROLLED_RADIAL_20261008/R0_half/H1_Y_HALF.mesh',
        r0,6.0,BASE/'seed_log_R0_A64_Z2_V060/nodes.csv'))
    run('audit_hv_mesh.py',['--levels',r0,prerequisite,BASE/'r25_log_A64_Z2_V060',
        '--voltage','6','--output',OUT/'radial_three_levels_V060.json'])
    if json.loads((prerequisite/'summary.json').read_text())['state']=='SOLVE_FLUX_PASS_MESH_PENDING':
        mapping=OUT/'radial_R1_R2_mapping'
        rc=run('transfer_hv_initial_guess.py',['--source-mesh',
            ROOT/'cases/HV_CONTROLLED_RADIAL_20261008/R1_half/H1_Y_HALF.mesh',
            '--source-nodes',prerequisite/'nodes.csv','--target-mesh',MESH,'--output',mapping])
        if rc==0:
            run('audit_hv_field_shape.py',['--mapping',mapping,'--source',prerequisite,
                '--target',BASE/'r25_log_A64_Z2_V060','--output',OUT/'radial_local_R1_R2_V060.json'])
    for v in [i/2 for i in range(1,13)]:
        while True:
            try:
                nodes,load=checkpoint_seed(PATH_RUN,v)
                break
            except ValueError:
                if (PATH_RUN/'execution.json').exists() or time.monotonic()>deadline:
                    raise RuntimeError(f'No usable accepted checkpoint at {v:g} V')
                time.sleep(15)
        dest=OUT/f'V{int(v*10):03d}'
        record=dict(voltage_V=v,directory=str(dest),checkpoint_nodes=str(nodes),
            checkpoint_nodes_sha256=sha(nodes),checkpoint_load=str(load),
            checkpoint_load_sha256=sha(load),status='RUNNING')
        state['points'].append(record);save()
        rc=run('run_hv_point.py',point_args(MESH,dest,v,nodes))
        summary=dest/'summary.json'
        record.update(status=json.loads(summary.read_text())['state'] if summary.exists() else 'FAILED_NO_SUMMARY',exit_code=rc)
        save()
        if rc not in (0,1):raise RuntimeError('Replay process failed unexpectedly')
    state['status']='FINITE_QUEUE_FINISHED_NOT_INTERVAL_CERTIFICATE'
except Exception as exc:
    state.update(status='STOPPED_WITH_ERROR',error=f'{type(exc).__name__}: {exc}')
    raise
finally:
    save()
