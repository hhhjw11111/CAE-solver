"""Isolated fixed-physics JE sensitivity solves; frozen source and old runs unchanged."""
from pathlib import Path
import argparse,ctypes,json,re,subprocess,time,sys
from run_hv_sweep import sha,classify,audit_nodes
from run_hv_regression import runtime_env
from hv_runtime import configure_storage,recorded_runtime

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'runs/JE_6_8_HALF_V_20261013'
EXE=ROOT/'build/ARC_JE_FINAL_20261012/TPSolverV2.exe'
MESH=ROOT/'cases/HV_SURFACE_20261008/XY0_Z2/H1_XY0_Z2.mesh'
STUDY=ROOT/'runs/JE_SENSITIVITY_20261015'


def short(path):
    path=Path(path).resolve();buf=ctypes.create_unicode_buffer(32768)
    if not ctypes.windll.kernel32.GetShortPathNameW(str(path),buf,len(buf)):
        raise OSError('Cannot resolve existing path '+str(path))
    value=buf.value.replace('\\','/');value.encode('ascii');return value


def replace(text,key,value):
    text,count=re.subn(r'(?<![A-Za-z0-9_])'+key+r"\s*=\s*(?:'[^']*'|[^,\n/]+)",key+'='+value,text)
    if count!=1:raise ValueError((key,count))
    return text


def run(voltage,scale,mesh,seed,out,threads=4,timeout_seconds=1800):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    mesh=Path(mesh).resolve();seed=Path(seed).resolve()
    text=(BASE/f'input_V{voltage*10:03d}.nml').read_text(encoding='utf-8-sig')
    changes={'case_id':"'JE_SENS_V%d_S%g'"%(voltage,scale),'mesh_file':"'"+short(mesh)+"'",
             'restart_nodes_file':"'"+short(seed)+"'",'restart_voltage_factor':'1.0','supg_scale':str(scale)}
    for k,v in changes.items():text=replace(text,k,v)
    assert 'heat_source_current_work=.true.' in text
    assert 'residual_tolerance=1e-9,update_tolerance=1e-9' in text
    assert 'target_joule_factor=1.0' in text and 'restart_joule_factor=1.0' in text
    inp=out/'input.nml';inp.write_text(text,encoding='ascii')
    assert threads in [4,8]
    assert timeout_seconds>0
    env=runtime_env();env.update(OMP_NUM_THREADS=str(threads),MKL_NUM_THREADS=str(threads))
    env=configure_storage(env,Path(short(out)),'OOC',4000,False)
    manifest={'executable':str(EXE),'executable_sha256':sha(EXE),'mesh':str(mesh),'mesh_sha256':sha(mesh),
              'model':'JE_CURRENT_WORK_DIAGNOSTIC_ZNO_NT','voltages':[voltage],'supg_scale':scale,
              'domain_symmetry':'FULL','log_weighted_quadrature_order':3,'runtime_env':recorded_runtime(env),
              'seed_nodes':str(seed),'seed_sha256':sha(seed),'initial_guess_status':'REQUIRES_REEQUILIBRATION',
              'heat_source_unchanged':True,'input_template':str(BASE/f'input_V{voltage*10:03d}.nml'),
              'changed_input_keys':list(changes),'wall_time_limit_seconds':timeout_seconds}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    start=time.monotonic();timed_out=False
    print('START '+out.name,flush=True)
    with (out/'stdout.txt').open('xb') as log:
        process=subprocess.Popen([short(EXE),'run-f1',short(inp),short(out)],cwd=short(ROOT),env=env,stdout=log,stderr=subprocess.STDOUT)
        (out/'process.json').write_text(json.dumps({'pid':process.pid}),encoding='utf-8')
        try:rc=process.wait(timeout=timeout_seconds)
        except subprocess.TimeoutExpired:process.kill();rc=process.wait();timed_out=True
    execution={'exit_code':rc,'seconds':time.monotonic()-start,'timed_out':timed_out,
               'input_sha256':sha(inp),'executable_sha256':sha(EXE)}
    (out/'execution.json').write_text(json.dumps(execution,indent=2),encoding='utf-8')
    status,met=classify(out,rc)
    if (out/'nodes.csv').exists():
        fields=audit_nodes(out/'nodes.csv',mesh,voltage)
        (out/'field_audit.json').write_text(json.dumps(fields,indent=2),encoding='utf-8')
        if not fields['passed']:status='FAILED_FIELD_AUDIT'
    if met and met.get('heat_source')!='E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE':
        raise RuntimeError('Incorrect heat model')
    summary={'state':status,'mesh_validated':False,'power_gap_validated':False,'elapsed_seconds':execution['seconds']}
    (out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print('END '+out.name+' '+status+' '+str(round(execution['seconds'],1))+'s',flush=True)
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base-batch',action='store_true');p.add_argument('--mesh-batch',action='store_true');p.add_argument('--voltage',type=int,choices=[6,8])
    p.add_argument('--scale',type=float,default=1.);p.add_argument('--mesh',type=Path,default=MESH)
    p.add_argument('--seed',type=Path);p.add_argument('--output',type=Path)
    p.add_argument('--study-root',type=Path,help='Fresh root for a reproducible batch; never overwrites existing runs')
    p.add_argument('--threads',type=int,choices=[4,8],default=4)
    p.add_argument('--timeout-seconds',type=int,default=1800,help='Resource limit only; never a numerical acceptance tolerance')
    a=p.parse_args()
    if a.study_root is not None:STUDY=a.study_root.resolve()
    if a.base_batch:
        for scale in [1.,.5,2.]:
            for v in [6,8]:
                run(v,scale,MESH,BASE/f'V{v*10:03d}'/'nodes.csv',STUDY/f'A16_S{scale:g}_V{v*10:03d}',threads=a.threads,timeout_seconds=a.timeout_seconds)
    elif a.mesh_batch:
        for count,previous in [(24,16),(32,24)]:
            mesh=STUDY/f'mesh_A{count}'/'H1_XY0_Z2.mesh'
            for v in [6,8]:
                source=STUDY/f'A{previous}_S1_V{v*10:03d}'
                assert json.loads((source/'summary.json').read_text())['state']=='SOLVE_FLUX_PASS_MESH_PENDING'
                source_manifest=json.loads((source/'manifest.json').read_text())
                mapping=STUDY/f'transfer_A{previous}_A{count}_V{v*10:03d}'
                if not (mapping/'transfer_report.json').exists():
                    subprocess.run([sys.executable,str(ROOT/'scripts/transfer_hv_initial_guess.py'),
                                    '--source-mesh',source_manifest['mesh'],'--source-nodes',str(source/'nodes.csv'),
                                    '--target-mesh',str(mesh),'--output',str(mapping)],check=True)
                report=json.loads((mapping/'transfer_report.json').read_text())
                assert report['source_nodes_sha256']==sha(source/'nodes.csv') and report['target_mesh_sha256']==sha(mesh)
                run(v,1.,mesh,mapping/'nodes.csv',STUDY/f'A{count}_S1_V{v*10:03d}',threads=a.threads,timeout_seconds=a.timeout_seconds)
    else:
        if a.voltage is None or a.seed is None or a.output is None or a.scale<=0:raise ValueError('Supply voltage, seed, output and positive scale')
        run(a.voltage,a.scale,a.mesh,a.seed,a.output,threads=a.threads,timeout_seconds=a.timeout_seconds)
