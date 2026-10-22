"""Small one-element scheduling check, not an A64/high-voltage validation."""
import csv,json,subprocess,sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).absolute().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from run_hv_point import apply_load_path,apply_trust_region,apply_heat_restart
from run_hv_sweep import make_input,read_metrics,sha
from run_hv_regression import runtime_env

OUT=HERE/'heat_restart_native_check_v2';OUT.mkdir(exist_ok=False)
EXE=ROOT/'build/HV_TRUST_REGION_20261008_R26C/TPSolverV2.exe'
source=ROOT/'cases/validation/G08_reference_bar/mesh.dat'
lines=source.read_text().splitlines();generated=[lines[0]]
for line in lines[1:21]:
    node,x,y,z=line.split();xx,yy,zz=float(x)*20-1e-5,float(y)*20-1e-5,float(z)-.5e-6
    generated.append(f'{node} {xx:.17e} {yy:.17e} {zz:.17e}')
generated+=lines[21:];mesh=OUT/'small_plain_plate.mesh';mesh.write_text('\n'.join(generated)+'\n')

def run(name,nodes=None,heat=None):
    d=OUT/name;d.mkdir()
    text=make_input(mesh,.5,nodes.parent if nodes else None,.5 if nodes else 0.,1.,False)
    text=apply_trust_region(apply_load_path(text,True),True)
    text=text.replace("problem_type='CENTER_HOLE_PLATE'","problem_type='PLAIN_PLATE_CONTROL'").replace('hole_radius=3e-6','hole_radius=0')
    text=text.replace('transport_supg=.true.,','transport_supg=.false.,transport_log_weighted=.true.,')
    if nodes:
        text=text.replace((nodes.parent/'nodes.csv').as_posix(),nodes.as_posix())
        text=apply_heat_restart(text,heat,.125,.25)
    inp=d/'input.nml';inp.write_text(text)
    env=runtime_env();env.update(OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    with (d/'stdout.txt').open('wb') as stream:
        rc=subprocess.run([str(EXE),'run-f1',str(inp),str(d)],env=env,stdout=stream,stderr=subprocess.STDOUT).returncode
    (d/'execution.json').write_text(json.dumps(dict(exit_code=rc,executable_sha256=sha(EXE),input_sha256=sha(inp))))
    return d,rc

baseline,rc0=run('continuous')
seeds=[]
for path in sorted(baseline.glob('accepted_*_load.txt')):
    load=read_metrics(path)
    if abs(abs(load['applied_voltage_V'])-.5)<1e-12 and 0<load['joule_factor']<1:
        seeds.append((path,load))
assert seeds,'No converged intermediate heat checkpoint'
loadpath,load=seeds[0]
resumed,rc1=run('resumed',Path(str(loadpath).replace('_load.txt','_nodes.csv')),load['joule_factor'])
arrays=[np.genfromtxt(d/'nodes.csv',delimiter=',',names=True) for d in [baseline,resumed]]
diff={}
for key,scale in zip(('ux_m','uy_m','uz_m','phi_V','eta_n','eta_p','theta_K'),(1e-10,1e-10,1e-10,1.,.1,.1,35.)):
    a,b=(arr[key] for arr in arrays)
    diff[key]=float(np.max(np.abs(a-b))/scale)
metrics=[read_metrics(d/'metrics.txt') for d in [baseline,resumed]]
complete=all(m.get('final_voltage_factor')==1. and m.get('final_joule_factor')==1. and
             m.get('maximum_accepted_residual',float('inf'))<=1e-9 and
             m.get('maximum_accepted_update',float('inf'))<=1e-9 for m in metrics)
report=dict(scope='ONE_ELEMENT_SCHEDULING_REGRESSION_NOT_HIGH_VOLTAGE_EVIDENCE',nodes=20,elements=1,
            exit_codes=[rc0,rc1],source_heat_factor=load['joule_factor'],scaled_field_difference=diff,
            native_complete=complete,native_metrics=metrics,passed=complete and max(diff.values())<1e-8)
(OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False))
print(json.dumps({k:v for k,v in report.items() if k!='native_metrics'}))
raise SystemExit(0 if report['passed'] else 1)
