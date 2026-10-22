"""Stream archived mechanical observables without changing frozen metrics."""
import csv, hashlib, json, math, time
from pathlib import Path
ROOT=Path(r'C:\Users\Administrator\Desktop\codex-workspaces\Console3\solver_v2')
OUT=Path(__file__).resolve().parent
RUNS=['r17_scale1_V060','r17_XY1_Z3_V060','r17_XY2_Z3_V060_OOC']
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
results=[]; start=time.monotonic()
for run in RUNS:
    src=ROOT/'runs/HV_SUPG_20261007'/run/'field_samples.csv'
    total=ring=0.; vu=vs=ru=rs=0.; vmax=0.; count=0
    with src.open(newline='') as f:
        reader=csv.reader(f); header=next(reader)
        for row in reader:
            x,y,z,w=map(float,row[:4]); u=math.sqrt(sum(float(v)**2 for v in row[11:14]))
            xx,yy,zz,yz,zx,xy=map(float,row[15:21])
            vm=math.sqrt(.5*((xx-yy)**2+(yy-zz)**2+(zz-xx)**2)+3*(yz*yz+zx*zx+xy*xy))
            total+=w; vu+=w*u; vs+=w*vm; vmax=max(vmax,vm); count+=1
            if 3e-6<=math.hypot(x,y)<=3.5e-6:
                ring+=w;ru+=w*u;rs+=w*vm
    results.append(dict(run=run,path=str(src),sha256=sha(src),sample_count=count,
        volume_mean_displacement_m=vu/total,volume_mean_von_mises_Pa=vs/total,
        ring_mean_displacement_m=ru/ring,ring_mean_von_mises_Pa=rs/ring,
        raw_sample_von_mises_max_Pa=vmax))
    print(run,results[-1]['ring_mean_von_mises_Pa'],flush=True)
report=dict(scope='ARCHIVED_MECHANICAL_DIAGNOSTICS_NOT_MESH_CERTIFICATE',script_sha256=sha(Path(__file__)),
    elapsed_seconds=time.monotonic()-start,runs=results,
    limitations=['R17 runs share full-domain rigid-body anchors; do not compare raw displacement directly to R21 half-domain gauge.',
    'Stress raw maxima are diagnostic, not a converged stress-concentration certificate.',
    'Mechanical integral stability never overrides the failed carrier mesh/SUPG gate.'])
with (OUT/'mechanical_diagnostic.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
