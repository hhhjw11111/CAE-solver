"""Read-only streamed diagnostics of archived fields; never changes acceptance."""
import csv, hashlib, json, math, time
from pathlib import Path

ROOT = Path(r'C:\Users\Administrator\Desktop\codex-workspaces\Console3\solver_v2')
OUT = Path(__file__).resolve().parent
RUNS = ['r17_scale1_V060', 'r17_XY1_Z3_V060', 'r17_XY2_Z3_V060_OOC']
DEPTH = [0, 2, 5, 10, 20, 50, 100, 250, 500.000001]

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()

def metrics(path):
    out = {}
    for line in path.read_text().splitlines():
        if '=' in line:
            k, v = line.split('=', 1)
            try: out[k] = float(v)
            except ValueError: out[k] = v
    return out

def add(d, key, w, p, n, q, theta):
    a = d.setdefault(key, [0., 0., 0., 0., 0., 0.])
    for i, v in enumerate([w, w*abs(p-1), w*(p-1), w*abs(n-1), w*q, w*theta]): a[i] += v

started = time.perf_counter()
reports = []
for run in RUNS:
    d = ROOT/'runs/HV_SUPG_20261007'/run
    src = d/'field_samples.csv'
    original = metrics(d/'metrics.txt')
    accum = {}; count = ring_count = 0
    with src.open(newline='') as f:
        reader = csv.reader(f)
        header = next(reader)
        if header[:4] != ['x_m','y_m','z_m','weight_m3']: raise ValueError('Bad header')
        for row in reader:
            x,y,z,w = map(float,row[:4]); q,theta,n,p = map(float,row[7:11])
            if not all(math.isfinite(v) for v in [x,y,z,w,q,theta,n,p]) or w <= 0 or min(n,p) <= 0:
                raise ValueError('Invalid sample')
            count += 1
            add(accum, 'all', w,p,n,q,theta)
            radius = math.hypot(x,y)
            if 3e-6 <= radius <= 3.5e-6:
                ring_count += 1; add(accum,'annulus',w,p,n,q,theta)
                depth = (.5e-6-abs(z))*1e9
                bi = next((i for i in range(len(DEPTH)-1) if DEPTH[i] <= depth < DEPTH[i+1]),None)
                if bi is None: raise ValueError('Depth outside domain')
                add(accum,f'depth_{DEPTH[bi]}_{DEPTH[bi+1]}nm',w,p,n,q,theta)
                add(accum,'side_top' if z >= 0 else 'side_bottom',w,p,n,q,theta)
                angle = (math.atan2(y,x)+2*math.pi)%(2*math.pi)
                add(accum,f'angular_sector_{int(angle/(math.pi/8)):02d}',w,p,n,q,theta)
                add(accum,'radial_inner' if radius < 3.25e-6 else 'radial_outer',w,p,n,q,theta)
    ring = accum['annulus']; analytic = math.pi*((3.5e-6)**2-(3e-6)**2)*1e-6
    rows = []
    for region,a in accum.items():
        rows.append(dict(region=region, volume_m3=a[0], mean_abs_p_relative=a[1]/a[0],
            mean_signed_p_relative=a[2]/a[0], mean_abs_n_relative=a[3]/a[0],
            Q_mean_W_m3=a[4]/a[0],theta_mean_K=a[5]/a[0],
            fraction_of_ring_p_integral=a[1]/ring[1] if region != 'all' else None))
    report = dict(run=run, source=str(src), source_sha256=sha(src),
        metrics_sha256=sha(d/'metrics.txt'), manifest=json.loads((d/'manifest.json').read_text()),
        sample_count=count,annulus_sample_count=ring_count, original_frozen_metrics=original,
        diagnostic_sampling='existing field_samples 5x5x5; frozen ring metric uses 5x5x3; NOT interchangeable',
        analytic_annulus_volume_m3=analytic, sampled_annulus_volume_relative_error=abs(ring[0]-analytic)/analytic,
        sampled_vs_frozen_ring_metric_relative=abs(ring[1]/ring[0]-original['ring_abs_p_relative_mean'])/original['ring_abs_p_relative_mean'],
        regions=rows)
    reports.append(report)
    print(run, 'samples', count, 'diagnostic_ring',ring[1]/ring[0], 'frozen_ring',original['ring_abs_p_relative_mean'],flush=True)

result = dict(scope='EXISTING_FIELDS_DIAGNOSTIC_ONLY_NOT_ACCEPTANCE', created=time.strftime('%Y-%m-%dT%H:%M:%S%z'),
    script_sha256=sha(Path(__file__)),elapsed_seconds=time.perf_counter()-started,runs=reports,
    limitations=['No matched-point differences across meshes; only weighted regional diagnostics.',
    'No new quadrature or PDE solve; archived 5x5x5 samples cannot replace frozen 5x5x3 ring metric.',
    'Mixed xy/z refinement does not uniquely identify directional discretization errors.',
    'Frozen 2% mesh and 1% stabilization gates remain unchanged.'])
with (OUT/'diagnostic.json').open('x',encoding='utf-8') as f: json.dump(result,f,indent=2,allow_nan=False)
with (OUT/'regional_metrics.csv').open('x',newline='',encoding='utf-8') as f:
    keys=['run']+list(reports[0]['regions'][0]); writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader()
    for r in reports:
        for row in r['regions']: writer.writerow(dict(run=r['run'],**row))
print('elapsed_seconds',result['elapsed_seconds'],flush=True)
