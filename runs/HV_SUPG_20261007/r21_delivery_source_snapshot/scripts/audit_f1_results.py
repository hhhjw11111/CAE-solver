"""Frozen comparisons on saved F1 results. Never invokes the FEM solver."""
from pathlib import Path
import argparse,json,math
import numpy as np

ROOT=Path(__file__).resolve().parents[1]/'runs/F1_20260906'

def read_metrics(case):
    path=ROOT/case/'metrics.txt';out={}
    for line in path.read_text(encoding='utf-8-sig').splitlines():
        if '=' not in line or line.startswith('stage='):continue
        k,v=line.split('=',1)
        try:out[k]=float(v)
        except ValueError:out[k]=v
    return out

def rel(a,b):return abs(a-b)/max(abs(a),abs(b),1e-300)

def normdiff(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-300))

def small(a,b):
    ma,mb=read_metrics(a),read_metrics(b)
    # 0.14.0 appended one empty trailing CSV field; the 13 labelled fields are intact.
    da=np.genfromtxt(ROOT/a/'nodes.csv',delimiter=',',names=True,usecols=range(13))
    db=np.genfromtxt(ROOT/b/'nodes.csv',delimiter=',',names=True,usecols=range(13))
    assert len(da)==len(db)
    assert all(np.array_equal(da[c],db[c]) for c in ['node','x_m','y_m','z_m'])
    ua=np.column_stack([da[c] for c in ['ux_m','uy_m','uz_m']])
    ub=np.column_stack([db[c] for c in ['ux_m','uy_m','uz_m']])
    n0=6e22;p0=(2.9e-4)**2/n0
    out={'cases':[a,b],'sampling':'matched physical nodal fields, same mesh; physical densities divided by n0/p0',
         'phi_relative_L2_nodes':normdiff(da['phi_V'],db['phi_V']),
         'u_relative_L2_nodes':normdiff(ua,ub),
         'n_reference_scaled_difference_RMS':float(np.sqrt(np.mean(((da['n_m3']-db['n_m3'])/n0)**2))),
         'p_reference_scaled_difference_RMS':float(np.sqrt(np.mean(((da['p_m3']-db['p_m3'])/p0)**2))),
         'maximum_density_shift':max(*(m['maximum_abs_n_relative'] for m in [ma,mb]),
            *(m['maximum_abs_p_relative'] for m in [ma,mb]),
            *(float(np.max(np.abs(d['n_m3']/n0-1))) for d in [da,db]),
            *(float(np.max(np.abs(d['p_m3']/p0-1))) for d in [da,db])),
         'Q_over_V2_relative_difference':rel(ma['joule_W']/ma['voltage_V']**2,mb['joule_W']/mb['voltage_V']**2),
         'theta_over_V2_relative_difference':normdiff(da['theta_K']/ma['voltage_V']**2,db['theta_K']/mb['voltage_V']**2),
         'density_degeneracy_note':'For this c||z, x-driven, traction-free extruded geometry, the isothermal first-order LQ density change is zero to roundoff. Do not divide by its near-zero perturbation norm or claim a finite first-order carrier response.'}
    fields=['phi_relative_L2_nodes','u_relative_L2_nodes','n_reference_scaled_difference_RMS','p_reference_scaled_difference_RMS']
    out['passed']=all(out[k]<.01 for k in fields) and out['maximum_density_shift']<=.01
    out['passed']=out['passed'] and all(m.get('case_acceptance')=='PASS_NOT_MESH_ACCEPTANCE' for m in [ma,mb])
    return out

def mesh(cases):
    data=[read_metrics(c) for c in cases]
    limits={'joule_W':.01,'theta_max_K':.01,'ring_theta_mean_K':.01,
            'ring_E_mean_V_m':.02,'ring_Q_mean_W_m3':.02,'ring_abs_n_relative_mean':.02,
            'ring_Q_p95_W_m3':.05,'E_p95_V_m':.05}
    supplements=[ROOT/c/'supplemental_metrics.json' for c in cases]
    if any(p.exists() for p in supplements):
        if not all(p.exists() for p in supplements):
            raise ValueError('Incomplete supplemental mesh metrics')
        for d,p in zip(data,supplements):
            d['ring_E_p95_V_m']=json.loads(p.read_text(encoding='utf-8-sig'))['ring_E_p95_V_m']
        limits['ring_E_p95_V_m']=.05
    changes={k:{'L0_L1':rel(data[0][k],data[1][k]),'L1_L2':rel(data[1][k],data[2][k]),'limit':lim}
             for k,lim in limits.items()}
    return {'cases':cases,'metric_changes':changes,
            'passed':all(v['L1_L2']<v['limit'] for v in changes.values()) and
                     all(d.get('case_acceptance')=='PASS_NOT_MESH_ACCEPTANCE' for d in data)}

def main():
    global ROOT
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['small','mesh'])
    p.add_argument('--run-root',type=Path,default=ROOT)
    p.add_argument('cases',nargs='*');args=p.parse_args()
    ROOT=args.run_root.resolve()
    if args.mode=='small':out=small(*(args.cases or ['SMALL_H1_LQ_L0','SMALL_H1_NT_L0']))
    else:out=mesh(args.cases or [f'MESH_H1_NT_L{i}' for i in range(3)])
    dest=ROOT/(args.mode+'_comparison.json')
    dest.write_text(json.dumps(out,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(out,indent=2,allow_nan=False))
    raise SystemExit(0 if out['passed'] else 1)

if __name__=='__main__':main()
