"""Check solved full/half symmetric states, with gauge-invariant displacement.

This checks equivalence of the sampled symmetric branch, not mesh accuracy,
asymmetric solutions, dynamical stability, or the full 0--8 V interval.
"""
import argparse,json
from pathlib import Path
import numpy as np
from generate_hv_half_y import rigid_gauge, U
from audit_hex20_boundary_geometry import read_mesh
from run_hv_sweep import classify,sha

FIELDS=U+['phi_V','eta_n','eta_p','theta_K']
LIMIT=1e-8
EXTENSIVE=['joule_W','electron_right_A','hole_right_A','total_volume_m3','ring_volume_m3']
INTENSIVE=['theta_max_K','ring_theta_mean_K','ring_E_mean_V_m','ring_Q_mean_W_m3',
           'ring_abs_n_relative_mean','ring_abs_p_relative_mean','stress_von_mises_max_Pa']


def audit(full,half):
    errors=[];manifest=[];metrics=[]
    for folder,domain in [(full,'FULL'),(half,'Y_HALF')]:
        m=json.loads((folder/'manifest.json').read_text());manifest.append(m)
        ex=json.loads((folder/'execution.json').read_text())
        state,physical=classify(folder,ex['exit_code']);metrics.append(physical)
        if state!='SOLVE_FLUX_PASS_MESH_PENDING':errors.append(domain+':SOLVE_OR_FLUX_FAIL')
        if json.loads((folder/'field_audit.json').read_text()).get('passed') is not True:
            errors.append(domain+':FIELD_AUDIT_FAIL')
        if m.get('domain_symmetry','FULL')!=domain or physical.get('domain_symmetry','FULL')!=domain:
            errors.append(domain+':DOMAIN_MISMATCH')
        for key in ['mesh','executable']:
            if sha(Path(m[key]))!=m[key+'_sha256']:errors.append(domain+':'+key.upper()+'_HASH_MISMATCH')
    for key in ['executable_sha256','supg_scale','model']:
        if manifest[0].get(key)!=manifest[1].get(key):errors.append('INCONSISTENT_'+key.upper())
    for key in ['voltage_V','model','heat_source','source','supg_scale','supg_directional_tau','transport_discretization']:
        if metrics[0].get(key)!=metrics[1].get(key):errors.append('INCONSISTENT_'+key.upper())
    if manifest[0].get('log_weighted_quadrature_order',3)!=manifest[1].get('log_weighted_quadrature_order',3):
        errors.append('INCONSISTENT_QUADRATURE_ORDER')
    fullxyz,_,_=read_mesh(Path(manifest[0]['mesh']));halfxyz,_,_=read_mesh(Path(manifest[1]['mesh']))
    # Restricted nodes retain their exact floating-point coordinates.
    lookup={tuple(x):i for i,x in enumerate(fullxyz)}
    ids=np.array([lookup[tuple(x)] for x in halfxyz]);distance=np.linalg.norm(fullxyz[ids]-halfxyz,axis=1)
    if distance.max()>1e-16 or len(np.unique(ids))!=len(ids):errors.append('HALF_NODE_MAPPING_ERROR')
    a=np.genfromtxt(full/'nodes.csv',delimiter=',',names=True);b=np.genfromtxt(half/'nodes.csv',delimiter=',',names=True)
    fullu=np.column_stack([a[k] for k in U]);fullu,gauge=rigid_gauge(fullxyz,fullu)
    scales=[1e-10]*3+[max(1.,metrics[0]['voltage_V']),.1,.1,35.]
    comparisons={}
    for j,(key,scale) in enumerate(zip(FIELDS,scales)):
        av=fullu[ids,j] if j<3 else a[key][ids]
        diff=float(np.max(abs(av-b[key]))/scale)
        comparisons[key]=dict(maximum_scaled_difference=diff,scale=scale,limit=LIMIT,passed=bool(diff<LIMIT))
    physical={}
    for key in EXTENSIVE+INTENSIVE:
        av=metrics[0][key];bv=metrics[1][key]*(2 if key in EXTENSIVE else 1)
        diff=abs(av-bv)/max(abs(av),abs(bv),1e-300)
        physical[key]=dict(full=av,half_reconstructed_full=bv,relative_difference=diff,limit=LIMIT,passed=diff<LIMIT)
    for key in ['heat_out_W','terminal_power_W']:
        av=metrics[0][key];bv=2*metrics[1].get('physical_exterior_'+key,float('nan'))
        diff=abs(av-bv)/max(abs(av),abs(bv),1e-300)
        physical[key]=dict(full=av,half_reconstructed_full=bv if np.isfinite(bv) else None,
            relative_difference=diff if np.isfinite(diff) else None,limit=LIMIT,passed=bool(diff<LIMIT))
    passed=not errors and all(x['passed'] for x in comparisons.values()) and all(x['passed'] for x in physical.values())
    return dict(scope='SAMPLED_SYMMETRIC_BRANCH_EQUIVALENCE',full_run=str(full.absolute()),half_run=str(half.absolute()),
        errors=errors,passed=passed,nodal=comparisons,physical=physical,gauge=gauge,
        mesh_validated=False,complete_quantitative_0_8V=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--full',type=Path,required=True);p.add_argument('--half',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=audit(a.full,a.half)
    with a.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps(result,indent=2,allow_nan=False));return 0 if result['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
