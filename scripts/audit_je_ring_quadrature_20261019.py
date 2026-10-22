"""Read-only quadrature diagnostic on fixed states; preserve native FAIL gates."""
from pathlib import Path
import argparse
import json
import numpy as np
import run_je_reconstruction_mesh_20261018 as r
import plot_je_8v_clouds_20261013 as fe

RULES=((5,5,3),(5,5,5),(7,7,7),(9,9,9))


def rule(orders):
    rules=[np.polynomial.legendre.leggauss(n) for n in orders]
    points=np.array([(x,y,z) for z in rules[2][0] for y in rules[1][0] for x in rules[0][0]])
    weights=np.array([a*b*c for c in rules[2][1] for b in rules[1][1] for a in rules[0][1]])
    return points,weights


def integrate(xyz,conn,eta,orders,inner,outer):
    pts,weights=rule(orders);shape=fe.shapes(pts);dn=fe.derivatives(pts)
    total=np.zeros(3);sample_count=0;minimum_det=float('inf')
    for start in range(0,len(conn),64):
        ids=conn[start:start+64];coords=xyz[ids]
        pos=np.einsum('qa,eaj->eqj',shape,coords)
        jac=np.einsum('qka,eaj->eqkj',dn,coords)
        det=np.linalg.det(jac);minimum_det=min(minimum_det,float(det.min()))
        if np.any(det<=0):raise ValueError('Nonpositive sampled Jacobian')
        radius=np.linalg.norm(pos[:,:,:2],axis=2)
        mask=(radius>=inner)&(radius<=outer)
        vals=np.einsum('qa,eaf->eqf',shape,eta[ids])
        with np.errstate(over='raise',invalid='raise'):
            perturbation=abs(np.expm1(vals))
        if not np.isfinite(perturbation).all():raise ValueError('Nonfinite carrier value')
        w=det*weights[None,:]*mask
        total[0]+=float(w.sum())
        total[1:]+=np.einsum('eq,eqf->f',w,perturbation)
        sample_count+=int(mask.sum())
    if total[0]<=0:raise ValueError('Empty annulus')
    return dict(orders=list(orders),ring_volume_m3=float(total[0]),
        ring_abs_n_relative_mean=float(total[1]/total[0]),
        ring_abs_p_relative_mean=float(total[2]/total[0]),
        annulus_sample_count=sample_count,all_domain_sample_count=len(conn)*len(pts),
        minimum_jacobian_m3=minimum_det)


def check():
    lengths=np.array([2e-6,4e-6,1e-6]);xyz=fe.NATURAL*lengths/2
    eta=np.tile([.13,-.17],(20,1));expected=abs(np.expm1(eta[0]))
    results=[]
    for orders in RULES:
        d=integrate(xyz,np.arange(20)[None,:],eta,orders,0,10e-6)
        np.testing.assert_allclose(d['ring_volume_m3'],np.prod(lengths),rtol=2e-12)
        np.testing.assert_allclose([d['ring_abs_n_relative_mean'],d['ring_abs_p_relative_mean']],expected,rtol=2e-12)
        results.append(d)
    return dict(status='CONSTANT_FIELD_BOX_CHECK_PASSED',rules=results)


def audit(run,output):
    manifest=r.accepted(run);_,metrics=r.classify(run,r.read(run/'execution.json')['exit_code'])
    assert metrics['model']=='NT' and metrics['hole_radius_m']==3e-6
    xyz,conn=fe.read_mesh(Path(manifest['mesh']))
    nodes=np.loadtxt(run/'nodes.csv',delimiter=',',skiprows=1)
    assert np.max(abs(nodes[:,1:4]-xyz))<1e-19
    files=[run/n for n in ('nodes.csv','input.nml','manifest.json','execution.json','metrics.txt')]
    files += [Path(manifest['mesh']),Path(__file__),Path(fe.__file__),r.ROOT/'src/post/field_metrics_m.f90']
    hashes={str(p):r.sha(p) for p in files}
    results=[]
    for orders in RULES:
        d=integrate(xyz,conn,nodes[:,8:10],orders,3e-6,3.5e-6)
        results.append(d)
        if orders==(5,5,3):
            for key in ('ring_volume_m3','ring_abs_n_relative_mean','ring_abs_p_relative_mean'):
                assert abs(d[key]-metrics[key])/abs(metrics[key])<5e-11,(key,d[key],metrics[key])
        print(json.dumps(dict(run=run.name,**d)),flush=True)
    assert all(r.sha(Path(p))==h for p,h in hashes.items())
    reference=results[0]
    for d in results:
        d['p_change_vs_native_fraction']=abs(d['ring_abs_p_relative_mean']-reference['ring_abs_p_relative_mean'])/reference['ring_abs_p_relative_mean']
    report=dict(run=str(run),voltage_V=metrics['voltage_V'],native_rule_reproduced=True,
        rules=results,source_hashes=hashes,production_acceptance=False,
        limitation='Fixed-field postprocessing only. Native metrics, mesh thresholds and existing FAIL records are unchanged.')
    r.save(output/(run.name+'.json'),report)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs',nargs='+',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    r.save(args.output/'closed_form_check.json',check())
    reports=[audit(p.resolve(),args.output) for p in args.runs]
    r.save(args.output/'summary.json',dict(scope='FIXED_STATE_RING_QUADRATURE_DIAGNOSTIC',
        runs=[x['run'] for x in reports],native_rule=[5,5,3],
        production_acceptance=False,thresholds_changed=False,
        maximum_native_to_order9_p_change_fraction=max(x['rules'][-1]['p_change_vs_native_fraction'] for x in reports)))


if __name__=='__main__':main()
