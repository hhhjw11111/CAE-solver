"""Transfer a HEX20 field ONLY as an initial guess; never as a computed result.

Inverse isoparametric mapping with a bounded extension for the slightly different
quadratic representations of the same circular wall. Every target run must solve
its own equations and repeat all acceptance checks.
"""
from pathlib import Path
import argparse,csv,json
import numpy as np
from audit_hex20_boundary_geometry import read_mesh,NATURAL,derivatives
from run_hv_sweep import sha

FIELDS=['ux_m','uy_m','uz_m','phi_V','eta_n','eta_p','theta_K']

def shapes(points):
    points=np.asarray(points);f=1+points[:,None,:]*NATURAL[None,:,:]
    n=np.zeros((len(points),20))
    n[:,:8]=np.prod(f[:,:8],axis=2)*(np.sum(points[:,None,:]*NATURAL[None,:8,:],axis=2)-2)/8
    for a in range(8,20):
        z=int(np.flatnonzero(NATURAL[a]==0)[0]);other=[k for k in range(3) if k!=z]
        n[:,a]=(1-points[:,z]**2)*f[:,a,other[0]]*f[:,a,other[1]]/4
    return n

def transfer(coords,conn,values,targets,extension=2e-9):
    elements=coords[conn];lo=elements.min(axis=1)-extension;hi=elements.max(axis=1)+extension
    result=np.empty((len(targets),values.shape[1]));outside=0;worst=0.;max_residual=0.
    for i,point in enumerate(targets):
        ids=np.flatnonzero(np.all((point>=lo)&(point<=hi),axis=1))
        if not len(ids):raise ValueError(f'No source element near target node {i+1}')
        x=elements[ids];xi=np.zeros((len(ids),3))
        for _ in range(12):
            r=np.einsum('pa,paj->pj',shapes(xi),x)-point
            jac=np.einsum('pka,paj->pjk',derivatives(xi),x)
            step=np.linalg.solve(jac,r[...,None])[...,0]
            xi-=np.clip(step,-1.,1.)
            if np.any((np.max(abs(r),axis=1)<2e-17)&(np.max(abs(xi),axis=1)<=1+1e-10)):break
        r=np.einsum('pa,paj->pj',shapes(xi),x)-point
        error=np.linalg.norm(r,axis=1);excess=np.maximum(abs(xi)-1,0).max(axis=1)
        valid=(error<1e-14)&(excess<.3)
        if not valid.any():raise ValueError(f'Inverse mapping failed at target node {i+1}: {error.min()}')
        score=np.where(valid,excess+error/1e-6,np.inf);j=int(np.argmin(score))
        # Quantify the physical distance outside the selected source element.
        clipped=np.clip(xi[j],-1,1)
        outside_distance=float(np.linalg.norm(shapes(clipped[None,:])[0]@x[j]-point))
        if outside_distance>extension+1e-14:raise ValueError(f'Extrapolation too large at node {i+1}')
        outside+=int(excess[j]>1e-9);worst=max(worst,outside_distance);max_residual=max(max_residual,float(error[j]))
        result[i]=shapes(xi[j:j+1])[0]@values[conn[ids[j]]]
        if i>0 and i%2000==0:print(f'TRANSFER_INITIAL_GUESS {i}/{len(targets)}',flush=True)
    if not np.isfinite(result).all():raise ValueError('Nonfinite transferred field')
    return result,dict(target_nodes=len(targets),boundary_extension_nodes=outside,
        maximum_boundary_extension_m=worst,maximum_inverse_mapping_residual_m=max_residual)

def main():
    p=argparse.ArgumentParser();p.add_argument('--source-mesh',type=Path,required=True)
    p.add_argument('--source-nodes',type=Path,required=True);p.add_argument('--target-mesh',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    coords,conn,_=read_mesh(a.source_mesh);target,_,_=read_mesh(a.target_mesh)
    data=np.genfromtxt(a.source_nodes,delimiter=',',names=True)
    if len(data)!=len(coords):raise ValueError('Source mesh/field count mismatch')
    if not np.allclose(np.column_stack([data[k] for k in ['x_m','y_m','z_m']]),coords,atol=1e-16,rtol=0):
        raise ValueError('Source mesh/field coordinates mismatch')
    vals=np.column_stack([data[k] for k in FIELDS]);values,report=transfer(coords,conn,vals,target)
    n0=float(np.median(data['n_m3']/np.exp(data['eta_n'])))
    p0=float(np.median(data['p_m3']/np.exp(data['eta_p'])))
    a.output.mkdir(parents=True,exist_ok=False)
    with (a.output/'nodes.csv').open('x',newline='',encoding='ascii') as f:
        w=csv.writer(f);w.writerow(['node','x_m','y_m','z_m']+FIELDS+['n_m3','p_m3'])
        for i,(xyz,val) in enumerate(zip(target,values),1):
            w.writerow([i]+list(xyz)+list(val)+[n0*np.exp(val[4]),p0*np.exp(val[5])])
    report.update(status='INITIAL_GUESS_ONLY_NOT_A_SOLVED_FIELD',source_mesh_sha256=sha(a.source_mesh),
        source_nodes_sha256=sha(a.source_nodes),target_mesh_sha256=sha(a.target_mesh))
    (a.output/'transfer_report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)

if __name__=='__main__':main()
