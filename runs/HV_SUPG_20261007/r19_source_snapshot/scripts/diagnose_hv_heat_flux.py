"""Read-only independent boundary heat-flux localization for frozen ZnO k=60.

Uses no assembled reactions. This diagnostic cannot grant case acceptance.
"""
import argparse,json
from pathlib import Path
import numpy as np
from audit_hex20_boundary_geometry import read_mesh,NATURAL,derivatives
from run_hv_sweep import read_metrics,sha

def diagnose(mesh,nodes,metrics):
    coords,conn,_=read_mesh(mesh);data=np.genfromtxt(nodes,delimiter=',',names=True)
    if len(data)!=len(coords):raise ValueError('Node count mismatch')
    theta=data['theta_K'];faces={}
    for e,c in enumerate(conn):
        for axis in range(3):
            for sign in [-1,1]:
                ids=np.flatnonzero(NATURAL[:,axis]==sign);key=tuple(sorted(c[ids]))
                if key in faces:faces.pop(key)
                else:faces[key]=(e,axis,sign,ids)
    gp,gw=np.polynomial.legendre.leggauss(5);lo=coords.min(axis=0);hi=coords.max(axis=0)
    groups={k:dict(faces=0,heat_out_W=0.,absolute_heat_out_W=0.,area_m2=0.) for k in
        ['x_min_contact','x_max_contact','y_min','y_max','z_min','z_max','hole']}
    for e,axis,sign,ids in faces.values():
        x=coords[conn[e]];face=x[ids];group='hole'
        for k,label in enumerate(['x','y','z']):
            for bound,suffix in [(lo[k],'min'),(hi[k],'max')]:
                if np.max(abs(face[:,k]-bound))<1e-12:
                    group=label+'_'+suffix+('_contact' if k==0 else '')
        aa=(axis+1)%3;bb=(axis+2)%3;p=np.zeros((25,3));weights=[]
        for n,(i,j) in enumerate((i,j) for i in range(5) for j in range(5)):
            p[n,axis]=sign;p[n,aa]=gp[i];p[n,bb]=gp[j];weights.append(gw[i]*gw[j])
        dn=derivatives(p);jac=np.einsum('pka,aj->pkj',dn,x)
        if np.any(np.linalg.det(jac)<=0):raise ValueError('Nonpositive face Jacobian')
        grad=np.einsum('pjk,pka->pja',np.linalg.inv(jac),dn)
        heat=-60*np.einsum('pja,a->pj',grad,theta[conn[e]])
        area=sign*np.cross(jac[:,aa,:],jac[:,bb,:])*np.array(weights)[:,None]
        normal=np.einsum('pj,pj->p',heat,area);r=groups[group]
        r['faces']+=1;r['heat_out_W']+=float(normal.sum());r['absolute_heat_out_W']+=float(abs(normal).sum())
        r['area_m2']+=float(np.linalg.norm(area,axis=1).sum())
    total=sum(r['heat_out_W'] for r in groups.values());m=read_metrics(metrics)
    result=dict(status='DIAGNOSTIC_ONLY',groups=groups,total_heat_out_W=total,mesh_sha256=sha(mesh),nodes_sha256=sha(nodes))
    if 'heat_out_W' in m:
        result['difference_from_fortran_W']=total-m['heat_out_W']
        result['joule_W']=m['joule_W'];result['heat_balance_relative']=abs(total-m['joule_W'])/max(abs(total),abs(m['joule_W']),1e-15)
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--mesh',required=True,type=Path)
    p.add_argument('--nodes',required=True,type=Path);p.add_argument('--metrics',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    r=diagnose(a.mesh,a.nodes,a.metrics)
    with a.output.open('x') as f:json.dump(r,f,indent=2)
    print(json.dumps(r,indent=2))

if __name__=='__main__':main()
