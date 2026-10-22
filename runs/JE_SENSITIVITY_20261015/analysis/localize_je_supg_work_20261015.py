"""Localize phi-weighted SUPG work on an already solved field, using assembly quadrature."""
from pathlib import Path
import argparse,json
import numpy as np
import diagnose_je_power_gap_20261015 as audit


def run(directory):
    directory=Path(directory).resolve()
    m=json.loads((directory/'manifest.json').read_text());scale=m['supg_scale']
    xyz,conn=audit.fe.read_mesh(Path(m['mesh']))
    state=np.loadtxt(directory/'nodes.csv',delimiter=',',skiprows=1)[:,4:11]
    rr=audit.rule(3);total=0.;positive=0.;negative=0.;direction=np.zeros(3)
    zones={k:0. for k in ['z_surface_0_10nm','z_surface_10_50nm','z_surface_50_100nm',
                         'z_interior_over_100nm','hole_distance_under_0.5um','contact_distance_under_0.5um']}
    for nodes in conn:
        a=audit.sample(xyz[nodes],state[nodes],rr,True)
        w=rr[1]*a['det'];work=scale*np.einsum('pn,n->p',a['S'],state[nodes,3])
        total+=float(w@work)
        positive+=float(w@np.maximum(work,0));negative+=float(w@np.minimum(work,0))
        # v_n and v_p are parallel to E in this isotropic material. Their
        # sum S_a = C * grad(N_a).E; phi*S = -C*|E|^2, split by E_i^2.
        e2=np.sum(a['E']**2,axis=1)
        direction+=np.sum((w*work)[:,None]*a['E']**2/np.maximum(e2[:,None],1e-300),axis=0)
        depth=.5e-6-abs(a['xyz'][:,2]);r=np.linalg.norm(a['xyz'][:,:2],axis=1)
        masks=[depth<10e-9,(depth>=10e-9)&(depth<50e-9),(depth>=50e-9)&(depth<100e-9),
               depth>=100e-9,r<3.5e-6,abs(a['xyz'][:,0])>9.5e-6]
        for key,mask in zip(zones,masks):zones[key]+=float(w[mask]@work[mask])
    assert abs(direction.sum()-total)<1e-12
    assert abs(sum(list(zones.values())[:4])-total)<1e-12
    report={'run':directory.name,'supg_scale':scale,'weighted_SUPG_W':total,
            'positive_local_work_W':positive,'negative_local_work_W':negative,
            'absolute_local_work_W':positive-negative,
            'signed_to_absolute_ratio':total/(positive-negative),
            'direction_x_y_z_W':direction.tolist(),'direction_x_y_z_percent':(direction/total*100).tolist(),
            'zones_W':zones,'zones_percent_of_signed_total':{k:v/total*100 for k,v in zones.items()},
            'note':'Signed contributions; hole and contact zones overlap the thickness partition. Direction is an algebraic work decomposition, not an independent mesh-error estimator.',
            'mesh_sha256':audit.fe.sha256(Path(m['mesh'])),'nodes_sha256':audit.fe.sha256(directory/'nodes.csv'),
            'script_sha256':audit.fe.sha256(Path(__file__))}
    with (directory/'SUPG_localization_v2.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--runs',nargs='+',type=Path,required=True);a=p.parse_args()
    for d in a.runs:run(d)
