"""Volume-weighted directional cell Peclet diagnostics on saved JE fields."""
from pathlib import Path
import argparse,json
import numpy as np
import diagnose_je_power_gap_20261015 as audit


def run(directory):
    directory=Path(directory).resolve();manifest=json.loads((directory/'manifest.json').read_text())
    xyz,conn=audit.fe.read_mesh(Path(manifest['mesh']))
    state=np.loadtxt(directory/'nodes.csv',delimiter=',',skiprows=1)[:,4:11]
    rr=audit.rule(3);peclet=[];weights=[]
    for nodes in conn:
        a=audit.sample(xyz[nodes],state[nodes],rr)
        T=300+rr[2]@state[nodes,6];mu=audit.fe.MU_N0*(T/300)**-2
        D=audit.fe.KB*T/audit.fe.CHARGE*mu;v=mu[:,None]*a['E']
        ji=np.linalg.inv(a['jac']);metric=np.einsum('pia,pja->pij',ji,ji)
        trace=np.trace(metric,axis1=1,axis2=2);speed2=np.sum(v*v,axis=1)
        av=np.einsum('pi,pij,pj->p',v,metric,v);eps=1e-12*D*D*trace
        gdir=(av+eps*trace/3)/(speed2+eps)
        pe=np.sqrt(speed2)/(D*np.sqrt(gdir))
        peclet.extend(pe);weights.extend(rr[1]*a['det'])
    pe=np.array(peclet);w=np.array(weights);idx=np.argsort(pe);cum=np.cumsum(w[idx])/w.sum()
    report={'definition':'h_parallel=2/sqrt(gdir); Pe=|mu_n E| h_parallel/(2 D_n); gdir matches native directional metric regularization',
            'run':directory.name,'quadrature_order':3,'max':float(pe.max()),
            'volume_weighted_mean':float(w@pe/w.sum()),
            'volume_weighted_p50':float(np.interp(.5,cum,pe[idx])),
            'volume_weighted_p95':float(np.interp(.95,cum,pe[idx])),
            'volume_fraction_over':{str(t):float(w[pe>t].sum()/w.sum()) for t in [1,10,100]},
            'mesh_sha256':audit.fe.sha256(Path(manifest['mesh'])),'nodes_sha256':audit.fe.sha256(directory/'nodes.csv'),
            'limitation':'Mesh-dependent diagnostic, not a separate convergence or positivity certificate'}
    with (directory/'directional_pe.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--runs',nargs='+',type=Path,required=True);a=p.parse_args()
    for d in a.runs:run(d)
