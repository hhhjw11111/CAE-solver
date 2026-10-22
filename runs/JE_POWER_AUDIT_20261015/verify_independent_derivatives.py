"""Independent finite-difference check of current divergence on curved HEX20 cells."""
from pathlib import Path
import sys
import json
import numpy as np

SOLVER=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(SOLVER/'scripts'))
import diagnose_je_power_gap_20261015 as audit
import plot_je_8v_clouds_20261013 as fe


def at(points):
    dn=fe.derivatives(points)
    hn=np.empty((len(points),3,3,20))
    for axis in range(3):
        shift=np.eye(3)[axis]
        hn[:,:,axis,:]=(fe.derivatives(points+shift)-fe.derivatives(points-shift))/2
    return points,np.ones(len(points)),fe.shapes(points),dn,hn


xyz,conn=fe.read_mesh(audit.MESH)
nodes=np.loadtxt(audit.ROOT/'V080/nodes.csv',delimiter=',',skiprows=1)
state=nodes[:,4:11]
ids=[0,17,223,224,467,895,1200,1567,2000,2687,3000,3135]
points=np.array([[-.37,.21,.13],[.47,-.29,-.41],[.07,.61,.31]])
reports=[]
for step in [1e-3,1e-4]:
    errs=[];element_errors=[]
    for e in ids:
        ids_e=conn[e];coords=xyz[ids_e];xe=state[ids_e]
        base=audit.sample(coords,xe,at(points),True)
        dJ_nat=np.empty((len(points),3,3))
        for axis in range(3):
            shift=np.eye(3)[axis]*step
            plus=audit.sample(coords,xe,at(points+shift))['J']
            minus=audit.sample(coords,xe,at(points-shift))['J']
            dJ_nat[:,axis,:]=(plus-minus)/(2*step)
        dJ_phys=np.einsum('pia,paj->pij',np.linalg.inv(base['jac']),dJ_nat)
        fd=np.trace(dJ_phys,axis1=1,axis2=2)
        scales=np.maximum(np.linalg.norm(dJ_phys,axis=(1,2)),1.)
        err=abs(fd-base['divJ'])/scales
        errs.extend(err.tolist())
        element_errors.append({'element':e+1,'max_gradient_normalized_error':float(err.max())})
    reports.append({'natural_coordinate_step':step,'points':len(ids)*len(points),
                    'max_gradient_normalized_error':max(errs),'elements':element_errors})
passed=reports[1]['max_gradient_normalized_error']<1e-5
out=Path(__file__).with_name('finite_difference_attempt.json')
out.write_text(json.dumps({'status':'PASS' if passed else 'FAIL_STEP_SENSITIVITY','definition':'Error divided by Frobenius norm of independently differentiated physical current gradient; no floor on div J itself.','reports':reports},indent=2),encoding='utf-8')
print(json.dumps([{k:v for k,v in r.items() if k!='elements'} for r in reports],indent=2))
assert passed,'Finite difference small-step gate failed; retain this attempt.'
