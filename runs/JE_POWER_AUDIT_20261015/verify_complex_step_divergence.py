"""Subtraction-free independent differentiation of physical current on curved cells."""
from pathlib import Path
import sys,json
import numpy as np
SOLVER=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(SOLVER/'scripts'))
import diagnose_je_power_gap_20261015 as audit
import plot_je_8v_clouds_20261013 as fe


def polynomial(points):
    # Explicit serendipity polynomials evaluated in complex arithmetic.
    p=np.asarray(points,dtype=complex);natural=fe.NATURAL
    fac=1+p[:,None,:]*natural[None,:,:]
    linear=np.sum(p[:,None,:]*natural[None,:,:],axis=2)-2
    n=np.empty((len(p),20),complex);dn=np.zeros((len(p),3,20),complex)
    n[:,:8]=np.prod(fac[:,:8,:],axis=2)*linear[:,:8]/8
    for axis in range(3):
        other=[a for a in range(3) if a!=axis]
        dn[:,axis,:8]=natural[None,:8,axis]/8*fac[:,:8,other[0]]*fac[:,:8,other[1]]*(linear[:,:8]+fac[:,:8,axis])
    for k in range(8,20):
        zero=int(np.flatnonzero(natural[k]==0)[0]);other=[a for a in range(3) if a!=zero]
        n[:,k]=(1-p[:,zero]**2)*fac[:,k,other[0]]*fac[:,k,other[1]]/4
        dn[:,zero,k]=-p[:,zero]*fac[:,k,other[0]]*fac[:,k,other[1]]/2
        for a in other:
            b=next(x for x in other if x!=a)
            dn[:,a,k]=(1-p[:,zero]**2)*natural[k,a]*fac[:,k,b]/4
    return n,dn


def current(coords,state,points):
    n,dn=polynomial(points);jac=np.einsum('pan,ni->pai',dn,coords)
    grad=np.linalg.solve(jac,dn);vals=n@state;dv=np.einsum('pin,nf->pif',grad,state)
    E=-dv[:,:,3];T=300+vals[:,6];J=np.zeros_like(E)
    for field,ref,mu0,sign in [(4,fe.N0,fe.MU_N0,1),(5,fe.P0,fe.MU_P0,-1)]:
        rho=ref*np.exp(vals[:,field]);mu=mu0*(T/300)**-2;D=fe.KB*T/fe.CHARGE*mu
        J+=fe.CHARGE*rho[:,None]*(mu[:,None]*E+sign*D[:,None]*dv[:,:,field])
    return J,jac


def at(points):
    dn=fe.derivatives(points);hn=np.empty((len(points),3,3,20))
    for a in range(3):
        shift=np.eye(3)[a];hn[:,:,a,:]=(fe.derivatives(points+shift)-fe.derivatives(points-shift))/2
    return points,np.ones(len(points)),fe.shapes(points),dn,hn


xyz,conn=fe.read_mesh(audit.MESH)
state=np.loadtxt(audit.ROOT/'V080/nodes.csv',delimiter=',',skiprows=1)[:,4:11]
ids=[0,17,223,224,467,895,1200,1567,2000,2687,3000,3135]
points=np.array([[-.37,.21,.13],[.47,-.29,-.41],[.07,.61,.31]])
reports=[]
for step in [1e-15,1e-20]:
    errs=[];jerrs=[]
    for e in ids:
        coords=xyz[conn[e]];xe=state[conn[e]];base=audit.sample(coords,xe,at(points),True)
        j0,jac=current(coords,xe,points)
        jerrs.append(float(np.max(abs(j0.real-base['J']))/max(np.max(abs(base['J'])),1.)))
        dJ_nat=np.empty((len(points),3,3))
        for axis in range(3):
            shift=np.eye(3)[axis]*step*1j
            dJ_nat[:,axis,:]=current(coords,xe,points+shift)[0].imag/step
        dJ_phys=np.einsum('pia,paj->pij',np.linalg.inv(jac.real),dJ_nat)
        derivative=np.trace(dJ_phys,axis1=1,axis2=2)
        scale=np.maximum(np.linalg.norm(dJ_phys,axis=(1,2)),1.)
        errs.extend((abs(derivative-base['divJ'])/scale).tolist())
    reports.append({'imaginary_step':step,'points':len(ids)*len(points),'max_gradient_normalized_error':max(errs),'current_evaluator_peak_relative_difference':max(jerrs)})
passed=all(x['max_gradient_normalized_error']<1e-7 for x in reports)
Path(__file__).with_name('complex_step_validation.json').write_text(json.dumps({'status':'PASS' if passed else 'FAIL','reports':reports,'method':'Complex-step derivative of the complete current, including geometry inverse, exp(eta), and T-dependent mobility/diffusion; no analytic current divergence or shape Hessian used on reference side.','prior_attempt':'finite_difference_attempt.json retains failed small-step subtractive test.'},indent=2),encoding='utf-8')
print(json.dumps(reports,indent=2));assert passed
