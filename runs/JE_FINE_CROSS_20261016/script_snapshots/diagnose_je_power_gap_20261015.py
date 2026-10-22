"""Read-only electrical power audit of saved JE states; does not alter solver physics.

W - Pext = integral(phi div J) - sum_internal_faces integral(phi jump(J.n)).
Independently assemble Galerkin and SUPG continuity residuals to identify the
weighted stabilization contribution. All currents are conventional charge currents.
"""
from pathlib import Path
from functools import lru_cache
import argparse
import csv
import json
import shutil
import numpy as np
import plot_je_8v_clouds_20261013 as fe

ROOT = fe.SOLVER/'runs/JE_6_8_HALF_V_20261013'
MESH = fe.SOLVER/'cases/HV_SURFACE_20261008/XY0_Z2/H1_XY0_Z2.mesh'


@lru_cache(None)
def rule(order,axis=-1,side=0):
    gp,gw=np.polynomial.legendre.leggauss(order)
    if axis < 0:
        pts=np.array([[r,s,t] for t in gp for s in gp for r in gp])
        weights=np.array([a*b*c for c in gw for b in gw for a in gw])
    else:
        pts=np.zeros((order*order,3));pts[:,axis]=side
        a,b=(axis+1)%3,(axis+2)%3
        pts[:,a]=np.tile(gp,order);pts[:,b]=np.repeat(gp,order)
        weights=np.outer(gw,gw).ravel()
    dn=fe.derivatives(pts)
    hn=np.empty((len(pts),3,3,20))
    for b in range(3):
        plus=pts.copy();minus=pts.copy();plus[:,b]+=1.;minus[:,b]-=1.
        hn[:,:,b,:]=(fe.derivatives(plus)-fe.derivatives(minus))*.5
    return pts,weights,fe.shapes(pts),dn,hn


def sample(coords,state,rule_data,hessian=False):
    pts,weights,n,dn,hn=rule_data
    jac=np.einsum('pan,ni->pai',dn,coords)
    ji=np.linalg.inv(jac);det=np.linalg.det(jac)
    assert np.all(det>0)
    grad=np.einsum('pia,pan->pin',ji,dn)
    vals=n@state;dv=np.einsum('pin,nf->pif',grad,state)
    E=-dv[:,:,3];T=300+vals[:,6];gt=dv[:,:,6]
    ans={'xyz':n@coords,'phi':vals[:,3],'E':E,'jac':jac,'det':det,'grad':grad,
         'J':np.zeros_like(E),'divJ':np.zeros(len(pts)),'S':np.zeros((len(pts),20)),
         'G':np.zeros((len(pts),20))}
    if hessian:
        hx=np.einsum('pabn,nk->pabk',hn,coords)
        corrected=hn-np.einsum('pabk,pkn->pabn',hx,grad)
        hess=np.einsum('pia,pabn,pjb->pijn',ji,corrected,ji)
        hs=np.einsum('pijn,nf->pijf',hess,state)
        metric=np.einsum('pia,pja->pij',ji,ji)
        trace=np.trace(metric,axis1=1,axis2=2)
    for field,ref,mu0,sgn in [(4,fe.N0,fe.MU_N0,1.),(5,fe.P0,fe.MU_P0,-1.)]:
        rho=ref*np.exp(vals[:,field]);mu=mu0*(T/300)**-2
        diff=fe.KB*T/fe.CHARGE*mu;mut=-2*mu/T;dt=-diff/T
        g=dv[:,:,field];v=mu[:,None]*E
        J=fe.CHARGE*rho[:,None]*(v+sgn*diff[:,None]*g)
        ans['J']+=J
        ans['G']+=np.einsum('pin,pi->pn',grad,J)
        if hessian:
            divv=-mu*np.trace(hs[:,:,:,3],axis1=1,axis2=2)+mut*np.sum(E*gt,axis=1)
            strong=np.sum(v*g,axis=1)+divv+sgn*(diff*(np.trace(hs[:,:,:,field],axis1=1,axis2=2)+np.sum(g*g,axis=1))+dt*np.sum(gt*g,axis=1))
            divJ=fe.CHARGE*rho*strong
            ans['divJ']+=divJ
            speed2=np.sum(v*v,axis=1);advmetric=np.einsum('pi,pij,pj->p',v,metric,v)
            eps=1e-12*diff*diff*trace
            gdir=(advmetric+eps*trace/3)/(speed2+eps)
            tau=1/np.sqrt(advmetric+(4*diff*gdir)**2)
            adv=np.einsum('pin,pi->pn',grad,v)
            ans['S']+=sgn*(tau*divJ)[:,None]*adv
    ans['W']=np.sum(ans['J']*E,axis=1)
    return ans


def topology(xyz,conn):
    faces={}
    for e,nodes in enumerate(conn):
        for axis in range(3):
            for side in [-1,1]:
                ids=nodes[np.abs(fe.NATURAL[:,axis]-side)<.5]
                key=tuple(sorted(ids.tolist()))
                faces.setdefault(key,[]).append((e,axis,side))
    assert all(len(x) in [1,2] for x in faces.values())
    lookup={}
    for ids,owners in faces.items():
        if len(owners)==2:
            kind='internal'
        else:
            x=xyz[list(ids)]
            kind='hole'
            for axis,label in enumerate(['x','y','z']):
                for side,limit in [('-',xyz[:,axis].min()),('+',xyz[:,axis].max())]:
                    if np.max(np.abs(x[:,axis]-limit))<1e-18:kind=label+side
        for owner in owners:lookup[owner]=kind
    return faces,lookup


def audit_order(xyz,conn,state,lookup,order,out):
    N=len(xyz);G=np.zeros(N);S=np.zeros(N)
    vol=rule(order);w=vol[1]
    ext={kind:{'power_W':0.,'current_A':0.,'absolute_current_A':0.,'area_m2':0.,'faces':0}
         for kind in ['x-','x+','y-','y+','z-','z+','hole']}
    W=D=jump=divsum=0.;area_normal=np.zeros(3);elements=[]
    for e,nodes in enumerate(conn):
        coords=xyz[nodes];xe=state[nodes]
        a=sample(coords,xe,vol,True);weight=w*a['det']
        we=float(weight@a['W']);de=float(weight@(a['phi']*a['divJ']))
        W+=we;D+=de;divsum+=float(weight@a['divJ'])
        np.add.at(G,nodes,weight@a['G']);np.add.at(S,nodes,weight@a['S'])
        inter=external=0.
        for axis in range(3):
            for side in [-1,1]:
                rr=rule(order,axis,side);b=sample(coords,xe,rr)
                av=side*np.cross(b['jac'][:,(axis+1)%3,:],b['jac'][:,(axis+2)%3,:])*rr[1][:,None]
                flux=np.sum(b['J']*av,axis=1)
                surface=float(b['phi']@flux)
                kind=lookup[e,axis,side]
                if kind=='internal':inter+=surface
                else:
                    external+=surface;area_normal+=av.sum(axis=0)
                    t=ext[kind];t['power_W']-=surface;t['current_A']+=float(flux.sum())
                    t['absolute_current_A']+=float(np.abs(flux).sum());t['area_m2']+=float(np.linalg.norm(av,axis=1).sum());t['faces']+=1
        jump+=inter
        elements.append({'element':e+1,'center_x_um':float(coords[:,0].mean()*1e6),'center_y_um':float(coords[:,1].mean()*1e6),
                         'center_z_um':float(coords[:,2].mean()*1e6),'JE_W':we,'phi_divJ_W':de,
                         'internal_phi_flux_W':inter,'exterior_phi_flux_W':external,
                         'element_ibp_error_W':we+external+inter-de})
        if e%800==0:print(f'ORDER={order} ELEMENT={e}/{len(conn)}',flush=True)
    left=np.abs(xyz[:,0]-xyz[:,0].min())<1e-18
    right=np.abs(xyz[:,0]-xyz[:,0].max())<1e-18
    contact=left|right;free=~contact;phi=state[:,3];R=G+S
    Pext=sum(x['power_W'] for x in ext.values());Pcontact=ext['x-']['power_W']+ext['x+']['power_W']
    reaction=-float(phi[contact]@R[contact]);supg=float(phi@S);free_work=float(phi[free]@R[free])
    report={'order':order,'JE_W':W,'all_exterior_power_W':Pext,'contact_only_power_W':Pcontact,
            'gap_W':W-Pext,'gap_relative':(W-Pext)/W,'phi_divJ_volume_W':D,'internal_phi_flux_jump_W':jump,
            'integration_by_parts_remainder_W':W-Pext-D+jump,'volume_divJ_A':divsum,
            'reaction_terminal_power_W':reaction,'phi_weighted_SUPG_W':supg,'phi_weighted_free_residual_W':free_work,
            'algebraic_identity_remainder_W':W-reaction-supg+free_work,
            'reaction_minus_exterior_W':reaction-Pext,'free_residual_max_A':float(np.max(abs(R[free]))),
            'free_Galerkin_max_A':float(np.max(abs(G[free]))),'free_SUPG_max_A':float(np.max(abs(S[free]))),
            'reaction_left_A':float(R[left].sum()),'reaction_right_A':float(R[right].sum()),
            'exterior_area_normal_m2':area_normal.tolist(),'boundary_sets':ext,
            'element_ibp_max_abs_W':max(abs(x['element_ibp_error_W']) for x in elements)}
    with (out/f'elements_order{order}.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=list(elements[0]));writer.writeheader();writer.writerows(elements)
    np.savez_compressed(out/f'nodal_residual_order{order}.npz',xyz=xyz,G=G,S=S,R=R,contact=contact,phi=phi)
    fe.write_json(out/f'order{order}.json',report)
    return report


def basic_checks():
    # Affine HEX20 box, uniform E: constant J, no div J, exact boundary work.
    coords=fe.NATURAL*np.array([1e-6,2e-6,.5e-6]);xe=np.zeros((20,7));xe[:,3]=coords[:,0]*-1e5
    rr=rule(5);a=sample(coords,xe,rr,True)
    assert np.max(abs(a['divJ']))<1e-8*max(np.max(abs(a['J']))/1e-6,1)
    W=float((rr[1]*a['det'])@a['W']);P=0.
    for axis in range(3):
        for side in [-1,1]:
            r=rule(5,axis,side);b=sample(coords,xe,r)
            av=side*np.cross(b['jac'][:,(axis+1)%3,:],b['jac'][:,(axis+2)%3,:])*r[1][:,None]
            P-=float(b['phi']@np.sum(b['J']*av,axis=1))
    assert abs(W-P)/abs(W)<1e-12
    return {'affine_constant_current_power_relative_error':abs(W-P)/abs(W)}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--voltage',type=int,default=8)
    parser.add_argument('--orders',type=int,nargs='+',default=[3,5,7])
    parser.add_argument('--output',type=Path,default=fe.SOLVER/'runs/JE_POWER_AUDIT_20261015/V080')
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    run=ROOT/f'V{args.voltage*10:03d}'
    met=fe.read_metrics(run/'metrics.txt')
    assert met['heat_source']=='E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE'
    assert met['voltage_V']==args.voltage
    xyz,conn=fe.read_mesh(MESH);nodes=np.loadtxt(run/'nodes.csv',delimiter=',',skiprows=1);state=nodes[:,4:11]
    assert np.max(abs(nodes[:,1:4]-xyz))<1e-19
    faces,lookup=topology(xyz,conn)
    sources=[MESH,run/'nodes.csv',run/'metrics.txt',ROOT/f'input_V{args.voltage*10:03d}.nml',Path(fe.__file__),Path(__file__)]
    manifest=json.loads((fe.SOLVER/'build/ARC_JE_FINAL_20261012/build_manifest.json').read_text(encoding='utf-8-sig'))
    expected={x['path']:x['sha256'] for x in manifest['files']}
    matches={}
    for relative in ['src/physics/material_point_m.f90','src/physics/transport_m.f90','src/physics/zno_candidate_m.f90',
                     'src/fem/transport_supg_m.f90','src/fem/element_kernel_m.f90','src/post/independent_flux_audit_m.f90']:
        path=fe.SOLVER/relative;actual=fe.sha256(path);matches[relative]=actual==expected[relative];sources.append(path)
    assert all(matches.values()),matches
    fe.write_json(out/'provenance.json',{'sources':[{'path':str(p),'sha256':fe.sha256(p)} for p in sources],
                                       'source_matches_build':matches,'native_metrics':met,
                                       'nodes':len(xyz),'elements':len(conn),'exterior_faces':sum(len(x)==1 for x in faces.values()),
                                       'internal_faces':sum(len(x)==2 for x in faces.values()),'verification':basic_checks()})
    results=[]
    for order in args.orders:results.append(audit_order(xyz,conn,state,lookup,order,out))
    if 5 in args.orders:
        r=next(r for r in results if r['order']==5)
        for k,mk in [('JE_W','joule_W'),('all_exterior_power_W','terminal_power_W')]:
            assert abs(r[k]-met[mk])/abs(met[mk])<1e-8,(k,r[k],met[mk])
    fe.write_json(out/'summary.json',results)
    shutil.copy2(__file__,out/Path(__file__).name)
    print(json.dumps(results,indent=2),flush=True)
