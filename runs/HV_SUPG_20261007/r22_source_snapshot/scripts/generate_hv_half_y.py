"""Restrict an aligned full HEX20 mesh to y>=0, preserving element polynomials.

An optional full solution supplies a gauge-corrected INITIAL GUESS only. The
half problem must be solved and compared independently before it is used.
"""
import argparse, csv, json
from pathlib import Path
import numpy as np
from audit_hex20_boundary_geometry import read_mesh, derivatives
from run_hv_sweep import sha

TOL=1e-15
U=['ux_m','uy_m','uz_m']


def verify_reflection(xyz,conn,materials):
    keys=np.rint(xyz/1e-16).astype(np.int64)
    lookup={tuple(c):i for i,c in enumerate(keys)}
    if len(lookup)!=len(xyz):raise ValueError('Duplicate coordinates at reflection tolerance')
    reflected=keys.copy();reflected[:,1]*=-1
    try:mirror=np.array([lookup[tuple(c)] for c in reflected])
    except KeyError as exc:raise ValueError('Full mesh nodes lack y reflection symmetry') from exc
    if np.max(abs(xyz* [1,-1,1]-xyz[mirror]))>1e-16:
        raise ValueError('Full mesh reflection coordinate mismatch')
    elements={tuple(sorted(c)):i for i,c in enumerate(conn)}
    if len(elements)!=len(conn):raise ValueError('Duplicate elements')
    for i,c in enumerate(conn):
        partner=elements.get(tuple(sorted(mirror[c])))
        if partner is None or materials[i]!=materials[partner]:
            raise ValueError('Full element topology/materials lack reflection symmetry')
    return dict(matched_nodes=len(xyz),matched_elements=len(conn),
                maximum_coordinate_error_m=float(np.max(abs(xyz* [1,-1,1]-xyz[mirror]))))


def restrict_mesh(xyz,conn):
    ys=xyz[conn,1]
    if np.any((ys.min(axis=1)<-TOL)&(ys.max(axis=1)>TOL)):
        raise ValueError('An element crosses the symmetry plane')
    keep=np.flatnonzero(ys.min(axis=1)>=-TOL)
    if len(keep)*2!=len(conn):
        raise ValueError('Expected equal positive and negative half element counts')
    nodes=np.unique(conn[keep]); reverse=np.full(len(xyz),-1,dtype=int)
    reverse[nodes]=np.arange(len(nodes))
    # Do not snap coordinates: the existing polynomial geometry is retained.
    return nodes,keep,xyz[nodes],reverse[conn[keep]]


def rigid_gauge(xyz,u):
    """Remove the rigid motion fixed by half-plane uy and the A/B anchors."""
    plane=abs(xyz[:,1])<TOL
    # uy=ty+omega_z*x-omega_x*z. Scale columns for a conditioned fit.
    design=np.column_stack([np.ones(plane.sum()),xyz[plane,0]/1e-5,-xyz[plane,2]/1e-6])
    fit,_,rank,_=np.linalg.lstsq(design,u[plane,1],rcond=None)
    if rank!=3:raise ValueError('Insufficient symmetry-plane coordinates')
    def anchor(x):
        ids=np.flatnonzero(np.max(abs(xyz-[x,0,-.5e-6]),axis=1)<TOL)
        if len(ids)!=1:raise ValueError('Missing or duplicate half-domain anchor')
        return int(ids[0])
    a,b=anchor(-1e-5),anchor(1e-5)
    omega=np.array([fit[2]/1e-6,(u[a,2]-u[b,2])/(xyz[b,0]-xyz[a,0]),fit[1]/1e-5])
    trans=np.array([u[a,0]-omega[1]*xyz[a,2],fit[0],u[a,2]+omega[1]*xyz[a,0]])
    corrected=u-trans-np.cross(omega,xyz)
    return corrected,dict(translation_m=trans.tolist(),rotation_rad=omega.tolist(),
        maximum_symmetry_uy_m=float(abs(corrected[plane,1]).max()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seed-nodes',type=Path)
    p.add_argument('--perturb-hole',type=float,default=0.,help='Initial guess perturbation, vanishes at contacts')
    a=p.parse_args()
    if not np.isfinite(a.perturb_hole):raise ValueError('Nonfinite perturbation')
    mesh=a.mesh.absolute();xyz,conn,materials=read_mesh(mesh)
    reflection=verify_reflection(xyz,conn,materials)
    ids,eids,hxyz,hconn=restrict_mesh(xyz,conn)
    gp,_=np.polynomial.legendre.leggauss(3)
    pts=np.array([(x,y,z) for x in gp for y in gp for z in gp])
    det=np.linalg.det(np.einsum('qka,eaj->eqkj',derivatives(pts),hxyz[hconn]))
    if not det.min()>0:raise ValueError('Nonpositive sampled Jacobian')
    a.output.mkdir(parents=True,exist_ok=False);outmesh=a.output/'H1_Y_HALF.mesh'
    with outmesh.open('x',encoding='ascii') as f:
        f.write(f'{len(hxyz)} {len(hconn)}\n')
        for i,c in enumerate(hxyz,1):f.write(f'{i} '+ ' '.join(f'{v:.17e}' for v in c)+'\n')
        for i,(c,m) in enumerate(zip(hconn,materials[eids]),1):
            f.write(f'{i} {m} '+' '.join(str(v+1) for v in c)+'\n')
    mapping=a.output/'full_node_ids.npy';np.save(mapping,ids+1)
    manifest=dict(domain_symmetry='Y_HALF',status='MESH_AND_INITIAL_GUESS_ONLY',
        source_mesh=str(mesh),source_mesh_sha256=sha(mesh),mesh=str(outmesh.absolute()),mesh_sha256=sha(outmesh),
        node_map_sha256=sha(mapping),nodes=len(ids),elements=len(eids),minimum_sampled_jacobian=float(det.min()),
        reflection_audit=reflection)
    if a.seed_nodes:
        with a.seed_nodes.open(newline='',encoding='utf-8-sig') as f:
            reader=csv.DictReader(f);keys=reader.fieldnames;rows=list(reader)
        if len(rows)!=len(xyz) or [int(r['node']) for r in rows]!=list(range(1,len(xyz)+1)):
            raise ValueError('Seed node IDs do not match the full mesh')
        array=np.array([[float(r[k]) for k in keys[1:]] for r in rows])
        if not np.isfinite(array).all() or np.max(abs(array[:,:3]-xyz))>1e-16:
            raise ValueError('Seed geometry mismatch or nonfinite values')
        u=np.array([[float(r[k]) for k in U] for r in rows]);u,gauge=rigid_gauge(xyz,u)
        # A nonsymmetric solution is not a valid seed for an equivalence claim.
        if gauge['maximum_symmetry_uy_m']>1e-16:raise ValueError('Seed violates symmetry beyond rigid motion')
        with (a.output/'nodes.csv').open('x',newline='',encoding='ascii') as f:
            writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader()
            for i,old in enumerate(ids,1):
                r=dict(rows[old]);r['node']=i
                for k,v in zip(U,u[old]):r[k]=f'{v:.17e}'
                delta=a.perturb_hole*max(0.,1-(xyz[old,0]/1e-5)**2)
                r['eta_p']=f'{float(r["eta_p"])+delta:.17e}'
                r['p_m3']=f'{float(r["p_m3"])*np.exp(delta):.17e}'
                writer.writerow(r)
        manifest.update(seed_nodes=str(a.seed_nodes.absolute()),seed_sha256=sha(a.seed_nodes),
            seed_status='INITIAL_GUESS_REQUIRES_REEQUILIBRATION',rigid_gauge=gauge,
            hole_perturbation=a.perturb_hole,output_seed_sha256=sha(a.output/'nodes.csv'))
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(manifest,indent=2))


if __name__=='__main__':main()
