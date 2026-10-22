"""Independent in-plane and thickness refinement of the frozen F1 hole geometry.

Creates a new mesh, never edits frozen cases. The original (xy_level=z_level)
levels 0--2 reproduce the F1 generator exactly. A mesh is not a solved field or
an accuracy certificate; refinement directions must be assessed separately.
"""
from pathlib import Path
import argparse, json, math
import numpy as np
from generate_f1_main_cases import PAT, R, W, TH, skin, mids, write_mesh
from audit_hex20_boundary_geometry import derivatives
from run_hv_sweep import sha


def refine_surface_intervals(corners, subdivisions, depth):
    if subdivisions not in (1,2,3,4) or not math.isfinite(depth) or not 0<depth<=TH/2:
        raise ValueError('Invalid nested surface refinement')
    if subdivisions==1:return corners
    result=[]
    for a,b in zip(corners,corners[1:]):
        distance=min(abs(a+TH/2),abs(a-TH/2),abs(b+TH/2),abs(b-TH/2))
        count=subdivisions if distance<depth else 1
        result.extend(a+(b-a)*j/count for j in range(count))
    return result+[corners[-1]]


def generate(xy_level, z_level,angular_count=None,z_subdivisions=1,z_refine_depth=1e-7):
    if xy_level not in (0,1,2) or z_level not in range(5):
        raise ValueError('xy_level must be 0--2; z_level must be 0--4')
    first_xy=14e-9/2**xy_level; first_z=14e-9/2**z_level
    nt=[16,24,32][xy_level] if angular_count is None else angular_count
    if not isinstance(nt,int) or nt<16 or nt>256 or nt%8:
        raise ValueError('angular_count must be a multiple of 8 in [16,256]')
    layer=skin(first_xy,.5e-6)
    zskin=skin(first_z,.5e-6)
    zcorn=[-TH/2+d for d in zskin]+[TH/2-d for d in reversed(zskin[:-1])]
    zcorn=refine_surface_intervals(zcorn,z_subdivisions,z_refine_depth)
    zs=mids(zcorn); nz=len(zcorn)-1; nr=2*(len(layer)-1)+4
    ids={}; coords=[]; elements=[]
    def corners(angle):
        outer=W/max(abs(math.cos(angle)),abs(math.sin(angle)))
        a=R+layer[-1]; b=outer-layer[-1]
        return ([R+d for d in layer]+[a+(b-a)*i/4 for i in range(1,4)]+
                [outer-d for d in reversed(layer)])
    for ez in range(nz):
        for et in range(nt):
            for er in range(nr):
                conn=[]
                for pr,pt,pz in PAT:
                    ir,it,iz=2*er+pr,(2*et+pt)%(2*nt),2*ez+pz
                    key=(ir,it,iz)
                    if key not in ids:
                        angle=math.pi*it/nt; rad=mids(corners(angle))[ir]
                        ids[key]=len(coords)+1
                        coords.append((rad*math.cos(angle),rad*math.sin(angle),zs[iz]))
                    conn.append(ids[key])
                elements.append(conn)
    manifest=dict(xy_level=xy_level,z_level=z_level,first_xy_m=first_xy,first_z_m=first_z/z_subdivisions,base_first_z_m=first_z,
        z_subdivisions=z_subdivisions,z_refine_depth_m=z_refine_depth,
        counts=dict(circumferential=nt,radial=nr,thickness=nz),z_corners_m=zcorn,
        nodes=len(coords),elements=len(elements),status='MESH_DESIGN_NOT_NUMERICAL_ACCEPTANCE')
    return coords,elements,manifest


def audit_geometry(coords,elements):
    xyz=np.array(coords); conn=np.array(elements)-1
    gp,_=np.polynomial.legendre.leggauss(3)
    points=np.array([(x,y,z) for x in gp for y in gp for z in gp])
    jac=np.einsum('qka,eaj->eqkj',derivatives(points),xyz[conn])
    mindet=float(np.linalg.det(jac).min())
    targets=np.array([[-W,-W,-TH/2],[-W,W,-TH/2],[-W,-W,TH/2]])
    anchors=[int((np.max(abs(xyz-t),axis=1)<1e-15).sum()) for t in targets]
    contacts=[int((abs(xyz[:,0]-x)<1e-15).sum()) for x in [-W,W]]
    if not mindet>0 or anchors!=[1,1,1] or min(contacts)<8:
        raise ValueError(f'Invalid geometry: det={mindet}, anchors={anchors}, contacts={contacts}')
    return dict(minimum_jacobian_3x3x3=mindet,anchor_counts=anchors,contact_counts=contacts)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--xy-level',type=int,choices=range(3),required=True)
    p.add_argument('--z-level',type=int,choices=range(5),required=True)
    p.add_argument('--z-subdivisions',type=int,choices=[1,2,3,4],default=1)
    p.add_argument('--z-refine-depth-m',type=float,default=1e-7)
    p.add_argument('--angular-count',type=int,help='Optional independent circumferential refinement; radial and z layers unchanged')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();coords,elements,manifest=generate(a.xy_level,a.z_level,a.angular_count,a.z_subdivisions,a.z_refine_depth_m)
    manifest['geometry_audit']=audit_geometry(coords,elements)
    a.output.mkdir(parents=True,exist_ok=False)
    mesh=a.output/f'H1_XY{a.xy_level}_Z{a.z_level}.mesh'
    write_mesh(mesh,coords,elements);manifest['mesh_sha256']=sha(mesh)
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps(manifest,indent=2))


if __name__=='__main__':main()
