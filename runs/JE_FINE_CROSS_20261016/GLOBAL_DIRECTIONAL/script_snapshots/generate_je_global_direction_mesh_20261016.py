"""Refine the entire radial or thickness logical coordinate, retaining the annulus.

The level-zero mesh reproduces A32/R14/Z14 exactly. Radial refinement increases
both the five annular intervals and the nine intervals from annulus to electrode.
The r=3.5 um annular measurement boundary remains an element boundary. Thickness
refinement resamples all 14 original intervals. No solver or material is changed.
"""
from pathlib import Path
import argparse,json,math
import numpy as np
from generate_f1_main_cases import PAT,R,W,TH,skin,mids,write_mesh
from generate_hv_surface_mesh import audit_geometry
from audit_hex20_boundary_geometry import derivatives
from run_hv_sweep import sha

def resample(corners,count):
    if count==len(corners)-1:return list(corners)
    return np.interp(np.linspace(0,len(corners)-1,count+1),np.arange(len(corners)),corners).tolist()

def generate(direction,level):
    assert direction in ['baseline','radial','thickness'] and level in [0,1,2]
    assert direction!='baseline' or level==0
    rl=level if direction=='radial' else 0;zl=level if direction=='thickness' else 0
    layer=skin(14e-9,.5e-6);zskin=skin(3.5e-9,.5e-6)
    oldz=[-TH/2+d for d in zskin]+[TH/2-d for d in reversed(zskin[:-1])]
    zcorn=resample(oldz,14+2*zl);zs=mids(zcorn)
    nt=32;nr=14+2*rl;nz=14+2*zl
    def corners(angle):
        outer=W/max(abs(math.cos(angle)),abs(math.sin(angle)))
        a=R+layer[-1];b=outer-layer[-1]
        inner=[R+d for d in layer]
        rest=[a]+[a+(b-a)*i/4 for i in range(1,4)]+[outer-d for d in reversed(layer)]
        return resample(inner,5+rl)+resample(rest,9+rl)[1:]
    radial=[mids(corners(math.pi*it/nt)) for it in range(2*nt)]
    ids={};coords=[];elements=[]
    for ez in range(nz):
        for et in range(nt):
            for er in range(nr):
                conn=[]
                for pr,pt,pz in PAT:
                    ir,it,iz=2*er+pr,(2*et+pt)%(2*nt),2*ez+pz;key=(ir,it,iz)
                    if key not in ids:
                        angle=math.pi*it/nt;rad=radial[it][ir]
                        ids[key]=len(coords)+1;coords.append((rad*math.cos(angle),rad*math.sin(angle),zs[iz]))
                    conn.append(ids[key])
                elements.append(conn)
    profiles={'radial_axis_corners_m':corners(0.),'radial_diagonal_corners_m':corners(math.pi/4),'z_corners_m':zcorn}
    manifest=dict(direction=direction,level=level,counts={'circumferential':nt,'radial':nr,'thickness':nz},
                  nodes=len(coords),elements=len(elements),annulus_outer_radius_m=R+.5e-6,
                  annulus_corner_index=5+rl,annulus_intervals=5+rl,remaining_radial_intervals=9+rl,
                  profiles=profiles,status='MESH_DESIGN_NOT_NUMERICAL_ACCEPTANCE',
                  scope='Full-direction logical-coordinate refinement; other directions fixed; nonnested meshes')
    manifest['maximum_radial_axis_step_m']=float(np.diff(profiles['radial_axis_corners_m']).max())
    manifest['maximum_radial_diagonal_step_m']=float(np.diff(profiles['radial_diagonal_corners_m']).max())
    manifest['maximum_z_step_m']=float(np.diff(zcorn).max())
    manifest['geometry_audit']=audit_geometry(coords,elements)
    q,w=np.polynomial.legendre.leggauss(3);pts=np.array([(a,b,c) for a in q for b in q for c in q])
    weights=np.array([a*b*c for a in w for b in w for c in w])
    jac=np.einsum('qka,eaj->eqkj',derivatives(pts),np.asarray(coords)[np.asarray(elements)-1])
    manifest['volume_3x3x3_m3']=float(np.sum(np.linalg.det(jac)*weights))
    return coords,elements,manifest

def write(directory,direction,level):
    coords,elements,m=generate(direction,level)
    directory.mkdir(parents=True,exist_ok=False);mesh=directory/'H1.mesh';write_mesh(mesh,coords,elements)
    m['mesh_sha256']=sha(mesh);m['generator_sha256']=sha(Path(__file__))
    with (directory/'manifest.json').open('x',encoding='utf-8') as f:json.dump(m,f,indent=2)
    return mesh,m

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--direction',choices=['baseline','radial','thickness'],required=True)
    p.add_argument('--level',type=int,choices=[0,1,2],required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();mesh,m=write(a.output,a.direction,a.level);print(json.dumps(m,indent=2))
