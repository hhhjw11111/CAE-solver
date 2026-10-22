"""Generate only the frozen F1 geometry/model matrix; never overwrite inputs.

Layer corner distances are frozen before solves. HEX20 radial/Cartesian midsides
are physical midpoints. Circular angular midsides retain the circle geometry.
"""
from pathlib import Path
import json, math

ROOT = Path(__file__).resolve().parents[1] / 'cases' / 'F1_20260906'
PAT = [(0,0,0),(2,0,0),(2,2,0),(0,2,0),(0,0,2),(2,0,2),(2,2,2),(0,2,2),
       (1,0,0),(2,1,0),(1,2,0),(0,1,0),(1,0,2),(2,1,2),(1,2,2),(0,1,2),
       (0,0,1),(2,0,1),(2,2,1),(0,2,1)]
R, W, TH = 3e-6, 10e-6, 1e-6

def skin(first, depth):
    # Bound growth by 2.0; distribute the final interval without a tiny sliver.
    vals=[0.0]; step=first
    while vals[-1]+2.0*step < depth:
        vals.append(vals[-1]+step); step*=2.0
    vals.append(depth)
    return vals

def mids(corners):
    out=[]
    for a,b in zip(corners,corners[1:]): out += [a,(a+b)/2]
    return out+[corners[-1]]

def write_mesh(path,coords,elements):
    with path.open('x',encoding='ascii') as f:
        f.write(f'{len(coords)} {len(elements)}\n')
        for i,xyz in enumerate(coords,1): f.write(f'{i} '+' '.join(f'{v:.17e}' for v in xyz)+'\n')
        for i,conn in enumerate(elements,1): f.write(f'{i} 1 '+' '.join(map(str,conn))+'\n')

def generate(level,hole):
    first=14e-9/2**level; nt=[16,24,32][level]
    layer=skin(first,0.5e-6)
    zskin=skin(first,0.5e-6)
    zcorn=[-TH/2+d for d in zskin]+[TH/2-d for d in reversed(zskin[:-1])]
    zs=mids(zcorn); nz=len(zcorn)-1
    ids={}; coords=[]; elements=[]
    if hole:
        # Ring outer radius is an exact radial corner surface. Four bulk layers
        # connect its outer edge to the outer-contact skin.
        nskin=len(layer)-1; nr=2*nskin+4
        def corners(angle):
            outer=W/max(abs(math.cos(angle)),abs(math.sin(angle)))
            a=R+layer[-1]; b=outer-layer[-1]
            return [R+d for d in layer]+[a+(b-a)*i/4 for i in range(1,4)]+[outer-d for d in reversed(layer)]
        for ez in range(nz):
            for et in range(nt):
                for er in range(nr):
                    conn=[]
                    for pr,pt,pz in PAT:
                        ir,it,iz=2*er+pr,(2*et+pt)%(2*nt),2*ez+pz
                        key=(ir,it,iz)
                        if key not in ids:
                            angle=math.pi*it/nt; rs=mids(corners(angle)); rad=rs[ir]
                            ids[key]=len(coords)+1;coords.append((rad*math.cos(angle),rad*math.sin(angle),zs[iz]))
                        conn.append(ids[key])
                    elements.append(conn)
        details={'circumferential':nt,'radial':nr,'thickness':nz}
    else:
        nb=[4,6,8][level]
        corners=[-W+d for d in layer]+[-W+layer[-1]+(2*W-2*layer[-1])*i/nb for i in range(1,nb)]+[W-d for d in reversed(layer)]
        axis=mids(corners); nx=len(corners)-1
        for ez in range(nz):
            for ey in range(nx):
                for ex in range(nx):
                    conn=[]
                    for px,py,pz in PAT:
                        key=(2*ex+px,2*ey+py,2*ez+pz)
                        if key not in ids:
                            ids[key]=len(coords)+1;coords.append((axis[key[0]],axis[key[1]],zs[key[2]]))
                        conn.append(ids[key])
                    elements.append(conn)
        details={'x':nx,'y':nx,'thickness':nz}
    name=f'H{int(hole)}_L{level}'; write_mesh(ROOT/(name+'.mesh'),coords,elements)
    return {'name':name,'first_layer_m':first,'skin_corners_m':layer,'z_corners_m':zcorn,
            'counts':details,'nodes':len(coords),'elements':len(elements)}

def case(hole,level,model,voltage,tag):
    name=f'{tag}_H{int(hole)}_{model}_L{level}'
    text=f'''&case_control
 schema_version=2, case_id='{name}', model_level='M2_THERMO_TWO_CARRIER',
 mesh_file='H{int(hole)}_L{level}.mesh'
/
&problem_control
 problem_type='{'CENTER_HOLE_PLATE' if hole else 'PLAIN_PLATE_CONTROL'}',
 material_profile='ZNO_F1_{model}', target_voltage={-voltage:.17e},
 target_voltage_factor=1.0, target_joule_factor=1.0, hole_radius={R if hole else 0:.17e}
/
&solver_control
 field_scales=1e-10,1e-10,1e-10,3.0,0.1,0.1,10.0,
 max_newton=40, residual_tolerance=1e-9, update_tolerance=1e-9,
 continuation_initial_step=0.25, continuation_minimum_step=0.015625,
 continuation_maximum_step=0.5, continuation_maximum_attempts=100,
 continuation_maximum_newton_per_attempt=30
/
'''
    with (ROOT/(name+'.nml')).open('x',encoding='ascii') as f: f.write(text)
    return name+'.nml'

def main():
    ROOT.mkdir(parents=True,exist_ok=False)
    meshes=[generate(level,True) for level in range(3)]
    meshes.append(generate(2,False))
    cases=[case(True,0,m,0.1,'SMALL') for m in ['LQ','NT']]
    cases += [case(True,i,'NT',3.0,'MESH') for i in range(3)]
    cases += [case(h,2,m,3.0,'MAIN') for h in [False,True] for m in ['LQ','N0','NT'] if not(h and m=='NT')]
    manifest={'plan':'F1_20260906','material_evidence':'NUMERICAL_CANDIDATE_NOT_EXPERIMENTALLY_CALIBRATED',
              'mesh_design_not_mesh_acceptance':True,'meshes':meshes,'cases':cases,
              'reuse':'MESH_H1_NT_L2 is MAIN H1-NT; repeat only after all acceptance checks.'}
    with (ROOT/'inputs_manifest.json').open('x',encoding='utf-8') as f: json.dump(manifest,f,indent=2)
    print(json.dumps(manifest,indent=2))

if __name__=='__main__': main()
