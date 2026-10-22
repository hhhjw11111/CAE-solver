"""Independent r=2/3/4 um meshes for the approved JE paper scope.

Only geometry is generated and checked. There is no native solve, input NML,
solution acceptance, or mutation of the older generators. At r=3 um the node
and element ordering and coordinate arithmetic reproduce the frozen generator.
"""
from pathlib import Path
import argparse
import json
import math
import numpy as np
from generate_f1_main_cases import PAT, W, TH, skin, mids, write_mesh
from generate_je_global_direction_mesh_20261016 import resample
from audit_hex20_boundary_geometry import derivatives
from transfer_hv_initial_guess import shapes
from run_hv_sweep import sha

RING_WIDTH = .5e-6
COORDINATE_TOLERANCE_M = 1e-18
GEOMETRY_RELATIVE_TOLERANCE = 1e-10
ROOT = Path(__file__).resolve().parents[1]


def save(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def dependencies():
    names = ['generate_je_paper_radius_mesh_20261018.py',
             'generate_je_reconstruction_mesh_20261018.py',
             'generate_je_global_direction_mesh_20261016.py', 'generate_f1_main_cases.py',
             'audit_hex20_boundary_geometry.py', 'transfer_hv_initial_guess.py', 'run_hv_sweep.py']
    return {str(ROOT/'scripts'/name): sha(ROOT/'scripts'/name) for name in names}


def _volume(xyz, conn, order, annulus_elements):
    gp, gw = np.polynomial.legendre.leggauss(order)
    points = np.array([(r,s,t) for r in gp for s in gp for t in gp])
    weights = np.array([a*b*c for a in gw for b in gw for c in gw])
    dn = derivatives(points)
    volume = ring_volume = 0.
    minimum = math.inf
    for begin in range(0, len(conn), 128):
        part = conn[begin:begin+128]
        jac = np.einsum('qka,eaj->eqkj', dn, xyz[part])
        det = np.linalg.det(jac)
        minimum = min(minimum, float(det.min()))
        element_volume = det @ weights
        volume += float(element_volume.sum())
        ring_volume += float(element_volume[annulus_elements[begin:begin+128]].sum())
    return dict(volume_m3=volume, annulus_volume_m3=ring_volume, minimum_jacobian=minimum)


def audit_geometry(coords, elements, radius, angular, radial, thickness, annulus_index):
    """Audit the actual HEX20 geometry, including the quadratic-circle correction."""
    xyz, conn = np.asarray(coords), np.asarray(elements)-1
    if not np.isfinite(xyz).all() or len({tuple(x) for x in xyz}) != len(xyz):
        raise ValueError('Nonfinite or duplicate node coordinates')
    if any(len(set(row)) != 20 for row in conn):
        raise ValueError('Repeated element nodes')
    if not np.array_equal(np.unique(conn), np.arange(len(xyz))):
        raise ValueError('Unused or invalid mesh nodes')
    expected_bounds = np.array([[-W,-W,-TH/2],[W,W,TH/2]])
    bounds = np.vstack([xyz.min(axis=0),xyz.max(axis=0)])
    bound_error = float(np.max(abs(bounds-expected_bounds)))
    targets = np.array([[-W,-W,-TH/2],[-W,W,-TH/2],[-W,-W,TH/2]])
    anchors = [int((np.max(abs(xyz-t),axis=1)<COORDINATE_TOLERANCE_M).sum()) for t in targets]
    contacts = [int((abs(xyz[:,0]-x)<COORDINATE_TOLERANCE_M).sum()) for x in [-W,W]]
    pattern = np.asarray(PAT)
    face_masks = {(axis,side): np.flatnonzero(pattern[:,axis] == (0 if side == -1 else 2))
                  for axis in range(3) for side in (-1,1)}
    faces = {}
    for element, nodes in enumerate(conn):
        for (axis,side), mask in face_masks.items():
            key = tuple(sorted(nodes[mask]))
            faces.setdefault(key,[]).append((element,axis,side))
    if not all(len(owners) in (1,2) for owners in faces.values()):
        raise ValueError('Nonmanifold face ownership')
    exterior = {name: [] for name in ['x-','x+','y-','y+','z-','z+','hole']}
    annulus_faces = 0
    annulus_node_error = hole_node_error = 0.
    ring_radius = radius+RING_WIDTH
    for ids, owners in faces.items():
        nodes = xyz[list(ids)]
        distances = np.linalg.norm(nodes[:,:2],axis=1)
        if np.max(abs(distances-ring_radius)) < COORDINATE_TOLERANCE_M:
            if len(owners) != 2:
                raise ValueError('Annulus interface must have exactly two owners')
            annulus_faces += 1
            annulus_node_error = max(annulus_node_error,float(np.max(abs(distances-ring_radius))))
        if len(owners) == 2:
            continue
        kind = 'hole'
        for axis, label in enumerate(('x','y','z')):
            for side, limit in [('-',expected_bounds[0,axis]),('+',expected_bounds[1,axis])]:
                if np.max(abs(nodes[:,axis]-limit)) < COORDINATE_TOLERANCE_M:
                    kind = label+side
        if kind == 'hole':
            hole_node_error = max(hole_node_error,float(np.max(abs(distances-radius))))
        exterior[kind].append(owners[0])
    expected_faces = {'x-':angular//4*thickness,'x+':angular//4*thickness,
                      'y-':angular//4*thickness,'y+':angular//4*thickness,
                      'z-':angular*radial,'z+':angular*radial,'hole':angular*thickness}
    face_counts = {k:len(v) for k,v in exterior.items()}
    if face_counts != expected_faces or annulus_faces != angular*thickness:
        raise ValueError('Incorrect physical boundary or annulus-interface topology')
    if bound_error > COORDINATE_TOLERANCE_M or hole_node_error > COORDINATE_TOLERANCE_M:
        raise ValueError('Geometry coordinates disagree with specified radius or outer boundary')
    if anchors != [1,1,1] or min(contacts) < 8:
        raise ValueError('Missing rigid-body anchors or electrode nodes')
    ring_elements = np.arange(len(conn)) % radial < annulus_index
    volumes = {order:_volume(xyz,conn,order,ring_elements) for order in (3,5)}
    # A quadratic circular arc with endpoint angles +/-h and midpoint angle 0
    # is x/R=1-(1-cos(h))*t^2, y/R=sin(h)*t. Its sector area is
    # R^2*sin(h)*(4-cos(h))/3. This checks the actual FE boundary independently.
    h = math.pi/angular
    circle_area_factor = angular*math.sin(h)*(4-math.cos(h))/3
    expected_volume = (4*W*W-circle_area_factor*radius*radius)*TH
    expected_ring = circle_area_factor*(ring_radius*ring_radius-radius*radius)*TH
    ideal_volume = (4*W*W-math.pi*radius*radius)*TH
    ideal_ring = math.pi*(ring_radius*ring_radius-radius*radius)*TH
    for item in volumes.values():
        item['quadratic_boundary_volume_relative_error'] = abs(item['volume_m3']-expected_volume)/expected_volume
        item['quadratic_boundary_annulus_relative_error'] = abs(item['annulus_volume_m3']-expected_ring)/expected_ring
        if item['minimum_jacobian'] <= 0 or max(item['quadratic_boundary_volume_relative_error'],
                                               item['quadratic_boundary_annulus_relative_error']) > GEOMETRY_RELATIVE_TOLERANCE:
            raise ValueError('Nonpositive sampled Jacobian or incorrect FE volume/ring volume')
    # Independent surface-divergence volume and closure of outward area vectors.
    gp, gw = np.polynomial.legendre.leggauss(3)
    area_vector = np.zeros(3)
    area_sum = divergence_volume = 0.
    all_exterior = [face for group in exterior.values() for face in group]
    for axis in range(3):
        others = [(axis+1)%3,(axis+2)%3]
        for side in (-1,1):
            selected = [e for e,a,s in all_exterior if a==axis and s==side]
            points = np.zeros((9,3)); points[:,axis] = side
            points[:,others] = np.array([(a,b) for a in gp for b in gp])
            weights = np.array([a*b for a in gw for b in gw])
            dn, shape = derivatives(points), shapes(points)
            for begin in range(0,len(selected),128):
                x = xyz[conn[selected[begin:begin+128]]]
                jac = np.einsum('qka,eaj->eqkj',dn,x)
                av = side*np.cross(jac[:,:,others[0],:],jac[:,:,others[1],:])*weights[None,:,None]
                sample = np.einsum('qa,eaj->eqj',shape,x)
                area_vector += av.sum(axis=(0,1))
                area_sum += float(np.linalg.norm(av,axis=2).sum())
                divergence_volume += float(np.sum(sample*av))/3
    closure = float(np.linalg.norm(area_vector)/area_sum)
    surface_error = abs(divergence_volume-expected_volume)/expected_volume
    if max(closure,surface_error) > GEOMETRY_RELATIVE_TOLERANCE:
        raise ValueError('Exterior surface orientation, closure or divergence-volume mismatch')
    return dict(status='GEOMETRY_CHECKS_PASSED_NOT_SOLUTION_ACCEPTANCE',
        geometry_coordinate_tolerance_m=COORDINATE_TOLERANCE_M,
        numerical_geometry_relative_tolerance=GEOMETRY_RELATIVE_TOLERANCE,
        bounds_m=bounds.tolist(),bound_coordinate_error_m=bound_error,
        anchor_counts=anchors,contact_counts=contacts,exterior_face_counts=face_counts,
        annulus_interface_faces=annulus_faces,annulus_interface_node_radius_error_m=annulus_node_error,
        hole_node_radius_error_m=hole_node_error,quadrature={str(k):v for k,v in volumes.items()},
        expected_quadratic_boundary_volume_m3=expected_volume,expected_quadratic_annulus_volume_m3=expected_ring,
        ideal_circle_volume_m3=ideal_volume,ideal_circle_annulus_volume_m3=ideal_ring,
        ideal_circle_volume_relative_geometry_error=(expected_volume-ideal_volume)/ideal_volume,
        ideal_circle_annulus_relative_geometry_error=(expected_ring-ideal_ring)/ideal_ring,
        quadratic_arc_relative_max_inward_sag=1-math.sqrt(1-(1-math.cos(h))**2/4),
        exterior_area_vector_m2=area_vector.tolist(),exterior_area_vector_relative_closure=closure,
        divergence_volume_m3=divergence_volume,divergence_volume_relative_error=surface_error,
        limitation='Positive Jacobians sampled at tensor orders 3 and 5; geometry checks are not FEM solution or spatial-convergence acceptance.')


def generate(radius_um=3, angular=40, radial=14, thickness=14):
    if radius_um not in (2,3,4):
        raise ValueError('Only approved radii 2, 3 and 4 um are supported')
    if angular not in (32,40,48) or radial not in (14,16,18,20,22) or thickness not in (14,16,18):
        raise ValueError('Only declared FULL A32/40/48, R14/16/18/20/22, Z14/16/18 are supported')
    # Decimal conversion preserves exactly the legacy literal 3e-6 at radius 3.
    radius = float(f'{radius_um:g}e-6')
    level = (radial-14)//2
    layer = skin(14e-9,RING_WIDTH); zskin = skin(3.5e-9,.5e-6)
    oldz = [-TH/2+d for d in zskin]+[TH/2-d for d in reversed(zskin[:-1])]
    zcorn = resample(oldz,thickness); zs = mids(zcorn)
    def corners(angle):
        outer = W/max(abs(math.cos(angle)),abs(math.sin(angle)))
        a,b = radius+layer[-1],outer-layer[-1]
        inner = [radius+d for d in layer]
        rest = [a]+[a+(b-a)*i/4 for i in range(1,4)]+[outer-d for d in reversed(layer)]
        return resample(inner,5+level)+resample(rest,9+level)[1:]
    profiles = [mids(corners(math.pi*it/angular)) for it in range(2*angular)]
    ids,coords,elements = {},[],[]
    for ez in range(thickness):
        for et in range(angular):
            for er in range(radial):
                conn=[]
                for pr,pt,pz in PAT:
                    ir,it,iz = 2*er+pr,(2*et+pt)%(2*angular),2*ez+pz
                    key=ir,it,iz
                    if key not in ids:
                        angle,rad = math.pi*it/angular,profiles[it][ir]
                        ids[key]=len(coords)+1
                        coords.append((rad*math.cos(angle),rad*math.sin(angle),zs[iz]))
                    conn.append(ids[key])
                elements.append(conn)
    audit = audit_geometry(coords,elements,radius,angular,radial,thickness,5+level)
    metadata=dict(domain_symmetry='FULL',hole_radius_m=radius,radius_um=radius_um,
        ring_width_m=RING_WIDTH,annulus_inner_radius_m=radius,annulus_outer_radius_m=radius+RING_WIDTH,
        annulus_corner_index=5+level,counts=dict(circumferential=angular,radial=radial,thickness=thickness),
        nodes=len(coords),elements=len(elements),profiles=dict(radial_axis_corners_m=corners(0.),
        radial_diagonal_corners_m=corners(math.pi/4),z_corners_m=zcorn),
        maximum_radial_axis_step_m=float(np.diff(corners(0.)).max()),maximum_z_step_m=float(np.diff(zcorn).max()),
        volume_3x3x3_m3=audit['quadrature']['3']['volume_m3'],geometry_audit=audit,
        status='MESH_DESIGN_NOT_NUMERICAL_ACCEPTANCE',
        scope='Radius-parametric full-domain geometry; 0.5 um ring retained. Separate radius-specific field/SUPG evidence required.')
    return coords,elements,metadata


def write(directory, radius_um=3, angular=40, radial=14, thickness=14):
    directory=Path(directory).resolve()
    if directory.exists():
        raise FileExistsError('Refusing to overwrite mesh directory: '+str(directory))
    coords,elements,metadata=generate(radius_um,angular,radial,thickness)
    directory.mkdir(parents=True,exist_ok=False)
    mesh=directory/'H1.mesh';write_mesh(mesh,coords,elements)
    metadata.update(mesh_sha256=sha(mesh),generator_sha256=sha(Path(__file__)),dependencies=dependencies())
    save(directory/'manifest.json',metadata)
    return mesh,metadata


def validate_preparation(directory):
    """Small geometry-only matrix; three exact legacy reproductions plus r2/r4."""
    import generate_je_reconstruction_mesh_20261018 as legacy
    directory=Path(directory).resolve();directory.mkdir(parents=True,exist_ok=False)
    cases=[(3,32,14,14),(3,40,20,16),(3,48,22,18),(2,32,14,14),(4,32,14,14)]
    frozen=dependencies()
    save(directory/'protocol.json',dict(cases=cases,dependencies=frozen,
        coordinate_tolerance_m=COORDINATE_TOLERANCE_M,geometry_relative_tolerance=GEOMETRY_RELATIVE_TOLERANCE,
        checks=['r3 exact mesh byte comparison at three A/R/Z settings','orders 3/5 sampled positive Jacobians',
                'explicit quadratic boundary volume and ring volume','outer boundaries and anchors',
                'manifold face ownership and shared ring interface','outward surface closure and divergence volume'],
        scope='Geometry preparation only; no production inputs, native processes, solution states or production PASS'))
    records=[]
    for r,a,nr,nz in cases:
        tag=f'r{r}_A{a}_R{nr}_Z{nz}'
        mesh,metadata=write(directory/tag,r,a,nr,nz)
        result=dict(case=tag,mesh=str(mesh),mesh_sha256=sha(mesh),
                    geometry_status=metadata['geometry_audit']['status'],geometry_audit=metadata['geometry_audit'])
        if r==3:
            reference,_=legacy.write(directory/(tag+'_legacy_reference'),a,nr,nz)
            result.update(legacy_reference=str(reference),legacy_reference_sha256=sha(reference),
                          byte_identical=mesh.read_bytes()==reference.read_bytes())
            if not result['byte_identical']:
                save(directory/(tag+'_failure.json'),result)
                raise ValueError('r3 mesh does not byte-reproduce frozen generator')
        records.append(result)
        print(json.dumps({k:v for k,v in result.items() if k!='geometry_audit'}),flush=True)
    unchanged=all(sha(Path(p))==digest for p,digest in frozen.items())
    save(directory/'report.json',dict(status='GEOMETRY_PREPARATION_COMPLETE_NOT_PRODUCTION_ACCEPTANCE',
        geometry_checks_passed=unchanged and len(records)==5,dependencies_unchanged=unchanged,
        original_generators_modified=False,native_solver_started=False,production_inputs_generated=False,
        records=records,limitation='Geometry tests do not certify carrier, thermal, electric, stress, power or SUPG convergence.'))
    save(directory/'artifact_hashes.json',{str(p.relative_to(directory)):sha(p) for p in directory.rglob('*') if p.is_file()})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius-um',type=int,choices=(2,3,4))
    p.add_argument('--counts',nargs=3,type=int,metavar=('ANGULAR','RADIAL','THICKNESS'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--validate-preparation',action='store_true')
    args=p.parse_args()
    if args.validate_preparation:
        if args.radius_um is not None or args.counts is not None:
            p.error('Preparation matrix owns its cases; omit radius/count options')
        validate_preparation(args.output)
    else:
        if args.radius_um is None or args.counts is None:
            p.error('A single mesh requires --radius-um and --counts')
        mesh,metadata=write(args.output,args.radius_um,*args.counts)
        print(json.dumps(dict(mesh=str(mesh),**metadata),indent=2))


if __name__=='__main__':main()
