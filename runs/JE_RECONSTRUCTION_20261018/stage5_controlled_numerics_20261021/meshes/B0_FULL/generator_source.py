"""Independent geometry-preserving bulk/annulus refinement for the JE study.

The old radial profile jumps from the r=3.5 um audit interface to a >1 um
bulk interval.  This family separately allocates hole skin, intermediate bulk
and contact skin intervals and bounds the first bulk step.  No solved field or
acceptance threshold is used to construct the mesh.  Historical generators stay
unchanged.  The r=3.5 um interface remains a corner on every angular ray.
"""
from pathlib import Path
import argparse
import json
import math
import numpy as np
from generate_f1_main_cases import PAT, R, W, TH, skin, mids, write_mesh
from generate_je_global_direction_mesh_20261016 import resample
from audit_hex20_boundary_geometry import derivatives
from run_hv_sweep import sha

COUNTS = ((48, 26, 18), (64, 30, 20), (80, 34, 22))


def profiles(level, angle):
    if level not in (0, 1, 2):
        raise ValueError('Three predeclared levels only')
    outer = W/max(abs(math.cos(angle)), abs(math.sin(angle)))
    inner_n, bulk_n, outer_n = 9+level, 8+2*level, 9+level
    inner = resample([R+d for d in skin(14e-9, .5e-6)], inner_n)
    contact = resample([outer-d for d in reversed(skin(14e-9, .5e-6))], outer_n)
    # Cosine nodes give small adjacent intervals at both ends of the bulk.
    # The r=3.5 um annular audit boundary is an exact corner, not a fitted mask.
    t = np.linspace(0, 1, bulk_n+1)
    bulk = (R+.5e-6)+(outer-R-1e-6)*(1-np.cos(math.pi*t))/2
    corners = inner+bulk[1:-1].tolist()+contact
    if len(corners) != COUNTS[level][1]+1 or np.any(np.diff(corners) <= 0):
        raise AssertionError('Invalid predeclared radial profile')
    zskin = skin(3.5e-9, .5e-6)
    original_z = [-TH/2+d for d in zskin]+[TH/2-d for d in reversed(zskin[:-1])]
    zcorn = resample(original_z, COUNTS[level][2])
    return corners, zcorn, inner_n


def generate(level):
    nt, nr, nz = COUNTS[level]
    radial = [mids(profiles(level, math.pi*it/nt)[0]) for it in range(2*nt)]
    z = mids(profiles(level, 0.)[1])
    ids, coords, elements = {}, [], []
    for ez in range(nz):
        for et in range(nt):
            for er in range(nr):
                element = []
                for pr, pt, pz in PAT:
                    ir, it, iz = 2*er+pr, (2*et+pt)%(2*nt), 2*ez+pz
                    key = ir, it, iz
                    if key not in ids:
                        angle, radius = math.pi*it/nt, radial[it][ir]
                        ids[key] = len(coords)+1
                        coords.append((radius*math.cos(angle), radius*math.sin(angle), z[iz]))
                    element.append(ids[key])
                elements.append(element)
    xyz, conn = np.asarray(coords), np.asarray(elements)-1
    # Both assembly and denser geometry sampling are checked before any solve.
    min_det, volume = {}, None
    for order in (3, 5):
        q, w = np.polynomial.legendre.leggauss(order)
        points = np.array([(a,b,c) for a in q for b in q for c in q])
        weights = np.array([a*b*c for a in w for b in w for c in w])
        dn = derivatives(points)
        minimum, integral = float('inf'), 0.
        for start in range(0, len(conn), 128):
            det = np.linalg.det(np.einsum('qka,eaj->eqkj', dn, xyz[conn[start:start+128]]))
            minimum = min(minimum, float(det.min()))
            integral += float(np.sum(det*weights))
        if minimum <= 0:
            raise ValueError('Nonpositive geometry')
        min_det[str(order)] = minimum
        if order == 3:
            volume = integral
    axis, zcorn, annulus_index = profiles(level, 0.)
    diagonal = profiles(level, math.pi/4)[0]
    exact_volume = (4*W*W-math.pi*R*R)*TH
    metadata = dict(scope='BALANCED_BULK_ANNULUS_JOINT_REFINEMENT', level=level,
        domain_symmetry='FULL', counts=dict(circumferential=nt, radial=nr, thickness=nz),
        nodes=len(coords), elements=len(elements), annulus_outer_radius_m=R+.5e-6,
        annulus_corner_index=annulus_index,
        profiles=dict(radial_axis_corners_m=axis, radial_diagonal_corners_m=diagonal, z_corners_m=zcorn),
        minimum_sampled_jacobian=min_det, volume_3x3x3_m3=volume,
        exact_volume_m3=exact_volume, volume_relative_error=abs(volume-exact_volume)/exact_volume,
        maximum_radial_axis_step_m=float(np.diff(axis).max()),
        maximum_z_step_m=float(np.diff(zcorn).max()),
        first_bulk_axis_step_m=axis[annulus_index+1]-axis[annulus_index],
        status='GEOMETRY_ONLY_NOT_NUMERICAL_ACCEPTANCE')
    return coords, elements, metadata


def write(directory, level):
    coords, elements, metadata = generate(level)
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    mesh = directory/'H1.mesh'
    write_mesh(mesh, coords, elements)
    metadata.update(mesh_sha256=sha(mesh), generator_sha256=sha(Path(__file__)))
    with (directory/'manifest.json').open('x', encoding='utf-8') as stream:
        json.dump(metadata, stream, indent=2, allow_nan=False)
    return mesh, metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--level', choices=(0,1,2), type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    _, metadata = write(args.output, args.level)
    print(json.dumps(metadata, indent=2))
