"""New full-domain meshes for the frozen JE reconstruction study.

This extends the audited GLOBAL_DIRECTIONAL coordinate profiles without changing
the historical generator. The r=3.5 um measurement interface remains a corner.
"""
from pathlib import Path
import argparse
import json
import math
import numpy as np
from generate_f1_main_cases import PAT, R, W, TH, skin, mids, write_mesh
from generate_je_global_direction_mesh_20261016 import resample
from generate_hv_surface_mesh import audit_geometry
from audit_hex20_boundary_geometry import derivatives
from run_hv_sweep import sha


def generate(angular=40, radial=14, thickness=14):
    if angular not in (32, 40, 48) or radial not in (14, 16, 18, 20, 22) or thickness not in (14, 16, 18):
        raise ValueError('Only predeclared FULL A32/40/48, R14/16/18/20/22, Z14/16/18 are supported')
    level = (radial - 14) // 2
    layer = skin(14e-9, .5e-6)
    zskin = skin(3.5e-9, .5e-6)
    oldz = [-TH/2+d for d in zskin] + [TH/2-d for d in reversed(zskin[:-1])]
    zcorn = resample(oldz, thickness)
    zs = mids(zcorn)

    def corners(angle):
        outer = W/max(abs(math.cos(angle)), abs(math.sin(angle)))
        a, b = R+layer[-1], outer-layer[-1]
        inner = [R+d for d in layer]
        rest = [a] + [a+(b-a)*i/4 for i in range(1, 4)] + [outer-d for d in reversed(layer)]
        return resample(inner, 5+level) + resample(rest, 9+level)[1:]

    profiles = [mids(corners(math.pi*it/angular)) for it in range(2*angular)]
    ids, coords, elements = {}, [], []
    for ez in range(thickness):
        for et in range(angular):
            for er in range(radial):
                conn = []
                for pr, pt, pz in PAT:
                    ir, it, iz = 2*er+pr, (2*et+pt) % (2*angular), 2*ez+pz
                    key = ir, it, iz
                    if key not in ids:
                        angle, radius = math.pi*it/angular, profiles[it][ir]
                        ids[key] = len(coords)+1
                        coords.append((radius*math.cos(angle), radius*math.sin(angle), zs[iz]))
                    conn.append(ids[key])
                elements.append(conn)
    q, w = np.polynomial.legendre.leggauss(3)
    points = np.array([(a,b,c) for a in q for b in q for c in q])
    weights = np.array([a*b*c for a in w for b in w for c in w])
    dn = derivatives(points)
    xyz, conn = np.asarray(coords), np.asarray(elements)-1
    volume = 0.
    for begin in range(0, len(conn), 256):
        jac = np.einsum('qka,eaj->eqkj', dn, xyz[conn[begin:begin+256]])
        volume += float(np.sum(np.linalg.det(jac)*weights))
    metadata = dict(domain_symmetry='FULL', counts=dict(circumferential=angular, radial=radial, thickness=thickness),
                    nodes=len(coords), elements=len(elements), annulus_outer_radius_m=R+.5e-6,
                    annulus_corner_index=5+level, profiles=dict(radial_axis_corners_m=corners(0.),
                    radial_diagonal_corners_m=corners(math.pi/4), z_corners_m=zcorn),
                    volume_3x3x3_m3=volume, geometry_audit=audit_geometry(coords,elements),
                    status='MESH_DESIGN_NOT_NUMERICAL_ACCEPTANCE',
                    scope='Global logical-coordinate refinement; other requested directions fixed; nonnested')
    metadata['maximum_radial_axis_step_m'] = float(np.diff(corners(0.)).max())
    metadata['maximum_z_step_m'] = float(np.diff(zcorn).max())
    return coords, elements, metadata


def write(directory, angular=40, radial=14, thickness=14):
    directory = Path(directory).resolve()
    coords, elements, metadata = generate(angular, radial, thickness)
    directory.mkdir(parents=True, exist_ok=False)
    mesh = directory/'H1.mesh'
    write_mesh(mesh, coords, elements)
    metadata.update(mesh_sha256=sha(mesh), generator_sha256=sha(Path(__file__)))
    with (directory/'manifest.json').open('x', encoding='utf-8') as stream:
        json.dump(metadata, stream, indent=2, allow_nan=False)
    return mesh, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts', nargs=3, type=int, metavar=('ANGULAR','RADIAL','THICKNESS'), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    _, metadata = write(args.output, *args.counts)
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
