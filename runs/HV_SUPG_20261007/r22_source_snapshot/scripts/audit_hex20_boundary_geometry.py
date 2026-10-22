"""Independent HEX20 Jacobian sampling; optional new radial-midpoint candidate.

Never changes the input mesh. Candidate output directories must not exist.
Sampling is a necessary geometry check, not a proof of positivity everywhere.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

NATURAL = np.array([
    [-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
    [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1],
    [0,-1,-1],[1,0,-1],[0,1,-1],[-1,0,-1],
    [0,-1,1],[1,0,1],[0,1,1],[-1,0,1],
    [-1,-1,0],[1,-1,0],[1,1,0],[-1,1,0],
], dtype=float)


def read_mesh(path):
    with path.open(encoding="utf-8-sig") as handle:
        nn, ne = map(int, handle.readline().split())
        coords = np.empty((nn, 3))
        for node in range(nn):
            row = handle.readline().split()
            if int(row[0]) != node + 1:
                raise ValueError("Non-sequential node IDs")
            coords[node] = list(map(float, row[1:]))
        connectivity = np.empty((ne, 20), dtype=int)
        materials = np.empty(ne, dtype=int)
        for element in range(ne):
            row = list(map(int, handle.readline().split()))
            if len(row) != 22 or row[0] != element + 1:
                raise ValueError("Invalid HEX20 element record")
            materials[element] = row[1]
            connectivity[element] = np.array(row[2:]) - 1
    return coords, connectivity, materials


def derivatives(points):
    """Independent analytic derivatives, with point and node vectorization."""
    points = np.asarray(points)
    factors = 1 + points[:, None, :] * NATURAL[None, :, :]
    result = np.zeros((len(points), 3, 20))
    linear = np.sum(points[:, None, :] * NATURAL[None, :, :], axis=2) - 2
    for axis in range(3):
        others = [a for a in range(3) if a != axis]
        result[:, axis, :8] = (
            .125 * NATURAL[None, :8, axis]
            * factors[:, :8, others[0]] * factors[:, :8, others[1]]
            * (linear[:, :8] + factors[:, :8, axis])
        )
    for node in range(8, 20):
        zero = int(np.flatnonzero(NATURAL[node] == 0)[0])
        others = [a for a in range(3) if a != zero]
        result[:, zero, node] = (
            -.5 * points[:, zero]
            * factors[:, node, others[0]] * factors[:, node, others[1]]
        )
        for axis in others:
            remaining = next(a for a in others if a != axis)
            result[:, axis, node] = (
                .25 * (1 - points[:, zero] ** 2)
                * NATURAL[node, axis] * factors[:, node, remaining]
            )
    return result


def sample(coords, connectivity, points):
    d = derivatives(points)
    minimum = float("inf")
    minimum_ratio = float("inf")
    invalid_elements = []
    count = 0
    worst = None
    for element, nodes in enumerate(connectivity):
        jac = np.einsum("pkn,nj->pkj", d, coords[nodes])
        det = np.linalg.det(jac)
        scale = np.max(np.abs(jac), axis=(1, 2))
        ratio = det / np.maximum(scale ** 3, np.finfo(float).tiny)
        invalid = (~np.isfinite(det)) | (ratio <= 1e-12)
        if invalid.any():
            invalid_elements.append(element + 1)
        count += int(invalid.sum())
        q = int(np.argmin(det))
        if det[q] < minimum:
            minimum = float(det[q])
            worst = {"element": element + 1, "natural_point": list(points[q])}
        minimum_ratio = min(minimum_ratio, float(ratio.min()))
    return {
        "status": "FAIL" if count else "PASS_SAMPLED",
        "invalid_points": count, "invalid_elements": invalid_elements,
        "minimum_det_m3": minimum, "minimum_det_over_max_jacobian_cubed": minimum_ratio,
        "worst_location": worst, "points_per_element": len(points),
    }


def audit(coords, connectivity):
    gp3 = np.polynomial.legendre.leggauss(3)[0]
    gp5 = np.polynomial.legendre.leggauss(5)[0]
    faces = []
    near = []
    for axis in range(3):
        others = [a for a in range(3) if a != axis]
        for sign in [-1., 1.]:
            for u, v in itertools.product(gp5, repeat=2):
                pt = np.zeros(3)
                pt[axis] = sign
                pt[others] = [u, v]
                faces.append(pt.copy())
                pt[axis] = sign * .999
                near.append(pt.copy())
    result = {
        "volume_3": sample(coords, connectivity, list(itertools.product(gp3, repeat=3))),
        "volume_5": sample(coords, connectivity, list(itertools.product(gp5, repeat=3))),
        "all_faces_5": sample(coords, connectivity, faces),
        "near_faces_0_999": sample(coords, connectivity, near),
        "closed_grid_9": sample(coords, connectivity, list(itertools.product(np.linspace(-1, 1, 9), repeat=3))),
    }
    return result


def candidate_midpoints(coords, connectivity):
    """Only for the scoped meshes whose local r edges are radial."""
    corrected = coords.copy()
    proposed = {}
    for nodes in connectivity:
        for mid, a, b in [(8, 0, 1), (10, 3, 2), (12, 4, 5), (14, 7, 6)]:
            identifier = nodes[mid]
            position = .5 * (coords[nodes[a]] + coords[nodes[b]])
            if identifier in proposed and not np.allclose(proposed[identifier], position, rtol=0, atol=1e-21):
                raise ValueError("Inconsistent shared-edge midpoint proposal")
            proposed[identifier] = position
    for identifier, position in proposed.items():
        corrected[identifier] = position
    return corrected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--candidate-dir", type=Path)
    args = parser.parse_args()
    coords, connectivity, materials = read_mesh(args.mesh)
    result = {
        "schema": "hex20-boundary-geometry-v1",
        "mesh_path": str(args.mesh.resolve()),
        "mesh_sha256": hashlib.sha256(args.mesh.read_bytes()).hexdigest(),
        "node_count": len(coords), "element_count": len(connectivity),
        "sampling_not_global_positivity_proof": True,
        "original": audit(coords, connectivity),
    }
    if args.candidate_dir:
        if args.candidate_dir.exists():
            raise FileExistsError(args.candidate_dir)
        corrected = candidate_midpoints(coords, connectivity)
        result["candidate"] = audit(corrected, connectivity)
        result["maximum_node_change_m"] = float(np.linalg.norm(corrected - coords, axis=1).max())
        result["changed_nodes"] = int(np.count_nonzero(np.linalg.norm(corrected - coords, axis=1) > 1e-21))
        if any(check["status"] != "PASS_SAMPLED" for check in result["candidate"].values()):
            raise ValueError("Candidate geometry failed")
        args.candidate_dir.mkdir(parents=True)
        mesh = args.candidate_dir / "mesh.dat"
        with mesh.open("w", encoding="ascii", newline="\n") as handle:
            handle.write(f"{len(coords)} {len(connectivity)}\n")
            for node, xyz in enumerate(corrected, 1):
                handle.write(f"{node} " + " ".join(f"{x:.16e}" for x in xyz) + "\n")
            for element, (mat, nodes) in enumerate(zip(materials, connectivity), 1):
                handle.write(f"{element} {mat} " + " ".join(str(n + 1) for n in nodes) + "\n")
        case = args.mesh.with_name("case.nml").read_text(encoding="utf-8-sig")
        case = case.replace("G10_H02_R090_", "G14_CANDIDATE_R090_")
        (args.candidate_dir / "case.nml").write_text(case, encoding="utf-8")
        result["candidate_mesh_sha256"] = hashlib.sha256(mesh.read_bytes()).hexdigest()
    args.report.parent.mkdir(parents=True, exist_ok=True)
    if args.report.exists():
        raise FileExistsError(args.report)
    args.report.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
