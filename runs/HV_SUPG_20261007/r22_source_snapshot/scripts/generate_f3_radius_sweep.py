"""Generate the user-authorized ZnO hole-radius sweep without overwriting evidence.

The 20 x 20 x 1 um specimen, 3 V excitation, 0.5 um physical sampling ring,
and L0 discretization policy are held fixed. Only the centered hole radius and
the LQ/N0/NT model branch vary. R=0 and R=3 um are reused from the preceding F1
matrix, so this generator creates only the new R=2 and R=4 um solves.
"""
from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "cases" / "F3_RADIUS_SWEEP_20260908"
PATTERN = (
    (0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0),
    (0, 0, 2), (2, 0, 2), (2, 2, 2), (0, 2, 2),
    (1, 0, 0), (2, 1, 0), (1, 2, 0), (0, 1, 0),
    (1, 0, 2), (2, 1, 2), (1, 2, 2), (0, 1, 2),
    (0, 0, 1), (2, 0, 1), (2, 2, 1), (0, 2, 1),
)
HALF_WIDTH = 10.0e-6
THICKNESS = 1.0e-6
FIRST_LAYER = 14.0e-9
RING_WIDTH = 0.5e-6
RADII_UM = (2, 4)
MODELS = ("LQ", "N0", "NT")


def graded_skin(first: float, depth: float) -> list[float]:
    values = [0.0]
    step = first
    while values[-1] + 2.0 * step < depth:
        values.append(values[-1] + step)
        step *= 2.0
    values.append(depth)
    return values


def mids(corners: list[float]) -> list[float]:
    values: list[float] = []
    for left, right in zip(corners, corners[1:]):
        values.extend((left, 0.5 * (left + right)))
    return values + [corners[-1]]


def write_mesh(path: Path, coordinates: list[tuple[float, float, float]], elements: list[list[int]]) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(f"{len(coordinates)} {len(elements)}\n")
        for index, xyz in enumerate(coordinates, 1):
            stream.write(f"{index} " + " ".join(f"{value:.17e}" for value in xyz) + "\n")
        for index, connectivity in enumerate(elements, 1):
            stream.write(f"{index} 1 " + " ".join(map(str, connectivity)) + "\n")


def generate_mesh(radius: float) -> dict[str, object]:
    circumferential = 16
    skin = graded_skin(FIRST_LAYER, RING_WIDTH)
    zskin = graded_skin(FIRST_LAYER, 0.5 * THICKNESS)
    zcorners = [-0.5 * THICKNESS + value for value in zskin]
    zcorners += [0.5 * THICKNESS - value for value in reversed(zskin[:-1])]
    zvalues = mids(zcorners)
    thickness_elements = len(zcorners) - 1
    radial_elements = 2 * (len(skin) - 1) + 4
    identifiers: dict[tuple[int, int, int], int] = {}
    coordinates: list[tuple[float, float, float]] = []
    elements: list[list[int]] = []

    def radial_corners(angle: float) -> list[float]:
        outer = HALF_WIDTH / max(abs(math.cos(angle)), abs(math.sin(angle)))
        inner_bulk = radius + skin[-1]
        outer_bulk = outer - skin[-1]
        return (
            [radius + value for value in skin]
            + [inner_bulk + (outer_bulk - inner_bulk) * i / 4.0 for i in range(1, 4)]
            + [outer - value for value in reversed(skin)]
        )

    for ez in range(thickness_elements):
        for et in range(circumferential):
            for er in range(radial_elements):
                connectivity: list[int] = []
                for pr, pt, pz in PATTERN:
                    ir = 2 * er + pr
                    it = (2 * et + pt) % (2 * circumferential)
                    iz = 2 * ez + pz
                    key = (ir, it, iz)
                    if key not in identifiers:
                        angle = math.pi * it / circumferential
                        radius_at_node = mids(radial_corners(angle))[ir]
                        identifiers[key] = len(coordinates) + 1
                        coordinates.append(
                            (radius_at_node * math.cos(angle), radius_at_node * math.sin(angle), zvalues[iz])
                        )
                    connectivity.append(identifiers[key])
                elements.append(connectivity)

    tag = f"R{round(radius * 1.0e6):02d}"
    mesh_name = f"F3_{tag}_L0.mesh"
    write_mesh(DEST / mesh_name, coordinates, elements)
    return {
        "radius_m": radius,
        "radius_over_width": radius / (2.0 * HALF_WIDTH),
        "mesh": mesh_name,
        "nodes": len(coordinates),
        "elements": len(elements),
        "first_layer_m": FIRST_LAYER,
        "ring_width_m": RING_WIDTH,
        "circumferential_elements": circumferential,
        "radial_elements": radial_elements,
        "thickness_elements": thickness_elements,
    }


def write_case(mesh: dict[str, object], model: str) -> dict[str, object]:
    radius = float(mesh["radius_m"])
    tag = f"R{round(radius * 1.0e6):02d}"
    case_id = f"F3_{tag}_{model}_L0"
    text = f"""&case_control
 schema_version=2, case_id='{case_id}', model_level='M2_THERMO_TWO_CARRIER',
 mesh_file='{mesh['mesh']}'
/
&problem_control
 problem_type='CENTER_HOLE_PLATE',
 material_profile='ZNO_F1_{model}', target_voltage=-3.00000000000000000e+00,
 target_voltage_factor=1.0, target_joule_factor=1.0, hole_radius={radius:.17e}
/
&solver_control
 field_scales=1e-10,1e-10,1e-10,3.0,0.1,0.1,10.0,
 max_newton=40, residual_tolerance=1e-9, update_tolerance=1e-9,
 continuation_initial_step=0.25, continuation_minimum_step=0.015625,
 continuation_maximum_step=0.5, continuation_maximum_attempts=100,
 continuation_maximum_newton_per_attempt=30
/
"""
    filename = f"{case_id}.nml"
    (DEST / filename).write_text(text, encoding="ascii")
    return {"case_id": case_id, "radius_m": radius, "model": model, "input": filename}


def main() -> None:
    if DEST.exists():
        raise FileExistsError(f"Refusing to overwrite existing input directory: {DEST}")
    DEST.mkdir(parents=True)
    meshes = [generate_mesh(radius_um * 1.0e-6) for radius_um in RADII_UM]
    cases = [write_case(mesh, model) for mesh in meshes for model in MODELS]
    manifest = {
        "study": "F3_ZnO_center_hole_radius_sweep",
        "evidence_level": "L0_TREND_ONLY_NOT_MESH_ACCEPTANCE",
        "geometry": {"width_m": 20.0e-6, "height_m": 20.0e-6, "thickness_m": THICKNESS},
        "sweep": {"all_radii_um": [0, 2, 3, 4], "newly_solved_radii_um": list(RADII_UM)},
        "excitation": {"left_phi_V": 0.0, "right_phi_V": -3.0, "nominal_field_V_m": 1.5e5},
        "contacts": "left/right eta_n=eta_p=theta=0; T=300 K",
        "mechanics": "minimal rigid-body anchors; otherwise zero natural traction",
        "other_surfaces": "zero natural electric-displacement, carrier-current and heat-flux terms",
        "hole_wall": "ideal insulating natural boundary; exterior dielectric and surface charge omitted",
        "models": {
            "LQ": "sequential linearized reference; no return feedback",
            "N0": "full-density nonlinear model with mu and D frozen at 300 K",
            "NT": "full-density nonlinear model with mu(T) and D(T)",
        },
        "reuse": {
            "R0": "runs/F1_FAST_HOLE_L0_20260907/FAST_H0_{model}_L0",
            "R3": "runs/F1_FAST_HOLE_L0_20260907/FAST_H1_{model}_L0",
        },
        "meshes": meshes,
        "new_cases": cases,
    }
    (DEST / "inputs_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"F3_INPUTS_GENERATED={len(cases)}")
    print(f"F3_INPUT_ROOT={DEST}")


if __name__ == "__main__":
    main()
