"""Generate an explicitly labelled L0 six-case presentation matrix."""
from pathlib import Path
import json
import generate_f1_main_cases as base


def main():
    root = Path(__file__).resolve().parents[1] / "cases" / "F1_FAST_HOLE_L0"
    if root.exists():
        raise FileExistsError(root)
    root.mkdir(parents=True)
    base.ROOT = root
    base.generate(0, True)
    base.generate(0, False)
    cases = []
    for hole in [False, True]:
        for model in ["LQ", "N0", "NT"]:
            name = base.case(hole, 0, model, 3.0, "FAST")
            cases.append(name)
    (root/"inputs_manifest.json").write_text(json.dumps({
        "purpose": "Fast presentation matrix; not final L2 mesh acceptance.",
        "mesh_level": "L0", "models": ["LQ", "N0", "NT"],
        "geometries": ["H0", "H1"], "voltage_V": 3.0, "cases": cases,
        "same_mesh_within_geometry": True,
        "material_evidence": "NUMERICAL_CANDIDATE_NOT_EXPERIMENTALLY_CALIBRATED",
    }, indent=2), encoding="utf-8")
    print(root)


if __name__ == "__main__":
    main()
