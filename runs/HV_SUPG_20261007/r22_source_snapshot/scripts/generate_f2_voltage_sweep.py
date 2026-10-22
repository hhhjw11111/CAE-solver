"""Generate the frozen H1-NT L0 voltage sweep inputs without touching F1 evidence."""
from __future__ import annotations

import json
import re
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "cases" / "F1_FAST_HOLE_L0"
DEST = ROOT / "cases" / "F2_H1_VOLTAGE_SWEEP_20260908"
VOLTAGES = (0.5, 1.0, 1.5, 2.0, 2.5, 3.5, 4.0, 4.5, 5.0)


def tag(voltage: float) -> str:
    return f"V{round(voltage * 10):03d}"


def main() -> None:
    if DEST.exists():
        raise FileExistsError(f"Refusing to overwrite existing input directory: {DEST}")
    DEST.mkdir(parents=True)
    shutil.copy2(SOURCE / "H1_L0.mesh", DEST / "H1_L0.mesh")
    template = (SOURCE / "FAST_H1_NT_L0.nml").read_text(encoding="utf-8")
    cases = []
    for voltage in VOLTAGES:
        case_id = f"F2_H1_NT_L0_{tag(voltage)}"
        text = re.sub(r"case_id='[^']+'", f"case_id='{case_id}'", template, count=1)
        text = re.sub(
            r"target_voltage=-?[0-9.Ee+\-]+",
            f"target_voltage={-voltage:.17e}",
            text,
            count=1,
        )
        text = re.sub(
            r"field_scales=1e-10,1e-10,1e-10,[0-9.]+,",
            f"field_scales=1e-10,1e-10,1e-10,{max(voltage, 0.5):.6g},",
            text,
            count=1,
        )
        path = DEST / f"{case_id}.nml"
        path.write_text(text, encoding="utf-8")
        cases.append({"case_id": case_id, "voltage_V": voltage, "input": path.name})
    manifest = {
        "study": "F2_H1_NT_L0_voltage_sweep",
        "geometry": "20um x 20um x 1um, centered insulated hole r=3um",
        "model": "NT",
        "mesh_source": str((SOURCE / "H1_L0.mesh").relative_to(ROOT)),
        "omitted_existing_voltage_V": 3.0,
        "zero_voltage": "exact zero state; no solver launch",
        "cases": cases,
    }
    (DEST / "inputs_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"F2_INPUTS_GENERATED={len(cases)}")
    print(f"F2_INPUT_ROOT={DEST}")


if __name__ == "__main__":
    main()
