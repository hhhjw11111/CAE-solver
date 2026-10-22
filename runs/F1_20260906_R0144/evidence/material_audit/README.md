# ZnO material input audit

This directory is **not a production material database**. Its candidate contract
is not registered in the Fortran production entry point.

Read `../docs/zno_material_input_audit_v2.md` for the Chinese audit and the
distinction between numerical consistency and experimental calibration.

## Contents

- `zno_300k_candidate_v2.json`: explicit units, constitutive conventions, candidate
  values, source tiers, open calibration issues and protected-file hashes.
- `sources/wiley_zno_chapter1_2009.pdf`: public publisher sample chapter, unchanged.
- `evidence_20260905_v2/material_audit.json`: numerical audit, derived quantities,
  negative controls, bibliographic metadata and provenance.
- `evidence_20260905_v2/reference_pages/`: selected source pages rendered for review.

## Repeat

With Python and NumPy available, from the repository root:

```text
python solver_v2/scripts/audit_zno_material.py
```

This is offline and read-only. To save a separate evidence run:

```text
python solver_v2/scripts/audit_zno_material.py --output solver_v2/material_audit/evidence_repeat
```

The output directory must not exist. `--refresh-metadata` optionally checks the
public DOIs against Crossref; `--render-evidence` requires pypdfium2 and the public
source snapshot. No manuscript, private notes or solver data are uploaded.

A numerical `PASS` is not an experimental-calibration pass. The audit deliberately
keeps `experimental_calibration_status=PENDING` and `production_enabled=false`.
It does not execute or prove the Fortran material-point or finite-element solver.
