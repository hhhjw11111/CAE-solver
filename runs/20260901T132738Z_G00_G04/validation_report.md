## Material Passport

- Artifact ID: TPS-V2-G00-G04-20260901T132738Z_G00_G04
- Artifact type: deterministic code verification
- Verification Status: VERIFIED
- Solver: TPSolverV2 0.3.0-g04
- Theory baseline: 2.1
- Architecture: TPS-V2-ARCH-001

## Executed checks

- G00 source/executable manifest: present
- Native input inspection: True
- V01 HEX20 geometry and interpolation: True
- G02 material-point response and tangent: True
- G03 element residual and Jacobian: True
- G04 global CSR, scatter, Dirichlet and Jacobian: True
- V01 exact-repeat comparison: True
- G02 exact-repeat comparison: True
- G03 exact-repeat comparison: True
- G04 exact-repeat comparison: True
- Visual Studio 2022 Debug x64 rebuild: 0 errors, 0 warnings
- IDE-built V01/G02/G03/G04 exit codes: 0/0/0/0
- Legacy source/project/EXE hashes unchanged: True

## Verified numerical levels

- V01 maximum interpolation/geometry/quadrature error: `1.78e-15`.
- G02 maximum material-point tangent group error: `1.79e-11`.
- G03 maximum element residual-row tangent error: `1.18e-10`.
- G04 maximum global residual-row tangent error: `1.68e-10`.
- G04 CSR structure, incremental Dirichlet rows and shared-node scatter errors: exactly zero in the deterministic test.
- G04 unconstrained electron/hole balance errors: `2.04e-17` and `2.06e-16`.

## Scope boundary

This verdict covers G00-G04. It does not yet validate PARDISO, scaling, Newton, line search, continuation, analytic M0 solution cases, natural boundary loads, C04, or production hole results.

## Verdict

VERIFIED
