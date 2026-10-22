## Material Passport

- Artifact ID: TPS-V2-G00-G03-20260901T131345Z_G00_G03
- Artifact type: deterministic code verification
- Verification Status: VERIFIED
- Solver: TPSolverV2 0.2.0-g03
- Theory baseline: 2.1
- Architecture: TPS-V2-ARCH-001

## Executed checks

- G00 source/executable manifest: present
- Native input inspection: True
- V01 HEX20 geometry and interpolation: True
- G02 material-point response and tangent: True
- G03 element residual and Jacobian: True
- V01 exact-repeat comparison: True
- G02 exact-repeat comparison: True
- G03 exact-repeat comparison: True
- Visual Studio 2022 Debug x64 rebuild: 0 errors, 0 warnings
- IDE-built V01 exit code: 0
- IDE-built G02 exit code: 0
- IDE-built G03 exit code: 0
- Legacy source/project/EXE hashes unchanged: True

## Verified numerical levels

- V01 interpolation, partition, gradient, quadrature, Jacobian and volume errors: floating-point level, maximum `1.78e-15`.
- G02 material-point tangent group errors: maximum `1.79e-11`; equilibrium neutrality, mass action and Einstein checks passed.
- G03 five residual-row tangent errors: maximum `1.18e-10`; electron/hole constant-test balances were below `8e-17`.

## Scope boundary

This verdict covers G00-G03 only. Global CSR assembly, boundary-condition elimination, PARDISO, Newton, continuation, analytic M0 cases, and production hole cases are not yet verified.

## Verdict

VERIFIED
