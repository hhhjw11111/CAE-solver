## Material Passport

- Artifact ID: TPS-V2-G00-G07-L01-20260902T025306Z_G00_G07_L01
- Artifact type: deterministic code verification
- Verification Status: VERIFIED
- Solver: TPSolverV2 0.7.0-g07-newton
- Theory baseline: 2.1
- Architecture: TPS-V2-ARCH-001

## Executed checks

- G00 source/executable manifest: present
- Native input inspection: True
- V01 HEX20 geometry and interpolation: True
- G02 material-point response and tangent: True
- G03 element residual and Jacobian: True
- G04 global CSR, scatter, Dirichlet and Jacobian: True
- L01 static PARDISO, scaling and manufactured global solve: True
- G05 M0 V02-V06 analytic solutions and convergence: True
- G05 error integration independent of solver quadrature: 5x5x5 Gauss
- G06 M2 full global 5x7 block-Jacobian audit: True
- G07 true-residual safeguarded Newton and discrete conservation: True
- V01 exact-repeat comparison: True
- G02 exact-repeat comparison: True
- G03 exact-repeat comparison: True
- G04 exact-repeat comparison: True
- L01 exact-repeat comparison: True
- G05 exact-repeat comparison: True
- G06 exact-repeat comparison: True
- G07 exact-repeat comparison: True
- Visual Studio rebuild has 0 errors and 0 warnings: True
- CLI and Visual Studio outputs match exactly: True

## Verdict

VERIFIED
