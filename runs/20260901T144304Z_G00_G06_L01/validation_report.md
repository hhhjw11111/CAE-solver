## Material Passport

- Artifact ID: TPS-V2-G00-G06-L01-20260901T144304Z_G00_G06_L01
- Artifact type: deterministic code verification
- Verification Status: VERIFIED
- Solver: TPSolverV2 0.6.0-g06-jacobian
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
- V01 exact-repeat comparison: True
- G02 exact-repeat comparison: True
- G03 exact-repeat comparison: True
- G04 exact-repeat comparison: True
- L01 exact-repeat comparison: True
- G05 exact-repeat comparison: True
- G06 exact-repeat comparison: True
- Visual Studio rebuild succeeded: True
- CLI and Visual Studio outputs match exactly: True

## Verdict

VERIFIED
