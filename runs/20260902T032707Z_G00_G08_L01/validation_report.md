## Material Passport

- Artifact ID: TPS-V2-G00-G08-L01-20260902T032707Z_G00_G08_L01
- Artifact type: deterministic code verification
- Verification Status: VERIFIED
- Solver: TPSolverV2 0.8.0-g08-continuation
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
- G08 voltage/Joule continuation, cutback rollback and restart: True
- G08 production case execution: True
- G08 production result exact-repeat comparison: True
- V01 exact-repeat comparison: True
- G02 exact-repeat comparison: True
- G03 exact-repeat comparison: True
- G04 exact-repeat comparison: True
- L01 exact-repeat comparison: True
- G05 exact-repeat comparison: True
- G06 exact-repeat comparison: True
- G07 exact-repeat comparison: True
- G08 exact-repeat comparison: True
- Visual Studio rebuild has 0 errors and 0 warnings: True
- CLI and Visual Studio outputs match exactly: True
- CLI and Visual Studio production results match exactly: True

## Verdict

VERIFIED
