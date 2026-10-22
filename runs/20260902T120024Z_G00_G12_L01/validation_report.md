## Material Passport

- Artifact ID: TPS-V2-G00-G12-L01-20260902T120024Z_G00_G12_L01
- Artifact type: deterministic code verification
- Verification Status: VERIFIED
- Solver: TPSolverV2 0.12.0-g12-radius-matrix
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
- G09 matched plain/hole geometry, topology, Jacobian and continuation: True
- G10 three-level hole-edge mesh convergence with independent 5x5x3 metrics: True
- G11 matched-boundary H00/H02 comparison with 8x8 cut-cell integration and electron, hole, thermal and displacement metrics: True
- G12 60/90/120 nm radius matrix with matched H00 controls: True
- G08 reference production case execution: True
- G09 H00 production case execution: True
- G09 H02 production case execution: True
- All three production results exact-repeat: True
- V01 exact-repeat comparison: True
- G02 exact-repeat comparison: True
- G03 exact-repeat comparison: True
- G04 exact-repeat comparison: True
- L01 exact-repeat comparison: True
- G05 exact-repeat comparison: True
- G06 exact-repeat comparison: True
- G07 exact-repeat comparison: True
- G08 exact-repeat comparison: True
- G09 exact-repeat comparison: True
- G10 exact-repeat comparison: True
- G11 exact-repeat comparison: True
- G12 exact-repeat comparison: True
- Visual Studio rebuild has 0 errors and 0 warnings: True
- CLI and Visual Studio outputs match exactly: True
- G08 reference CLI/Visual Studio production result match: True
- G09 H00 CLI/Visual Studio production result match: True
- G09 H02 CLI/Visual Studio production result match: True

## Verdict

VERIFIED
