## Material Passport

- Artifact ID: TPS-V2-G13-ENDPOINTS-20260902T132814Z_G13_ENDPOINTS
- Artifact type: deterministic endpoint mesh-sensitivity verification
- Verification Status: VERIFIED
- Solver: TPSolverV2 0.12.0-g12-radius-matrix
- Frozen protocol: g13_endpoint_mesh_sensitivity_spec_v1.md
- Executable SHA-256: $exeHash

## Executed checks

- R060 three-level G10 convergence: True
- R060 exact-repeat output: True
- R120 three-level G10 convergence: True
- R120 exact-repeat output: True
- Overall verdict: VERIFIED

## Interpretation boundary

G13 reuses the frozen G10 convergence thresholds for the two G12 endpoint radii. It establishes endpoint mesh sensitivity under the current discretization, not experimental validity, a monotonic radius law, or publication-level material calibration. The G00–G12 cross-build gate already verifies the solver executable; this G13 artifact is the endpoint mesh-sensitivity supplement.
