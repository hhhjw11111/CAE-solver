# TPSolverV2

> **当前唯一程序基准（2026-09-08 用户指定）**：旧 Console3 已废弃。后续所有验证、程序研究和流程图以 `solver_v2/` 当前源码为准。当前执行范围及结果状态优先读取 [current_execution_plan.md](docs/current_execution_plan.md)；本文件下方旧阶段清单保留为历史说明。废弃范围见 [Console3 废弃说明](../LEGACY_CONSOLE3_DEPRECATED.md)。

`TPSolverV2` is the clean implementation of the frozen two-carrier thermo-piezoelectric-semiconductor theory v2.1.

Current verification scope:

- Phase 1 scaffold and versioned input;
- seven degrees of freedom per node;
- HEX20 shape functions and 3x3x3 Gauss quadrature;
- two-carrier thermo-piezoelectric material-point response and consistent tangent;
- HEX20 five-row weak residual and 140x140 consistent element Jacobian;
- seven-field scalar CSR graph, shared-node global assembly and incremental Dirichlet rows;
- statically linked sequential MKL PARDISO and explicit seven-field scaling;
- M0 V02-V06 analytic benchmarks with independent error quadrature;
- M2 global 5x7 block-Jacobian audit with a frozen multi-step finite-difference scan;
- scaled PARDISO Newton increments, log-density trust regions and true-residual Armijo backtracking;
- G07 dual-start nonlinear convergence, solution agreement and discrete conservation audit;
- voltage-then-Joule continuation with deterministic cutback, rollback and restart;
- versioned M2 reference-bar production input and deterministic nodal/result output;
- explicit plain-plate/center-hole problem profiles and boundary-topology audit;
- matched G09 H00/H02 coarse meshes with production-path Jacobian and continuation checks;
- three-level G10 hole-edge local refinement with independent 5x5x3 integral/percentile metrics;
- G11 matched-boundary H00/H02 comparison on the converged L2 resolution, with 8x8 cut-cell subdivision for the common annulus;
- G12 matched H00 controls for the 60/90/120 nm radius matrix;
- G13 independent three-level endpoint mesh sensitivity for R060 and R120;
- G00-G12 and L01 verification entry points.

The integrated verified report, publication-ready figures, source CSV files and figure hashes are
under `final_results/`.

The legacy `../Console3.f90` is not compiled into this program.

## Commands

```text
TPSolverV2.exe version
TPSolverV2.exe inspect cases\validation\V01_patch\case.nml
TPSolverV2.exe verify V01
TPSolverV2.exe verify G02
TPSolverV2.exe verify G03
TPSolverV2.exe verify G04
TPSolverV2.exe verify L01
TPSolverV2.exe verify V02
TPSolverV2.exe verify V03
TPSolverV2.exe verify V04
TPSolverV2.exe verify V05
TPSolverV2.exe verify V06
TPSolverV2.exe verify G05
TPSolverV2.exe verify G06
TPSolverV2.exe verify G07
TPSolverV2.exe verify G08
TPSolverV2.exe verify G09 cases\validation\G09_center_hole\H00_plain\case.nml cases\validation\G09_center_hole\H02_r090\case.nml
TPSolverV2.exe verify G10 cases\validation\G10_mesh_convergence\L0\case.nml cases\validation\G10_mesh_convergence\L1\case.nml cases\validation\G10_mesh_convergence\L2\case.nml
TPSolverV2.exe verify G11 cases\validation\G11_matched_control\H00_L2\case.nml cases\validation\G10_mesh_convergence\L2\case.nml
TPSolverV2.exe verify G12 cases\validation\G11_matched_control\H00_L2\case.nml cases\validation\G12_radius_matrix\R060\case.nml cases\validation\G12_radius_matrix\R090\case.nml cases\validation\G12_radius_matrix\R120\case.nml
TPSolverV2.exe run cases\validation\G08_reference_bar\case.nml result.dat
TPSolverV2.exe run cases\validation\G09_center_hole\H02_r090\case.nml hole_result.dat
```

Use `scripts\build_debug.ps1` from PowerShell to create a checked x64 Debug build with Intel oneAPI ifx.

The `run` command is released for the explicitly named M2 `REFERENCE_BAR`,
`PLAIN_PLATE_CONTROL`, and `CENTER_HOLE_PLATE` profiles with the
`TPS_G08_REFERENCE` material. G09 validates one matched coarse H00/H02 pair, including natural
hole-wall boundaries. G10 establishes mesh stability for six frozen integral/percentile metrics at
the r=90 nm, -5 V point. It does not establish radius trends, C04 physics, or convergence of raw
point maxima. G11 establishes a numerically matched H00/H02 geometry contrast at that one point,
including electron, hole, thermal and displacement observables. The G11-v2 common-annulus volume
errors are below 1% for both geometries. G12 reports the frozen 60/90/120 nm numerical radius
matrix, and G13 verifies the two endpoint mesh sensitivities. These results do not establish
experimental validity, a universal monotonic radius law, C04 voltage-path physics, or raw-point
maximum convergence.
Unknown production problem/material profiles fail explicitly.

## ZnO material preparation (2026-09-05)

The separate `material_audit/` package records a candidate ZnO input contract,
numerical input checks, bibliographic evidence and unresolved experimental
calibration issues. See `docs/zno_material_input_audit_v2.md`.
This does not change the trusted solver checkpoint or enable ZnO production runs.
The next approved numerical stage is a reduced analytical benchmark through the
production assembly path with independently integrated boundary fluxes.

`scripts\generate_g09_meshes.ps1` deterministically regenerates the frozen matched G09 meshes.
`scripts\generate_g10_meshes.ps1` deterministically regenerates the three frozen graded G10 meshes.
`scripts\generate_g11_meshes.ps1` regenerates the matched-resolution H00 control mesh.
`scripts\generate_g12_meshes.ps1` regenerates the frozen 60/90/120 nm radius meshes.
`scripts\generate_g13_meshes.ps1` regenerates the frozen R060/R120 endpoint refinement meshes.
`scripts\postprocess_final_results.py` regenerates the final source tables and figures without
rerunning the solver.
