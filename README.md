# TPSolverV2

> **当前唯一程序基准（2026-09-08 用户指定）**：旧 Console3 已废弃。后续所有验证、程序研究和流程图以 `solver_v2/` 当前源码为准。当前执行范围及结果状态优先读取 [current_execution_plan.md](docs/current_execution_plan.md)；本文件下方旧阶段清单保留为历史说明。废弃范围见 [Console3 废弃说明](../LEGACY_CONSOLE3_DEPRECATED.md)。

Word 稿件统一遵循 [项目 Word 格式规范](docs/manuscript_word_format.md) 及其中冻结的期刊格式模板。

Current priority: obtain 6–8 V states on the existing A64 mesh. Further mesh and
integration comparisons are deferred at the user's request; their old queue was
stopped while the high-voltage solve remains running. Original model and gates
are retained. Automatic follow-up checks now continue this task after a run ends
or fails when enabled; the user has now paused the repeated follow-up and retired
the old parallel voltage route. Only `r26_log_A64_trust_heat_restart_V080` is active.
Preliminary reporting separates obtained full-load fields from failed
flux checks without certifying mesh accuracy. See the top of the current execution
plan for authoritative active run paths and status.

`TPSolverV2` is the clean implementation of the frozen two-carrier thermo-piezoelectric-semiconductor theory v2.1.

The current development version is `0.16.0-trust-region`. R26C adds opt-in hybrid
Newton trust-region globalization; original equations, load targets and tolerances
are unchanged. G07TR, G07S, TCTRL, 16 byte-identical frozen CLI comparisons and
79 Python HV tests pass, including the default-off fixed-voltage heat restart controls.
R26/R26B preliminary candidates failed the new cold-start
check and are retained as failures, not used for production. The active 8 V
zero-to-full-heat run is `r26_log_A64_trust_heat_restart_V080`; full heat is still required
before calling it an 8 V result. See the current execution plan for active paths.

R25 permits the verified log-weighted form
in the existing temperature-parameterized diagnostic runner. R24 fixes a Debug
array-assignment failure in optional order-5 integration
for the new transport form. R22/R23 each pass all 16 frozen CLI regressions;
three-dimensional mesh and integration studies remain in progress.
The new log-weighted transport option is off by default and retains the original
physical currents and heat source. See the [prototype checks](docs/hv_log_weighted_probe_20261008.md).
Its optional HEX20
residual-based SUPG and continuation improvements are disabled for existing inputs.
The 0–8 V study is **not yet fully validated**; follow the current execution plan and
[implementation/acceptance record](docs/hv_0_8v_completion_plan_20261007.md).
The latest completed evidence inventory is
[R25 stage 4](runs/HV_SUPG_20261007/delivery_r25_stage4/evidence.json).
Its 0.5 V thickness-only three-level study passes both scalar/ring and local-field
checks (last-level carrier difference 0.0662% or less). Its independent angular
32/48/64 study passes at the last pair (1.098% or less), and A64 integration
orders 3/5 pass both scalar and local-field checks (local carrier difference
0.00247% or less). These are separate directional studies, not an unrun joint mesh.
At 6 V, the angular last-pair local carrier difference is at most 1.058%,
the thickness last-pair difference is at most 0.06579%, and A64 quadrature
3/5 local differences are at most 0.003270%; all corresponding gates pass.
The independent radial 14/16/18 study also passes, with last-pair carrier
differences at most 0.09367%. These remain directional studies, not a joint-mesh
or interval certificate. The 1 V accepted checkpoint has independently passed
fixed-voltage solve and physical-flux replay.
The A64 continuation from 0.5 to 6 V has completed: all 12 nonzero half-volt
checkpoints passed independent fixed-voltage solve/flux replay. The fixed-mesh
6.25 V replay also passes (24.6163 K temperature rise; 0.90081% exterior heat error).
One original-model 8 V heat continuation is now running on the fixed A64 mesh.
The independently replayed 6.375 V field converges but fails the exterior heat-flux
gate (1.11589% versus 1%); it remains explicitly failed. Further mesh comparisons
are deferred; see [running high-voltage status](runs/HV_SUPG_20261007/r26_high_voltage_priority/results_restart/status_6_8V.csv).
Temperature-parameterized branch samples are diagnostic until fixed-voltage replay.
The sampled 6.46235 V turning-region point fails the L2 physical heat-flux gate;
the 6.39300 V descending sample passes after independent S2 thickness re-solve,
but lacks a complete spatial certificate. No 8 V state is certified.
Historical result files and executables remain separate from the experimental builds.
R21 adds an opt-in y-half symmetric-branch runner, checked against full-domain
3 V and 6 V states. It separately audits physical-exterior and symmetry-cut fluxes;
all raw integrals remain actual half-domain quantities. See the
[half-domain and directional-refinement record](docs/hv_half_y_refinement_20261008.md).
The frozen historical [R21 credibility assessment](runs/HV_SUPG_20261007/delivery_r21/可信范围核验.md)
lists all 17 requested voltage states and the separately replayed turning-branch samples.
Completed 0.5 V spatial/coefficient studies and 6 V angular refinement still fail
the frozen hole-field limits. Only the exact 0 V equilibrium is fully verified;
no nonzero interval is certified. Missing 6.5–8 V states are not interpolated.
The R17 delivery remains an immutable earlier-stage record.

The experimental checks are `verify SUPG` (high-Pe analytical transport),
`verify SUPGJ` / `verify SUPGJT` (directional / trace stabilization tangents),
`verify G07S` (stagnation rollback trigger), `verify G08P` / `verify G08G`
(predictor / fixed-voltage continuation checkpoints), and `verify TCTRL`
(inverse temperature control retaining every original PDE). Existing inputs retain all
new options disabled. See the implementation record before enabling them.
The separate `trace-f1-temperature` entry is an experimental inverse-load diagnostic;
it does not replace fixed-voltage acceptance or establish dynamical stability.

The R18 linear-solver experiment (`verify BLK`) forms a hole-block candidate and
checks its componentwise backward error against the complete original matrix
(limit `1e-11`); an unsuitable candidate falls back to the original full solve.
It is disabled by default. The point/temperature runners expose `--split-holes`
and `--storage OOC --ooc-core-mb <MB>` for explicitly selected new runs. OOC uses
an isolated run-owned scratch directory; its RAM budget excludes assembly and OS
memory. See [the block-solver note](docs/hv_block_solver_note_20261008.md).
R19 additionally omits only exactly zero off-diagonal entries from the candidate
blocks. Every diagonal and every nonzero coefficient is retained; the original
full-matrix check and fallback remain in place. This optimization also stays off
for existing inputs.

For a fresh 0–8 V run, `scripts/run_hv_trajectory.py --exe <experimental-exe>
--mesh <H1-mesh> --output <new-directory>` saves actual 0.5 V checkpoints and
independently audits each solved point. The output directory must not already
exist. By default the main trajectory and its flux audits run sequentially to
bound memory use. A `SOLVE_FLUX_PASS_MESH_PENDING` point is not a mesh certificate;
the script never labels an incomplete interval as validated.

`scripts/run_hv_temperature_trace.py` requires a warm full-heat NT checkpoint
on the same mesh. It records the voltage actually recovered at each temperature;
an input target of 8 V does not make those checkpoints 8 V solutions. Recheck
selected states with `run_hv_point.py` at their recorded voltages.
`audit_hv_stabilization.py` compares SUPG scales 0.5/1/2 at one common voltage,
mesh and executable, with a predeclared 1% metric tolerance.
`run_hv_scale_path.py` reaches a different SUPG coefficient through independently
solved steps at fixed voltage, with rollback and cutback after a rejected trial.
It first replays the seed in the selected executable and preserves every attempt.
This path does not itself grant a stabilization or mesh certificate.
`audit_hv_mesh.py` accepts three sweep roots or three single-point runs; it checks
both carrier annular metrics at the original 2% limit, alongside the other frozen
metrics, and rejects mismatched voltage/model/executable evidence.
`diagnose_hv_heat_components.py` reports drift-heat components and temperature/
field ranges without changing the frozen heat model or the acceptance gates.
`summarize_hv_credibility.py --plan <evidence-plan.json> --output <new-directory>`
rechecks native point evidence and recomputes registered mesh/coefficient studies.
It retains actual half-domain integrals alongside explicitly reconstructed full-domain
values and keeps off-grid branch samples separate from the 0.5 V grid. A successfully
written report does not mean its numerical studies passed.

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
