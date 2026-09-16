# TPSolverV2 refactor status — 2026-09-02

## Current trusted checkpoint

Version: `0.13.2-radius-control-audit`

Theory: `thermo_piezo_semiconductor_two_carrier_theory_v2_1`

Architecture: `TPS-V2-ARCH-001`

| Gate | Status | Evidence |
|---|---|---|
| G00 reproducible build/provenance | PASS | clean ifx build, build manifest, source and EXE SHA-256 |
| G01 HEX20 geometry/interpolation | PASS | V01, maximum error `1.78e-15`; Jacobian degeneracy check is scale-relative for SI microstructures |
| G02 material-point physics/tangent | PASS | maximum group directional-derivative error `1.79e-11` |
| G03 element residual/Jacobian | PASS | maximum five-row directional-derivative error `1.18e-10` |
| G04 global assembly/boundary conditions | PASS | two shared HEX20 elements; maximum global tangent error `1.68e-10` |
| L01 linear solve prerequisite | PASS | static sequential PARDISO plus seven-field scaling; 224-equation manufactured solve error `1.60e-12` |
| G05 M0 analytic solutions | PASS | V02/V06: approximately third-order L2 and second-order H1; V03/V04 zero-current signs; V05 heat balance `1.63e-15` |
| G06 global Jacobian | PASS | full M2 5x7 block audit; maximum error `4.90e-12`; all active blocks stable for 9/9 steps; structural-zero leakage `0` |
| G07 Newton/conservation | PASS | warm start 5 iterations, cold start 9 iterations; final residuals about `1.7e-16`; solution difference `1.02e-13`; maximum conservation defect `1.26e-15` |
| G08 continuation/production path | PASS | main, cutback-stress and restart paths; 2 real cutbacks with exact rollback; path difference `3.72e-14`; maximum balance defect `2.28e-16` |
| G09 center-hole geometry/path | PASS | matched H00/H02 outer boundaries; 32-element r=90 nm hole mesh; hole-wall BC overlap `0`; volume error `3.75e-6`; global Jacobian error `3.02e-12`; balance defect `4.08e-16` |
| G10 center-hole mesh convergence | PASS after 2026-09-06 geometry repair | Repaired 32/72/128-element graded sequence; all six frozen metrics pass, with L1→L2 changes: total Joule 0.282%, annulus Joule mean 0.205%, p95 0.346%, electric-field mean 0.587%, electron redistribution mean 0.780%, temperature-rise mean 0.248%; 5/6 changes decrease |
| G11 matched H00/H02 control | PASS after repair; exact repeat | 0.13.2 repaired core evidence; old physical metrics remain superseded |
| G12 radius matrix | PASS after repair; exact repeat | 0.13.2; all three H00 controls checked; 27 printed contrasts independently consistent |
| G13 endpoint mesh sensitivity | PASS after repair; exact repeat | 0.13.1; R060/R120 both pass; relevant production/G10 source matches 0.13.2 |
| A1J production Joule/flux audit | PASS | 13 synthetic production-path cases; independent 4/5-order flux, field/temperature analytic check, dual-start/continuation and negative controls |
| G12C per-radius control guard | PASS | 90 negative controls; independent output-audit self-test adds 119 negative controls |

## 2026-09-06 geometry correction

An independent boundary-trace audit found negative HEX20 Jacobians in the old graded-hole meshes. Corrected G09-G13 CLI gates and exact-repeat checks are complete. See `../runs/20260906_REPAIRED_CORE_0132/repaired_validation_report.md`. This does not complete the ZnO material or LQ/N0/NT study; the old hole figures remain superseded.

## Implemented theory

- Numerical nodal state: `[u1,u2,u3,phi,eta_n,eta_p,theta]`.
- Gauss-point concentrations: `n=n_ref*exp(eta_n)`, `p=p_ref*exp(eta_p)`.
- Linear thermo-piezoelectric constitutive response.
- Two-carrier low-field drift-diffusion with temperature-dependent mobility and Einstein diffusion.
- Frozen `R=0` recombination interface.
- Drift-conductivity Joule heat with voltage/Joule continuation factors.
- Full material-point chain derivatives with respect to strain, electric field, temperature, temperature gradient, both log concentrations and both log-concentration gradients.
- Five weak residual rows and consistent `140x140` HEX20 element Jacobian.
- Node-pair sorted seven-field scalar CSR graph without a dense node adjacency matrix.
- Shared-node global residual/Jacobian assembly.
- State projection and Newton-increment Dirichlet unit rows.
- LP64 `mtype=11` PARDISO wrapper with checked status propagation and explicit release.
- Variable scaling and Jacobian-row residual equilibration.
- Independent five-equation solve and 224-equation seven-field manufactured global solve.
- Registered M0 linearized variables `eta_n=Delta n/n0` and `eta_p=Delta p/p0`, isolated from the M1/M2 exponential full-density model.
- V02-V06 electromechanical, electron, hole, heat, and joint thermo-piezo-electron analytic benchmarks.
- Independent 5x5x5 error quadrature and an audited 1-2-4-8-16 element convergence sequence.
- Full M2 global 5x7 block-Jacobian audit using the production residual/Jacobian assembly path.
- Frozen nine-step central-difference scan, field/residual scaling, structural-zero leakage check, and fixed-DOF perturbation check.
- Single-load-point scaled PARDISO Newton with true-residual Armijo backtracking.
- Component trust region for exponential electron/hole log-density updates with the actual merit directional derivative.
- Dual-start convergence-to-the-same-root audit plus charge, electron current, hole current and heat balance checks.
- Voltage-first, Joule-second continuation with deterministic growth/cutback and accepted-state checkpoints.
- Independent failed-trial rollback audit and restart-from-`(1,0)` equality check.
- Schema-v2 `REFERENCE_BAR` production setup and deterministic continuation/nodal result writer.
- Explicit `PLAIN_PLATE_CONTROL` and `CENTER_HOLE_PLATE` profiles with declared hole radius.
- Deterministic matched H00/H02 mesh generation and exterior/internal/hole face topology audit.
- Hole-wall natural-boundary proof, distorted-mesh global Jacobian audit, and full continuation solve.
- Three-level deterministic hole-edge refinement and independent 5x5x3 field re-evaluation.
- Frozen global Joule-power and 30 nm hole-annulus weighted mean/p95 convergence criteria.
- Matched-boundary H00-L2 mesh and common physical-annulus H00/H02 contrast.
- Independent electron, hole, temperature-rise and displacement metrics for the matched comparison.
- G11-v2 8x8 cut-cell subdivision, replacing the v1 local metrics after its H00 annulus-volume diagnostic exposed 7.89% integration error.
- G12 fixed-voltage radius matrix with matched H00 virtual-annulus controls at 60/90/120 nm.

## Deliberately unavailable

The `run` command is available for the frozen schema-v2 reference bar, plain control plate, and
center-hole plate profiles with the `TPS_G08_REFERENCE` material. The standalone `jacobian`
command remains unavailable; G06 and G09 are the trusted Jacobian audits. G10 proves mesh stability
of the six frozen integral/percentile observables for the r=90 nm, -5 V case. G11 proves a controlled
H00/H02 numerical contrast at that point. G12 exercises a 60/90/120 nm matrix, and G13 verifies
the two endpoint mesh sensitivities; experimental validity remains open. Neither gate proves raw point maxima, C04
convergence, or publication-level quantitative material predictions.

## Requested validation objective

**REOPENED / NOT CURRENTLY VERIFIED as an integrated hole-result package.**
The 2026-09-06 geometry audit invalidated the old G10–G13 hole evidence.
Corrected-mesh double runs and output audits are complete for the declared CLI gates. The former integrated report
`../final_results/final_validation_report.md` is historical, not a current completion claim.
See `geometry_audit_findings_20260906.md` and `g12_radius_control_audit_20260906.md`.

Geometry/removed-volume decomposition, C04 voltage-path physics, experimental calibration,
nonzero recombination and flexoelectricity are separate future research units and are not required
to close the present correctness-validation objective.

## Evidence directory

The latest G00-G12 plus A1J and L01 evidence directory is generated by `scripts/run_verification.ps1`; G13 endpoint evidence remains in its separate run directory. Earlier evidence remains immutable in its original run directories.

Use each run's report to determine its actual coverage. A repaired CLI double run is not
automatically a new Visual Studio cross-build comparison. Incomplete evidence directories
must not be called VERIFIED.
