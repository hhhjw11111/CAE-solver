# TPSolverV2 build and provenance policy

Status: policy for the 2026-09-15 Git baseline. Version synchronization is **not yet implemented**.

## Current build facts

- Supported scripted build: `scripts/build_debug.ps1`; compiler is Intel oneAPI `ifx` (exact installed release: **not yet recorded**).
- Linkage uses Intel MKL PARDISO. Sequential and threaded MKL modes are selected by the script; exact MKL release: **not yet recorded**.
- PARDISO storage mode is runtime-selectable with `TPS_PARDISO_STORAGE=IC|OOC`; default is `IC`.
- Visual Studio project: `TPSolverV2.sln` / `TPSolverV2.vfproj`. CLI and IDE builds must be compared by version output, verification output, result hashes, and executable hashes.

## Version policy (decision)

`solver_v2/VERSION` will become the single human-maintained source of the release version. The build must generate `build_info_m.f90` into the build directory from that value (and must not edit tracked source). Until that generator is implemented, the checked-in `src/core/build_info_m.f90` and `VERSION` are an acknowledged drift and neither is changed by this baseline task.

The generated build information should report solver, theory, architecture, source commit, dirty state, compiler/MKL identity, and build timestamp. A future release task must update the IDE project to use the same generated module or an equivalent pre-build step.

## Source and run identity

Every build manifest records hashes of all source Fortran files, build scripts, project files, `VERSION`, and the executable. It must also record compiler options, MKL linkage, IC/OOC mode, and environment details when available. Every validation run records the executable hash, source/build manifest hash, input and mesh hashes, command line, UTC time, and result artifact hashes. A dirty Git worktree is recorded explicitly; dirty builds are exploratory unless the run is explicitly labeled otherwise.

## Git evidence policy

- Commit source, scripts, schemas, frozen validation inputs, project documentation, small reference data, manifests, acceptance summaries, and small key logs.
- Do not commit rebuildable objects, modules, executables, build directories, IDE caches, or Python caches.
- Large field outputs and full run payloads remain locally preserved or move to an external archive; commit their manifest, hashes, input identity, executable identity, and validation status.
- `final_results/` is historical and **SUPERSEDED** for the graded-hole physics package. Preserve its README and provenance, but do not present it as current validation evidence.

## Baseline verification

Required before the first commit: clean scripted build in a fresh build directory, `version`, and `verify V01`. Existing frozen evidence may support the wider status summary; the full production matrix is a later task. No core `.f90` implementation is changed by this policy.
