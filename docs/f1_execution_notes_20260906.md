# F1 execution notes

## 2026-09-07 fast-delivery update (supersedes old active-run statements)

The user authorized reduced meshes and workflow optimization for the same-day
six-case presentation. The formal stack-safe L1 result is retained under
`runs/F1_20260906_R0144`. Its following L2 process exited with access violation
(-1073741819); that failure remains preserved. No three-level mesh pass is claimed.

The fast matrix uses L0, a separately frozen threaded binary, in-core PARDISO,
the same material/boundaries and original residual, update, balance and quadrature
tolerances. The first three H0 runs were launched with eight threads; later runs
and the H1-NT repeat use four, MKL_DYNAMIC=FALSE and compact affinity. The earlier
settings are retrospective records, not invented contemporaneous manifests.
The H1-N0 initial/maximum continuation step is 1.0 with original cutback safeguards;
its pre-change input is preserved in evidence/input_history.

The fast build's original manifest says stack-safe although its embedded version
is fast-threaded. Preserve that original and add a clearly marked metadata
correction referencing the same executable hash. Do not relabel the old manifest
as though it had been accurate at build time.

Output reconstruction uses original HEX20 functions, deterministic one-sided
gradients, separate top/bottom evaluations and explicit gaps inside the hole.
The line-point locator passes affine potential/electric-field reproduction tests
on both meshes, including a 59-point hole gap on the 201-point transverse cut.
Final status is recorded in `runs/F1_FAST_HOLE_L0_20260907/fast_audit.json`.

All six L0 cases and the representative rerun have now completed. Six-case
convergence, balance, source reconstruction and quadrature checks pass.
The representative files are NOT byte-identical: maximum potential difference
3.109e-15 V, temperature difference 4.441e-15 K, maximum solver-scaled state
difference 1.862e-13. Keep the strict byte-identity gate failed; do not replace it
with a post-hoc tolerance or claim every original gate passed.
The original PARDISO wrapper sets general ordering/precision options, but no
explicit conditional numerical reproducibility option; this is a plausible
investigation direction, not a proven diagnosis of the observed last-bit differences.

Deliverables are in `output/F1_汇报成果_20260907`: 24 individual line plots,
72 field maps, four current composites, two inherited analytic composites,
all in PNG/SVG/PDF, source CSV/NPZ and a nine-page Chinese report.
Visual review of all current images and report pages completed.
H0's 5x5x3 indicator quadrature under-samples the specified annulus volume by
28.48%; annulus means/P95 are explicitly qualified as coarse discrete trends.
No silent adjustment to solver metrics was made.

## Historical implementation notes

## Directly relevant implementation corrections

1. A1/A2 command-line output directory forwarding was missing. Both initial numerical acceptance logs passed; after forwarding was fixed, the saved-field runs produced byte-identical logs. Original and export logs are retained.
2. N0 explicitly freezes both mobility and diffusion at 300 K. Setting beta to zero alone does not freeze Einstein diffusion. LQ uses linear relative concentrations, reference transport coefficients and a three-stage one-way solution, not one Newton iteration.
3. Candidate ZnO is registered through `run-f1`. The old reference material and the original candidate audit are retained, with calibration gaps unchanged. This is conditional numerical evidence, not experimental validation.
4. The first small-load LQ run used about 4.17 GB peak working memory. At the subsequent resource check, Windows reported only about 4.6 GB of remaining commit capacity. The 8,064-element final hole mesh cannot be assumed to fit in the same in-core configuration. Version 0.14.1 adds an opt-in disk-assisted PARDISO mode, leaving equations, mesh, linear precision, and nonlinear tolerances unchanged. Old callers default to in-core storage.

## Disk-assisted storage, not a new physical experiment

Fortran `iparm(60)=2` selects OOC; `TPS_PARDISO_STORAGE=OOC` is the solver opt-in. Main-case script sets `MKL_PARDISO_OOC_MAX_CORE_SIZE=2048` MB, swap allowance zero, and a task-local scratch-file prefix. One targeted L01 linear-solve check exercises this new storage path. The main matrix and its final repeat must use the same frozen executable and storage configuration.

Authoritative implementation reference: [Intel oneMKL PARDISO Fortran interface, Storage of Matrices](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2026-0/onemkl-pardiso-parallel-direct-sparse-solver-iface.html). This documents disk storage, memory controls, and the fact that a full OOC file prefix overrides the separate directory setting. No iterative approximate linear-solver mode is enabled.

## Small-load interpretation

With c along +z, an x-driven extruded geometry and free traction, the isothermal LQ solution has carrier changes at floating-point noise. In-plane electric fields can be accompanied by a compatible shear displacement while maintaining constant carrier densities. This is a result to explain, not a failure to hide. Consequently, the small-load audit uses relative norms for potential/displacement and fixed n0/p0 scales for physical carrier differences; it never divides by a numerically zero LQ concentration perturbation. It additionally enforces the frozen maximum 1% density-perturbation condition and reports V-squared-normalized thermal differences.

No radius, exponent, conductivity, material, or voltage sweep has been added. Any failed case and its saved fields remain distinguishable from accepted evidence.

## Confirmed zero-source balance issue (0.14.2)

The original small NT run converged, had heat-flux error 5.85e-6 and electron-flux error 5.88e-12, but the maximum path balance was 0.356. This was not accepted blindly. A targeted rerun records all four balance blocks plus the actual free heat residual sum and the original absolute residual sum. At the isothermal path points, the heat residual numerator is about 1e-38 W and its denominator about 1e-38 to 1e-37 W. Their ratio is dominated by roundoff.

The optional zero-source normalizer is the **pre-computable no-hole constant-conductivity power**, q(n0 mu_n0+p0 mu_p0) V^2 Wt/L, equal to 1.97067725982e-6 W for the 0.1 V small case. It is used only when the Joule load factor is exactly zero. Nonzero-source balance calculations and the 1e-8 threshold remain unchanged; legacy callers use the original default because the optional reference is zero. The governing residual/Jacobian and nonlinear solve are unchanged.

The rerun under version 0.14.2 passed, maximum path balance 3.812e-12. It retains the actual heat numerators and denominators in `balance_detail_point` records. The final power, peak temperature, and fields agree with the original IC run to floating-point differences. The old FAIL and old comparison JSON remain archived. A trailing empty field in the original nodal CSV is also fixed for new writes; readers explicitly select the 13 labelled columns when reading older files.

Small-load comparison passed: potential relative norm 3.45208e-5, displacement 1.25061e-3, largest physical carrier shift 2.90308e-4. This is not an additional voltage scan. Main S4-S6 uses the frozen 0.14.2 OOC executable, with the same mode for the representative final repeat.
