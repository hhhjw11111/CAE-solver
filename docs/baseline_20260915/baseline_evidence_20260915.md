# Validated source baseline evidence — 2026-09-15

Purpose: identify why the first Git commit is a reproducible **source baseline**, while separating release identity, validation evidence, generated output, historical runs, and superseded physics. This document does not claim a release baseline.

## Identity and known drift

The canonical physical root is `C:/Users/Administrator/Desktop/gemini修改/Console3/solver_v2` (the Codex workspace path is a junction alias). Git commit is `UNBORN` until the first commit. `VERSION` is `0.14.5-f3-radius-threaded`; the clean executable used for the current smoke checks reports `0.14.9-cauchy-restart`. This known drift is retained verbatim and blocks a release/version baseline, but does not prevent recording today’s source snapshot.

## Current verification evidence

The clean scripted build in `build/BASELINE_AUDIT_20260915` generated a build manifest and compiler options. Its manifest hashes all 54 current Fortran source files; those hashes match the current `src/` tree. The executable is deliberately excluded from Git. A fresh bounded smoke set (`version`, V01, G02, G03, G04, L01, G06, G07, G08) passed; logs and hashes are in this directory. A1/A2 and repaired geometry/mesh evidence are frozen historical evidence and are labelled as such in `verification_evidence.json`; archived manifests differ from current source in listed files.

## Inputs

Validation inputs are retained. F1 fast canonical L0 inputs are retained; F1 main, F2 voltage and F3 radius generators were executed in isolated temporary directories and matched existing generated files byte-for-byte. Generated duplicates are excluded from the first commit while their generator scripts and comparison report are retained. F2 high-voltage restart definitions remain experimental; restart payloads are outside Git and are not claimed runnable from the commit.

## Results and archives

No complete `runs/` directory is staged. Every run payload larger than 10 MB has a SHA-256 and local archive decision in `large_run_archive_manifest.json`; associated small provenance/manifest sidecars are retained only where explicitly selected. Full field samples, nodal fields, binaries and dense dumps remain local; no external archive was created by this audit.

`final_results/README.md` and its figure manifest are retained for traceability. Old graded-hole quantitative outputs are SUPERSEDED by the geometry audit. Analytic assets are marked historical and are not promoted to current coupled-physics evidence.

## Reproduction boundary

Rebuild with `scripts/build_debug.ps1`, then run the commands in the smoke logs. Record the resulting Git commit SHA, clean/dirty state, executable SHA-256, compiler/MKL exact versions, input hashes, and result hashes in the next provenance update.

## Final audit, 2026-09-16

`include_paths.txt` is the sole final set, superseding earlier chat counts. F1 main, F2 standard voltage sweep and F3 generated inputs are excluded (35 exact regeneration matches). No case is declared obsolete without evidence. Experimental restart definitions depend on local payload identities in `cases_classification.json`; they are not smoke prerequisites. Exact compiler/MKL releases remain unrecorded and IDE equivalence was not revalidated. Historical plotting scripts may require parent output archives; source preservation does not assert standalone figure reproduction. README links outside this independent repository are historical references.
