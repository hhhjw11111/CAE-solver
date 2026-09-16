# Historical result package — geometry audit superseded prior hole results

Status: **NOT CURRENTLY VERIFIED for the curved-hole physics results**. On 2026-09-06 an independent HEX20 boundary-trace audit found negative Jacobians in the old graded-hole meshes. The old G10–G13 hole results are retained for traceability but are superseded.

- `final_validation_report.md`: historical Chinese report; read together with `../docs/geometry_audit_findings_20260906.md`.
- `figures/`: editable SVG, PDF, 600 dpi TIFF and preview PNG figures.
- `data/`: source CSV files and the R090-L2 deterministic nodal result used for visualization.
- `qa/figure_manifest.json`: figure/source SHA-256 manifest.

Do not run the old `postprocess_final_results.py` to regenerate current paper figures:
its pinned evidence and nodal input still refer to superseded meshes. The entry point is
temporarily guarded until the repaired results and source bindings are audited.

The old quantitative hole metrics must not be used as current physical evidence until G12 and G13 are rerun on the corrected meshes. The top-surface nodal maps are explanatory visualizations and must not replace independent geometry and balance checks.
