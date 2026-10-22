# Final figure QA

Date: 2026-09-03

Status: **PASS**

## Reproducibility

- `postprocess_final_results.py` compiles with `py_compile` and runs without warnings.
- Figure manifest hashes match all 12 generated outputs.
- Source hashes are recorded for G05, G12 and the R090-L2 nodal result.
- Source-table row counts are 20 for analytic convergence and 3 for the radius matrix.

## Vector and raster checks

- All three SVG files contain editable `<text>` elements.
- All three PDF files contain one page.
- TIFF outputs are 600 dpi:
  - Figure 1: 4392 × 1723 px;
  - Figure 2: 4391 × 3377 px;
  - Figure 3: 4392 × 2290 px.

## Visual inspection

- Figure 1: convergence curves, reference slopes, residual bars, labels and legends are readable; no clipping or overlap was found.
- Figure 2: the hole is visible in the deformed three-dimensional surface and all four top-surface maps; units, colorbars, the ×1000 displacement statement and the continuation history are explicit; no clipped data were found.
- Figure 3: the four radius-response groups use consistent axes, markers and zero references; labels and legends are readable.
- Color is not the only encoding in the line figures; markers and line styles distinguish series.
- Figure 2 nodal contours are identified as explanatory interpolation in the report. Quantitative claims use the independent G10–G13 metrics.

## Scientific interpretation guardrails

- The figures do not claim experimental calibration.
- The three-dimensional displacement is explicitly exaggerated by 1000.
- Local Joule concentration and mean temperature response are kept as separate quantities.
- The radius matrix is described as a result for the frozen model and boundary conditions, not as a universal physical law.
