# Final verification figure contract

## Core conclusion

The rebuilt solver recovers the expected analytic convergence orders and produces a reproducible, converged three-dimensional nonlinear hole solution with coupled potential, temperature, carrier, and displacement redistribution.

## Figure archetypes

- Analytic verification: quantitative grid.
- Nonlinear hole fields: asymmetric mixed-modality figure with a 3D hero panel.
- Radius dependence: quantitative trend grid.

## Backend and exports

- Backend: Python (`matplotlib`), used exclusively for drawing and QA.
- Final width: 183 mm double-column.
- Primary editable format: SVG.
- Secondary formats: PDF, 600 dpi TIFF, and PNG preview.
- Font: sans serif; SVG text remains editable.

## Panel map

### Analytic verification

- a: V02 potential L2/H1 error versus mesh size, with observed orders.
- b: V06 coupled potential L2/H1 error versus mesh size, with observed orders.
- c: zero-current and thermal analytic residuals on a logarithmic scale.

### Nonlinear hole fields

- a: three-dimensional top-surface displacement geometry, exaggerated and coloured by temperature rise.
- b: top-surface electric potential.
- c: electron concentration redistribution.
- d: hole concentration redistribution.
- e: displacement magnitude.
- f: continuation convergence history.

### Radius dependence

- a: global Joule-power contrast.
- b: local annulus Joule mean and p95 contrasts.
- c: electron and hole redistribution contrasts.
- d: annulus temperature-rise and displacement contrasts.

## Evidence hierarchy

- Hero evidence: theoretical convergence orders and converged R090-L2 coupled fields.
- Validation evidence: zero-current/heat-balance residuals and continuation residual/update histories.
- Robustness: G12 radius matrix and G13 endpoint mesh-sensitivity evidence.

## Source data and integrity

- Every quantitative panel is generated from archived verification logs or a solver-produced nodal result file.
- Parsed source data are exported as CSV next to the figures.
- Spatial panels use only global colour mapping and triangulated interpolation of solver nodes.
- The hole is masked geometrically; no data are invented inside the void.
- The 3D displacement is visually exaggerated and the scale factor is printed in the panel.

## Reviewer risks

- Analytic benchmarks verify discretisation and coupled implementation, not experimental material validity.
- Nodal spatial interpolation is for visualisation; G11/G12 cut-cell integral metrics remain the quantitative evidence.
- Radius trends are conditional on the frozen model, material profile, boundary conditions, and `-5 V` loading.
