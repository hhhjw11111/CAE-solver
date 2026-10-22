"""Targeted checks of line locations and HEX20 affine-field reproduction."""
from pathlib import Path
import numpy as np
import plot_f1_main_results as fe
from plot_fast_f1_gallery import locate_points, sample_line


def main():
    root = Path(__file__).resolve().parents[1]/"cases/F1_FAST_HOLE_L0"
    t = np.linspace(-10, 10, 201)*1e-6
    for geometry in ("H0", "H1"):
        xyz, conn, _ = fe.read_mesh(root/(geometry+"_L0.mesh"))
        state = np.zeros((len(xyz), 7))
        state[:, 3] = -1.5 - 1.5e5*xyz[:, 0]
        for cut, points in [
            ("x", np.column_stack((t, np.full_like(t, 4e-6), np.zeros_like(t)))),
            ("y", np.column_stack((np.zeros_like(t), t, np.full_like(t, .5e-6)))),
        ]:
            located = locate_points(xyz, conn, points, geometry == "H1")
            mask = located[0] >= 0
            data = sample_line(xyz, conn, state, "LQ", located)
            np.testing.assert_allclose(data["phi"][mask], -1.5-1.5e5*points[mask, 0], atol=2e-10)
            np.testing.assert_allclose(data["E"][mask], 1.5e5, rtol=1e-10)
            if geometry == "H1" and cut == "y":
                assert np.count_nonzero(~mask) == 59
                assert np.isnan(data["phi"][~mask]).all()
            else:
                assert mask.all()
            print(f"PASS {geometry}/{cut}: {mask.sum()} material points, map error {located[2]:.3g} m")


if __name__ == "__main__":
    main()
