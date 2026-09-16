"""Draw the minimal A1 boundary-condition schematic from a beam template."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT.parent / "output" / "figures" / "A1_MINIMAL_BOUNDARY_TEMPLATE_20260908"
STEM = OUT / "A1_minimal_boundary_model"

INK = "#1D2832"
MUTED = "#62727F"
BEAM = "#38A9D6"
BEAM_DARK = "#1682B3"
ELECTRODE = "#E7A428"
COLD = "#2F78B7"
HOT = "#D6453D"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial", "DejaVu Sans"],
        "font.size": 8.0,
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def arrow(ax, start, end, color=INK, style="-|>", lw=1.0, size=9):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle=style,
            mutation_scale=size,
            linewidth=lw,
            color=color,
            shrinkA=0,
            shrinkB=0,
            clip_on=False,
        )
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.20, 3.45))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    x0, x1 = 0.170, 0.830
    y0, y1 = 0.345, 0.555
    electrode_w = 0.012

    # Beam and the two short-circuit electrodes.
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=BEAM, edgecolor=BEAM_DARK, lw=1.05))
    ax.add_patch(Rectangle((x0, y0), electrode_w, y1 - y0, facecolor=ELECTRODE, edgecolor=INK, lw=0.7))
    ax.add_patch(Rectangle((x1 - electrode_w, y0), electrode_w, y1 - y0, facecolor=ELECTRODE, edgecolor=INK, lw=0.7))
    ax.text(
        0.490,
        0.455,
        "热—压电双载流子半导体",
        ha="center",
        va="center",
        fontsize=10.0,
        color=INK,
        weight="bold",
    )
    ax.text(0.490, 0.408, "Thermo-piezoelectric semiconductor (n–p)", ha="center", va="center", fontsize=6.6, color=INK)

    # Axial constraint at the left, retained from the reference visual language.
    wall_x = x0 - 0.018
    ax.plot([wall_x, wall_x], [y0 - 0.080, y1 + 0.095], color=INK, lw=1.3)
    for y in [y0 - 0.055 + i * 0.045 for i in range(8)]:
        ax.plot([wall_x - 0.025, wall_x], [y - 0.028, y], color=INK, lw=0.9)
    ax.plot([wall_x, x0], [(y0 + y1) / 2, (y0 + y1) / 2], color=INK, lw=1.0)

    # Coordinate axes.
    origin = (wall_x - 0.030, y1 + 0.095)
    arrow(ax, origin, (origin[0] + 0.095, origin[1]), color=INK, lw=1.0)
    arrow(ax, origin, (origin[0], origin[1] + 0.110), color=INK, lw=1.0)
    arrow(ax, origin, (origin[0] - 0.040, origin[1] - 0.050), color=MUTED, lw=0.9)
    ax.text(origin[0] + 0.102, origin[1] - 0.006, r"$x_1$", fontsize=9.0, va="center")
    ax.text(origin[0] - 0.005, origin[1] + 0.120, r"$x_2$", fontsize=9.0, ha="center")
    ax.text(origin[0] - 0.052, origin[1] - 0.065, r"$x_3$", fontsize=8.4, color=MUTED)

    # Essential and natural endpoint conditions.
    ax.text(
        x0 + 0.025,
        0.955,
        r"$x_1=0$" "\n" r"$u_1=0,\ \theta=0$" "\n" r"$\eta_n=\eta_p=0$",
        ha="left",
        va="top",
        fontsize=7.7,
        color=COLD,
        weight="bold",
        linespacing=1.20,
    )
    arrow(ax, (x0 + 0.025, 0.790), (x0 + 0.025, y1 + 0.012), color=COLD, lw=1.10)

    ax.text(
        x1 - 0.025,
        0.955,
        r"$x_1=L$" "\n" r"$\sigma_{11}=0,\ \theta=\Delta T=10\ \mathrm{K}$" "\n" r"$J_1^n=J_1^p=0$",
        ha="right",
        va="top",
        fontsize=7.7,
        color=HOT,
        weight="bold",
        linespacing=1.20,
    )
    arrow(ax, (x1 - 0.025, 0.790), (x1 - 0.025, y1 + 0.012), color=HOT, lw=1.10)

    ax.text(
        0.500,
        0.695,
        r"侧面自然边界：$\boldsymbol{t}=0,\ D\!\cdot\boldsymbol{\nu}=J_n\!\cdot\boldsymbol{\nu}=J_p\!\cdot\boldsymbol{\nu}=\boldsymbol{h}\!\cdot\boldsymbol{\nu}=0$",
        ha="center",
        va="center",
        fontsize=7.3,
        color=MUTED,
    )

    # Short-circuit wiring: both endpoint electrodes share phi=0.
    wire_y = 0.155
    ax.plot([x0 + electrode_w / 2, x0 + electrode_w / 2], [y0, wire_y], color=INK, lw=1.0)
    ax.plot([x1 - electrode_w / 2, x1 - electrode_w / 2], [y0, wire_y], color=INK, lw=1.0)
    ax.plot([x0 + electrode_w / 2, x1 - electrode_w / 2], [wire_y, wire_y], color=INK, lw=1.0)
    ground_x = 0.500
    ax.plot([ground_x, ground_x], [wire_y, wire_y - 0.036], color=INK, lw=1.0)
    for width, yy in [(0.052, wire_y - 0.036), (0.034, wire_y - 0.052), (0.016, wire_y - 0.068)]:
        ax.plot([ground_x - width / 2, ground_x + width / 2], [yy, yy], color=INK, lw=1.0)
    ax.text(0.500, 0.040, r"短路电极：$\varphi(0)=\varphi(L)=0$", ha="center", va="center", fontsize=8.0, color=INK)

    # Geometry dimensions.
    dimension_y = 0.275
    arrow(ax, (x0, dimension_y), (x1, dimension_y), style="<->", lw=0.85, size=8)
    ax.text(0.500, dimension_y - 0.020, r"$L=100\ \mathrm{nm}$", ha="center", va="top", fontsize=7.8)
    arrow(ax, (0.855, y0), (0.855, y1), style="<->", lw=0.85, size=8)
    ax.text(0.873, (y0 + y1) / 2, r"$H=10\ \mathrm{nm}$", rotation=90, ha="center", va="center", fontsize=7.6)
    ax.text(0.812, y0 + 0.020, r"$B=10\ \mathrm{nm}$（面外）", ha="right", va="bottom", fontsize=7.2, color=INK)

    fig.savefig(STEM.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.035)
    fig.savefig(STEM.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.035)
    fig.savefig(STEM.with_suffix(".png"), dpi=450, bbox_inches="tight", pad_inches=0.035)
    fig.savefig(
        STEM.with_suffix(".tiff"),
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.035,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)

    contract = {
        "core_conclusion": "A1 is an axially constrained, traction-free one-dimensional thermo-piezoelectric two-carrier bar with short-circuit electrical endpoints and asymmetric carrier boundary conditions.",
        "template_changes": [
            "Remove the transverse point load Q.",
            "Replace flexoelectric semiconductor with thermo-piezoelectric two-carrier semiconductor.",
            "Retain the left support symbol but label it as u1=0 rather than a full 3D clamp.",
            "Add endpoint temperature conditions theta(0)=0 and theta(L)=10 K.",
            "Use electrodes on both endpoints and connect them to the same ground: phi(0)=phi(L)=0.",
            "Add eta_n=eta_p=0 at x1=0 and zero electron/hole normal current at x1=L.",
            "Add natural insulation/no-flux conditions on lateral surfaces.",
        ],
        "geometry_nm": {"L": 100.0, "H": 10.0, "B": 10.0},
        "backend": "Python/matplotlib only",
        "editability": "SVG retains editable text and vector objects; PDF embeds TrueType text.",
    }
    (OUT / "A1_minimal_boundary_model_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"A1_MINIMAL_BOUNDARY_MODEL_OK={STEM}")


if __name__ == "__main__":
    main()
