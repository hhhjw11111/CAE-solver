"""Create an editable vector schematic for the frozen A1 analytic model."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "figures" / "A1_analytic_geometry_editable_20260907"
STEM = OUT / "A1_analytic_geometry_boundary_conditions"

MM = 1.0 / 25.4
INK = "#26343D"
MUTED = "#62727D"
GRID = "#B7C2C9"
FACE = "#F5F7F8"
TOP = "#E5EAED"
BLUE = "#4D86B5"
BLUE_FILL = "#DDEBF6"
RED = "#C85F55"
RED_FILL = "#F5DDD8"
TEAL = "#3E8D7E"


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial", "DejaVu Sans"],
        "font.size": 7.0,
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def arrow(ax, start, end, color=INK, style="-|>", width=1.0, size=9, zorder=5):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle=style,
        mutation_scale=size,
        linewidth=width,
        color=color,
        shrinkA=0,
        shrinkB=0,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def callout(ax, xy, width, height, title, lines, edge, fill):
    ax.add_patch(
        Rectangle(
            xy,
            width,
            height,
            facecolor=fill,
            edgecolor=edge,
            linewidth=1.0,
            zorder=8,
        )
    )
    ax.text(
        xy[0] + 0.012,
        xy[1] + height - 0.025,
        title,
        ha="left",
        va="top",
        fontsize=7.6,
        fontweight="bold",
        color=INK,
        zorder=9,
    )
    ax.text(
        xy[0] + 0.012,
        xy[1] + height - 0.070,
        "\n".join(lines),
        ha="left",
        va="top",
        fontsize=6.55,
        linespacing=1.35,
        color=INK,
        zorder=9,
    )


def ground(ax, x, y, color):
    ax.plot([x, x], [y, y - 0.020], color=color, lw=0.9, zorder=7)
    for dy, half in [(0.020, 0.012), (0.028, 0.008), (0.035, 0.004)]:
        ax.plot([x - half, x + half], [y - dy, y - dy], color=color, lw=0.9, zorder=7)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(183 * MM, 92 * MM))
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")

    ax.text(
        0.018,
        0.965,
        "A1  一维热–压电–双载流子耦合解析模型",
        ha="left",
        va="top",
        fontsize=10.5,
        fontweight="bold",
        color=INK,
    )
    ax.text(
        0.982,
        0.960,
        "几何：100 nm × 10 nm × 10 nm   |   ΔT = 10 K",
        ha="right",
        va="top",
        fontsize=7.0,
        color=MUTED,
    )

    # Main prism, in the same visual direction as the supplied reference.
    x0, x1 = 0.205, 0.790
    y0, y1 = 0.365, 0.565
    dx, dy = 0.040, 0.070
    p0, p1, p2, p3 = (x0, y0), (x1, y0), (x1, y1), (x0, y1)
    q0, q1, q2, q3 = (
        (x0 + dx, y0 + dy),
        (x1 + dx, y0 + dy),
        (x1 + dx, y1 + dy),
        (x0 + dx, y1 + dy),
    )

    ax.add_patch(Polygon([p3, p2, q2, q3], closed=True, facecolor=TOP, edgecolor=INK, lw=1.0, zorder=2))
    ax.add_patch(Polygon([p0, p1, p2, p3], closed=True, facecolor=FACE, edgecolor=INK, lw=1.1, zorder=3))
    ax.add_patch(Polygon([p0, p3, q3, q0], closed=True, facecolor=BLUE_FILL, edgecolor=BLUE, lw=1.5, zorder=6))
    ax.add_patch(Polygon([p1, q1, q2, p2], closed=True, facecolor=RED_FILL, edgecolor=RED, lw=1.5, zorder=6))
    for a, b in [(q0, q1), (q1, q2), (q2, q3), (q3, q0), (p0, q0), (p1, q1), (p2, q2), (p3, q3)]:
        ax.plot([a[0], b[0]], [a[1], b[1]], color=INK, lw=0.9, zorder=4)

    # Twenty axial HEX20 cells; the light lines keep the model legible.
    for i in range(1, 20):
        t = i / 20.0
        xf = x0 + t * (x1 - x0)
        xb = x0 + dx + t * (x1 - x0)
        ax.plot([xf, xf], [y0, y1], color=GRID, lw=0.38, alpha=0.75, zorder=4)
        ax.plot([xf, xb], [y1, y1 + dy], color=GRID, lw=0.34, alpha=0.65, zorder=4)

    # One-dimensional centre line and domain label.
    ax.plot([x0 + 0.010, x1 + dx - 0.010], [(y0 + y1) / 2 + 0.015] * 2, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=5)
    ax.text(
        (x0 + x1 + dx) / 2,
        0.475,
        r"$\Omega=[0,L]\times[0,H]\times[0,B]$",
        ha="center",
        va="center",
        fontsize=7.5,
        fontweight="bold",
        color=INK,
        zorder=7,
    )
    ax.text(
        (x0 + x1 + dx) / 2,
        0.425,
        r"解析场仅随 $x_1$ 变化",
        ha="center",
        va="center",
        fontsize=6.5,
        color=MUTED,
        zorder=7,
    )

    # Dimensions.
    arrow(ax, (x0, 0.290), (x1, 0.290), style="<->", width=0.9, size=8)
    ax.plot([x0, x0], [0.300, y0 - 0.010], color=INK, lw=0.65)
    ax.plot([x1, x1], [0.300, y0 - 0.010], color=INK, lw=0.65)
    ax.text((x0 + x1) / 2, 0.260, "L = 100 nm", ha="center", va="center", fontsize=7.0, color=INK)
    ax.text((x0 + x1) / 2, 0.225, r"20 × HEX20，$h_x=5$ nm", ha="center", va="center", fontsize=6.2, color=MUTED)

    ax.text(x1 + 0.014, (y0 + y1) / 2, "H = 10 nm", ha="center", va="center", fontsize=5.8, color=RED, rotation=90, zorder=8)
    arrow(ax, (x0 + 0.006, y1 + 0.014), (x0 + dx + 0.006, y1 + dy + 0.014), style="<->", width=0.8, size=7)
    ax.text(x0 + dx + 0.016, y1 + dy + 0.020, "B = 10 nm", ha="left", va="bottom", fontsize=6.2, color=INK)

    # Temperature loading and electrical grounding at both prescribed-potential ends.
    arrow(ax, (x1 + dx + 0.010, y1 + 0.010), (x1 + dx + 0.010, y1 + 0.105), color=RED, width=1.2, size=10)
    ax.text(x1 + dx - 0.006, y1 + 0.105, "ΔT = 10 K", ha="right", va="center", fontsize=6.8, fontweight="bold", color=RED)
    ground(ax, x0 - 0.010, y0 - 0.002, BLUE)
    ground(ax, x1 + dx + 0.010, y0 - 0.002, RED)

    # Coordinate triad.
    origin = (0.835, 0.180)
    arrow(ax, origin, (0.910, 0.180), color=RED, width=1.1, size=9)
    arrow(ax, origin, (0.835, 0.245), color=TEAL, width=1.1, size=9)
    arrow(ax, origin, (0.795, 0.140), color=BLUE, width=1.1, size=9)
    ax.text(0.916, 0.180, r"$x_1$", ha="left", va="center", fontsize=7.0, color=RED)
    ax.text(0.835, 0.251, r"$x_2$", ha="center", va="bottom", fontsize=7.0, color=TEAL)
    ax.text(0.787, 0.132, r"$x_3$", ha="right", va="top", fontsize=7.0, color=BLUE)

    # Boundary-condition callouts.
    callout(
        ax,
        (0.015, 0.345),
        0.168,
        0.300,
        r"$\Gamma_0:\ x_1=0$",
        [r"$u_1=0$", r"$\varphi=0$", r"$\eta_n=\eta_p=0$", r"$\theta=0$ K"],
        BLUE,
        BLUE_FILL,
    )
    arrow(ax, (0.183, 0.485), (x0 - 0.004, 0.485), color=BLUE, width=1.0, size=9, zorder=10)

    callout(
        ax,
        (0.850, 0.365),
        0.142,
        0.275,
        r"$\Gamma_L:\ x_1=L$",
        [r"$\sigma_{11}=0$", r"$\varphi=0$", r"$J_n\!\cdot\!n=J_p\!\cdot\!n=0$", r"$\theta=10$ K"],
        RED,
        RED_FILL,
    )
    arrow(ax, (0.850, 0.490), (x1 + dx + 0.004, 0.490), color=RED, width=1.0, size=9, zorder=10)

    callout(
        ax,
        (0.260, 0.705),
        0.465,
        0.170,
        r"$\Gamma_s$：四个侧面",
        [r"$\mathbf{t}=0$；$\mathbf{D}\!\cdot\!\mathbf{n}=0$；$\mathbf{J}_n\!\cdot\!\mathbf{n}=\mathbf{J}_p\!\cdot\!\mathbf{n}=0$；$\mathbf{q}_T\!\cdot\!\mathbf{n}=0$", r"一维约束：$u_2=u_3=0$，$\partial(\cdot)/\partial x_2=\partial(\cdot)/\partial x_3=0$"],
        GRID,
        "#F0F3F5",
    )
    arrow(ax, (0.493, 0.705), (0.493, y1 + dy + 0.006), color=MUTED, width=0.9, size=8, zorder=10)

    ax.text(
        0.018,
        0.070,
        r"$\theta=T-T_0$，$\eta_n=\Delta n/n_0$，$\eta_p=\Delta p/p_0$；焦耳热 $Q_J=0$，复合项 $R=0$",
        ha="left",
        va="center",
        fontsize=6.4,
        color=MUTED,
    )
    ax.text(
        0.982,
        0.070,
        r"验证变量：$u_1$、$\varphi$、$\eta_n$、$\eta_p$、$\theta$",
        ha="right",
        va="center",
        fontsize=6.4,
        color=MUTED,
    )

    fig.savefig(STEM.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.02)
    fig.savefig(STEM.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.02)
    fig.savefig(STEM.with_suffix(".png"), dpi=450, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(STEM.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

    with (OUT / "A1_analytic_geometry_parameters.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["item", "value", "unit_or_definition"])
        writer.writerows(
            [
                ("length_L", 100, "nm"),
                ("cross_section_side_b", 10, "nm"),
                ("temperature_difference", 10, "K"),
                ("axial_hex20_elements", 20, "elements"),
                ("axial_element_length", 5, "nm"),
                ("transverse_elements", "1 x 1", "HEX20"),
                ("joule_heating", 0, "disabled"),
                ("recombination", 0, "disabled"),
            ]
        )

    contract = {
        "core_conclusion": "The A1 code is tested on a 100 nm by 10 nm by 10 nm prism with one-dimensional fields and explicit thermal, electrical, mechanical, and carrier boundary conditions.",
        "archetype": "schematic-led validation figure",
        "backend": "Python/matplotlib",
        "editable_master": "SVG with live text and separate vector objects",
        "exports": ["SVG", "PDF", "PNG 450 dpi", "TIFF 600 dpi"],
        "geometry": {"L_nm": 100, "b_nm": 10, "axial_hex20_elements": 20},
        "temperature_difference_K": 10,
        "boundary_conditions_verified_against": "the current A1 finite-element implementation",
        "statistics": "not applicable; deterministic analytical verification schematic",
    }
    (OUT / "A1_analytic_geometry_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"A1_EDITABLE_GEOMETRY_OK={STEM}")


if __name__ == "__main__":
    main()
