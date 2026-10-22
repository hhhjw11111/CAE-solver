"""Create an editable A1 geometry schematic and FEM surface-cloud figure.

The plotted fields come only from the archived 20-element A1 nodal output.
The analytical A1 solution is one-dimensional, so each field is averaged over
nodes sharing the same x coordinate before it is mapped uniformly across the
10 nm x 10 nm cross-section.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "runs" / "A1_100NM_010NM_NE20_DT10K_20260907"
OUT = ROOT.parent / "output" / "figures" / "A1_GEOMETRY_AND_CLOUDS_20260908"
STEM = OUT / "A1_geometry_and_FEM_clouds_vertical"

INK = "#202A33"
MUTED = "#65727C"
BLUE = "#5B8DB8"
BLUE_LIGHT = "#DCEAF5"
RED = "#C85D52"
RED_LIGHT = "#F3DDD8"
FACE = "#EEF1F3"
EDGE = "#26343F"

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


def read_axial_fem() -> dict[str, np.ndarray]:
    raw = np.genfromtxt(DATA / "A1_ne20_nodes.csv", delimiter=",", names=True)
    x_unique = np.unique(np.round(raw["x_m"], 18))
    fields = {
        "ux": "u_m",
        "phi": "phi_V",
        "eta_n": "eta_n",
        "eta_p": "eta_p",
    }
    result: dict[str, np.ndarray] = {"x_m": x_unique}
    cross_section_spread = {}
    rounded = np.round(raw["x_m"], 18)
    for name, column in fields.items():
        values = np.asarray(
            [np.mean(raw[column][rounded == x_value]) for x_value in x_unique],
            dtype=float,
        )
        spreads = np.asarray(
            [np.ptp(raw[column][rounded == x_value]) for x_value in x_unique],
            dtype=float,
        )
        result[name] = values
        cross_section_spread[name] = float(np.max(spreads))
    result["cross_section_spread"] = cross_section_spread
    return result


def arrow2d(ax, start, end, color=INK, style="-|>", lw=1.0, size=10):
    patch = FancyArrowPatch(
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
    ax.add_patch(patch)


def draw_geometry(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    left, bottom, length, height = 0.205, 0.265, 0.585, 0.285
    dx, dy = 0.065, 0.100
    p0 = np.array([left, bottom])
    p1 = np.array([left + length, bottom])
    p2 = np.array([left + length, bottom + height])
    p3 = np.array([left, bottom + height])
    q0, q1, q2, q3 = p0 + [dx, dy], p1 + [dx, dy], p2 + [dx, dy], p3 + [dx, dy]

    ax.add_patch(Polygon([p3, p2, q2, q3], closed=True, facecolor="#F7F8F9", edgecolor=EDGE, lw=1.0))
    ax.add_patch(Polygon([p0, p1, p2, p3], closed=True, facecolor=FACE, edgecolor=EDGE, lw=1.15))
    ax.add_patch(Polygon([p0, p3, q3, q0], closed=True, facecolor=BLUE_LIGHT, edgecolor=BLUE, lw=1.25))
    ax.add_patch(Polygon([p1, q1, q2, p2], closed=True, facecolor=RED_LIGHT, edgecolor=RED, lw=1.25))
    for a, b in [(q0, q1), (q1, q2), (q2, q3), (q3, q0), (p0, q0), (p1, q1), (p2, q2), (p3, q3)]:
        ax.plot([a[0], b[0]], [a[1], b[1]], color=EDGE, lw=0.85)

    # Clamp hatch at x=0.
    ax.plot([left - 0.012, left - 0.012], [bottom - 0.025, bottom + height + 0.040], color=INK, lw=1.2)
    for yy in np.linspace(bottom - 0.005, bottom + height + 0.025, 10):
        ax.plot([left - 0.052, left - 0.012], [yy - 0.035, yy], color=INK, lw=0.8)

    # Axes.
    origin = np.array([0.805, 0.170])
    arrow2d(ax, origin, origin + [0.095, 0.0], RED, lw=1.1)
    arrow2d(ax, origin, origin + [0.0, 0.105], "#3F8D78", lw=1.1)
    arrow2d(ax, origin, origin + [-0.052, -0.050], BLUE, lw=1.1)
    ax.text(0.910, 0.166, r"$x$", color=RED, fontsize=9)
    ax.text(0.800, 0.286, r"$y$", color="#3F8D78", fontsize=9)
    ax.text(0.741, 0.105, r"$z$", color=BLUE, fontsize=9)

    # Dimensions.
    arrow2d(ax, (left, 0.120), (left + length, 0.120), style="<->", lw=0.9, size=8)
    ax.text(left + length / 2, 0.082, r"$L=100\ \mathrm{nm}$", ha="center", va="center", fontsize=8.2)
    arrow2d(ax, (0.825, bottom), (0.825, bottom + height), style="<->", lw=0.9, size=8)
    ax.text(0.847, bottom + height / 2, r"$H=10\ \mathrm{nm}$", rotation=90, ha="center", va="center", fontsize=8.0)
    ax.text(left + 0.38 * length, bottom + 0.48 * height, r"$B=10\ \mathrm{nm}$", color=MUTED, fontsize=7.6)

    # Boundary-condition callouts.
    ax.add_patch(Rectangle((0.010, 0.250), 0.165, 0.345, facecolor=BLUE_LIGHT, edgecolor=BLUE, lw=0.9))
    ax.text(0.022, 0.565, r"左端 $x=0$", weight="bold", fontsize=7.8, va="top")
    ax.text(
        0.022,
        0.520,
        "$u_x=0$\n$\\varphi=0$\n$\\eta_n=\\eta_p=0$\n$\\theta=0$",
        fontsize=7.2,
        linespacing=1.28,
        va="top",
    )
    arrow2d(ax, (0.176, 0.410), (left - 0.010, 0.410), BLUE, lw=0.9, size=8)

    ax.add_patch(Rectangle((0.865, 0.300), 0.130, 0.310, facecolor=RED_LIGHT, edgecolor=RED, lw=0.9))
    ax.text(0.878, 0.580, r"右端 $x=L$", weight="bold", fontsize=7.8, va="top")
    ax.text(
        0.878,
        0.535,
        "$\\varphi=0$\n$\\theta=\\Delta T=10$ K\n$\\sigma_{xx}=0$\n$J_n\\!\\cdot\\nu=J_p\\!\\cdot\\nu=0$",
        fontsize=6.8,
        linespacing=1.25,
        va="top",
    )
    arrow2d(ax, (0.864, 0.445), (q1[0] + 0.006, 0.445), RED, lw=0.9, size=8)

    ax.text(
        0.500,
        0.755,
        r"侧面：零牵引、热绝缘，$D\!\cdot\nu=J_n\!\cdot\nu=J_p\!\cdot\nu=0$",
        ha="center",
        va="center",
        color=MUTED,
        fontsize=7.4,
    )
    arrow2d(ax, (0.500, 0.725), (0.500, q3[1] + 0.005), MUTED, lw=0.8, size=7)
    ax.text(0.500, 0.935, "A1 解析验证模型：三维离散，一维多场响应", ha="center", va="center", fontsize=11.0, weight="bold")
    ax.text(0.500, 0.875, r"20 个 HEX20 单元；$\Delta T=10$ K；$Q_J=0$，$R=0$", ha="center", va="center", fontsize=8.0, color=MUTED)


def prism_cloud(ax, x_nodes, values, title, unit, scale, cmap="coolwarm"):
    x_dense = np.linspace(float(x_nodes[0]), float(x_nodes[-1]), 241)
    v_dense = np.interp(x_dense, x_nodes, values) * scale
    xmin, xmax = float(np.min(v_dense)), float(np.max(v_dense))
    if np.isclose(xmin, xmax):
        xmax = xmin + 1.0
    norm = Normalize(vmin=xmin, vmax=xmax)
    colormap = mpl.colormaps[cmap]

    # Reference-style oblique prism. Each colored strip is a vector polygon.
    t_dense = (x_dense - x_dense[0]) / (x_dense[-1] - x_dense[0])
    left, bottom, length, height = 0.070, 0.235, 0.755, 0.285
    dx, dy = 0.045, 0.080
    p0 = np.array([left, bottom])
    p1 = np.array([left + length, bottom])
    p2 = np.array([left + length, bottom + height])
    p3 = np.array([left, bottom + height])
    q0, q1, q2, q3 = p0 + [dx, dy], p1 + [dx, dy], p2 + [dx, dy], p3 + [dx, dy]

    for i in range(len(t_dense) - 1):
        ta, tb = t_dense[i], t_dense[i + 1]
        color = colormap(norm(0.5 * (v_dense[i] + v_dense[i + 1])))
        fa0 = p0 + ta * (p1 - p0)
        fa1 = p0 + tb * (p1 - p0)
        fb1 = p3 + tb * (p2 - p3)
        fb0 = p3 + ta * (p2 - p3)
        ta0 = p3 + ta * (p2 - p3)
        ta1 = p3 + tb * (p2 - p3)
        tb1 = q3 + tb * (q2 - q3)
        tb0 = q3 + ta * (q2 - q3)
        ax.add_patch(Polygon([fa0, fa1, fb1, fb0], closed=True, facecolor=color, edgecolor="none"))
        top_rgb = np.asarray(color[:3]) * 0.82 + 0.18
        ax.add_patch(Polygon([ta0, ta1, tb1, tb0], closed=True, facecolor=top_rgb, edgecolor="none"))

    # End faces and silhouette.
    ax.add_patch(Polygon([p0, p3, q3, q0], closed=True, facecolor=colormap(norm(v_dense[0])), edgecolor=EDGE, lw=0.65))
    ax.add_patch(Polygon([p1, q1, q2, p2], closed=True, facecolor=colormap(norm(v_dense[-1])), edgecolor=EDGE, lw=0.65))
    for a, b in [(p0, p1), (p1, p2), (p2, p3), (p3, p0), (q0, q1), (q1, q2), (q2, q3), (q3, q0), (p0, q0), (p1, q1), (p2, q2), (p3, q3)]:
        ax.plot([a[0], b[0]], [a[1], b[1]], color=EDGE, lw=0.65)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    ax.set_title(title, loc="left", fontsize=9.0, weight="bold", color=INK, pad=0)

    sm = ScalarMappable(norm=norm, cmap=colormap)
    sm.set_array([])
    cbar = ax.figure.colorbar(sm, ax=ax, fraction=0.030, pad=0.018, shrink=0.76, aspect=16)
    cbar.set_label(unit, fontsize=7.3, labelpad=5)
    cbar.ax.tick_params(labelsize=6.5, width=0.6, length=2.5)
    cbar.outline.set_linewidth(0.6)
    return xmin, xmax


def save_figure(fig, stem):
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.025)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.025)
    fig.savefig(stem.with_suffix(".png"), dpi=450, bbox_inches="tight", pad_inches=0.025)
    fig.savefig(
        stem.with_suffix(".tiff"),
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.025,
        pil_kwargs={"compression": "tiff_lzw"},
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = read_axial_fem()
    fields = [
        ("ux", r"轴向位移 $u_x$", "pm", 1.0e12),
        ("phi", r"电势 $\varphi$", "mV", 1.0e3),
        ("eta_n", r"电子相对浓度扰动 $\eta_n$", "%", 100.0),
        ("eta_p", r"空穴相对浓度扰动 $\eta_p$", "%", 100.0),
    ]

    fig = plt.figure(figsize=(7.20, 8.55))
    grid = fig.add_gridspec(
        5,
        1,
        height_ratios=[1.90, 1.0, 1.0, 1.0, 1.0],
        left=0.025,
        right=0.955,
        top=0.985,
        bottom=0.035,
        hspace=0.12,
    )
    geometry_ax = fig.add_subplot(grid[0])
    draw_geometry(geometry_ax)

    ranges = {}
    for row, (key, title, unit, scale) in enumerate(fields, start=1):
        ax = fig.add_subplot(grid[row])
        ranges[key] = prism_cloud(ax, data["x_m"], data[key], title, unit, scale)

    fig.text(
        0.025,
        0.012,
        "云图：A1 的 20 单元 HEX20 数值解（248 节点）。解析场仅随 x 变化，因此每个横截面上的颜色应保持一致。",
        ha="left",
        va="bottom",
        fontsize=6.7,
        color=MUTED,
    )
    save_figure(fig, STEM)
    plt.close(fig)

    # Individual field plates are convenient for PowerPoint placement.
    for key, title, unit, scale in fields:
        fig = plt.figure(figsize=(6.40, 2.05))
        ax = fig.add_subplot(111)
        prism_cloud(ax, data["x_m"], data[key], title, unit, scale)
        fig.subplots_adjust(left=0.02, right=0.92, top=0.88, bottom=0.03)
        save_figure(fig, OUT / f"A1_FEM_cloud_{key}")
        plt.close(fig)

    contract = {
        "core_conclusion": "The A1 FEM fields are uniform on every cross-section and vary only along x, as required by the one-dimensional analytical solution.",
        "archetype": "schematic-led composite",
        "backend": "Python/matplotlib only",
        "source": "Archived A1 Ne=20 nodal output, delta-T=10 K",
        "geometry": {"length_nm": 100.0, "cross_section_nm": [10.0, 10.0], "hex20_elements": 20, "nodes": 248},
        "fields": ["u_x", "phi", "eta_n", "eta_p"],
        "display_ranges": {key: [float(lo), float(hi)] for key, (lo, hi) in ranges.items()},
        "cross_section_spread_raw_units": data["cross_section_spread"],
        "integrity_note": "No synthetic spatial variation was added. Repeated cross-section values are the one-dimensional A1 FEM solution mapped onto the prism surface.",
        "validation_note": "Cloud maps show spatial distribution; analytic correctness remains established by analytical-curve/FEM-node comparisons and independently integrated error norms.",
        "exports": ["editable SVG", "PDF", "PNG 450 dpi", "TIFF 600 dpi"],
    }
    (OUT / "A1_geometry_and_clouds_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"A1_GEOMETRY_AND_CLOUDS_OK={STEM}")


if __name__ == "__main__":
    main()
