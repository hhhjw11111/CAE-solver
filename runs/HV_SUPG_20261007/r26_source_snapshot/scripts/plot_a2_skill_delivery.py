"""Build the A2 analytical-solution geometry and FEM-comparison delivery.

Only frozen solver_v2 evidence is read.  The script does not run the solver and
does not modify historical run/output directories.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle


ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "solver_v2" / "runs" / "F1_20260906_R0144" / "A2"
OUT = ROOT / "output" / "figures" / "A2_ANALYTIC_SKILL_20260908"

MM = 1.0 / 25.4
L_UM = 10.0
A_UM = 3.0
T_UM = 1.0
E0 = 1.0e4
REF_BLUE = "#2676D2"

HEX20_R = np.array([
    -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0,
    0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0,
    -1.0, 1.0, 1.0, -1.0,
])
HEX20_S = np.array([
    -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0,
    -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0,
    -1.0, -1.0, 1.0, 1.0,
])
HEX20_T = np.array([
    -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0,
    -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0,
])

INK = "#222A30"
MUTED = "#63717B"
ANALYTIC = "#505D66"
FEM = "#236B8E"
RED = "#C95F52"
RED_LIGHT = "#F4DFDA"
BLUE = "#4F88B8"
BLUE_LIGHT = "#DCEAF5"
ORANGE = "#D8873B"
LIGHT = "#F2F5F6"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial", "DejaVu Sans"],
        "font.size": 8.0,
        "axes.titlesize": 9.0,
        "axes.labelsize": 8.0,
        "legend.fontsize": 7.0,
        "axes.linewidth": 0.75,
        "axes.unicode_minus": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_metrics(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if "=" not in line:
            continue
        key, raw = line.split("=", 1)
        try:
            values[key.strip()] = float(raw.strip())
        except ValueError:
            continue
    return values


def exact_phi(x_m: np.ndarray) -> np.ndarray:
    a_m = A_UM * 1.0e-6
    return -E0 * x_m * (1.0 + (a_m / x_m) ** 2)


def hex20_shape_numpy(r: float, s: float, t: float) -> tuple[np.ndarray, np.ndarray]:
    """Mirror solver_v2/src/mesh/hex20_shape_m.f90 for read-only post-processing."""
    shape = np.zeros(20)
    derivative = np.zeros((3, 20))
    for node in range(8):
        a, b, c = HEX20_R[node], HEX20_S[node], HEX20_T[node]
        ar, bs, ct = 1.0 + a * r, 1.0 + b * s, 1.0 + c * t
        linear = a * r + b * s + c * t - 2.0
        shape[node] = 0.125 * ar * bs * ct * linear
        derivative[0, node] = 0.125 * a * bs * ct * (linear + ar)
        derivative[1, node] = 0.125 * b * ar * ct * (linear + bs)
        derivative[2, node] = 0.125 * c * ar * bs * (linear + ct)
    for node in range(8, 20):
        a, b, c = HEX20_R[node], HEX20_S[node], HEX20_T[node]
        if a == 0.0:
            shape[node] = 0.25 * (1.0 - r * r) * (1.0 + b * s) * (1.0 + c * t)
            derivative[0, node] = -0.5 * r * (1.0 + b * s) * (1.0 + c * t)
            derivative[1, node] = 0.25 * (1.0 - r * r) * b * (1.0 + c * t)
            derivative[2, node] = 0.25 * (1.0 - r * r) * c * (1.0 + b * s)
        elif b == 0.0:
            shape[node] = 0.25 * (1.0 - s * s) * (1.0 + a * r) * (1.0 + c * t)
            derivative[0, node] = 0.25 * (1.0 - s * s) * a * (1.0 + c * t)
            derivative[1, node] = -0.5 * s * (1.0 + a * r) * (1.0 + c * t)
            derivative[2, node] = 0.25 * (1.0 - s * s) * c * (1.0 + a * r)
        else:
            shape[node] = 0.25 * (1.0 - t * t) * (1.0 + a * r) * (1.0 + b * s)
            derivative[0, node] = 0.25 * (1.0 - t * t) * a * (1.0 + b * s)
            derivative[1, node] = 0.25 * (1.0 - t * t) * b * (1.0 + a * r)
            derivative[2, node] = -0.5 * t * (1.0 + a * r) * (1.0 + b * s)
    return shape, derivative


def reconstruct_l2_connectivity(nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Rebuild the frozen nt=32, nr=4 A2 connectivity and audit node ordering."""
    nt, nr = 32, 4
    ids = np.full((2 * nr + 1, 2 * nt, 3), -1, dtype=int)
    coordinates: list[list[float]] = []
    connectivity: list[list[int]] = []
    for et in range(nt):
        for er in range(nr):
            element_nodes: list[int] = []
            for local in range(20):
                ir = 2 * er + int(round(HEX20_R[local])) + 1
                it = (2 * et + int(round(HEX20_S[local])) + 1) % (2 * nt)
                iz = int(round(HEX20_T[local])) + 1
                if ids[ir, it, iz] < 0:
                    angle = np.pi * it / nt
                    outer = L_UM * 1.0e-6 / max(abs(np.cos(angle)), abs(np.sin(angle)))
                    rho = (ir / (2 * nr)) ** 2
                    if ir % 2 == 1:
                        rho = 0.5 * (((ir - 1) / (2 * nr)) ** 2 + ((ir + 1) / (2 * nr)) ** 2)
                    radius = A_UM * 1.0e-6 + rho * (outer - A_UM * 1.0e-6)
                    ids[ir, it, iz] = len(coordinates)
                    coordinates.append([
                        radius * np.cos(angle),
                        radius * np.sin(angle),
                        0.5 * T_UM * 1.0e-6 * (iz - 1),
                    ])
                element_nodes.append(ids[ir, it, iz])
            connectivity.append(element_nodes)
    rebuilt = np.asarray(coordinates)
    saved = np.column_stack((nodes["x_m"], nodes["y_m"], nodes["z_m"]))
    if rebuilt.shape != saved.shape:
        raise RuntimeError(f"A2 L2 reconstruction node-count mismatch: {rebuilt.shape} != {saved.shape}")
    maximum_difference = float(np.max(np.abs(rebuilt - saved)))
    if maximum_difference > 1.0e-15:
        raise RuntimeError(f"A2 L2 reconstruction coordinate mismatch: {maximum_difference:.6e} m")
    return rebuilt, np.asarray(connectivity, dtype=int), maximum_difference


def centerline_fem_joule_samples(nodes: np.ndarray) -> tuple[np.ndarray, float]:
    """Evaluate |grad(phi_h)|^2/E0^2 inside L2 HEX20 elements on y=z=0."""
    coordinates, connectivity, coordinate_difference = reconstruct_l2_connectivity(nodes)
    phi = nodes["phi_V"]
    records: list[tuple[float, float, float, int, float, float, float]] = []
    radial_points = (-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0))
    for side, et in ((-1.0, 16), (1.0, 0)):
        for er in range(4):
            element = et * 4 + er
            conn = connectivity[element]
            xyz_nodes = coordinates[conn]
            for r_natural in radial_points:
                shape, derivative_natural = hex20_shape_numpy(r_natural, -1.0, 0.0)
                jacobian = derivative_natural @ xyz_nodes
                derivative_global = np.linalg.inv(jacobian) @ derivative_natural
                xyz = shape @ xyz_nodes
                electric_field = -(derivative_global @ phi[conn])
                q_fem = float(np.dot(electric_field, electric_field) / E0**2)
                x = float(xyz[0])
                q_exact = float((1.0 - (A_UM * 1.0e-6 / x) ** 2) ** 2)
                x_over_l = (x + L_UM * 1.0e-6) / (2.0 * L_UM * 1.0e-6)
                records.append((x_over_l, x, q_fem, element + 1, r_natural, q_exact, side))
    dtype = [
        ("x1_over_L", "f8"), ("x_m", "f8"), ("Q_fem_over_Qinf", "f8"),
        ("element_id", "i4"), ("r_natural", "f8"), ("Q_exact_over_Qinf", "f8"),
        ("side", "f8"),
    ]
    samples = np.array(records, dtype=dtype)
    samples.sort(order="x1_over_L")
    return samples, coordinate_difference


def save_all(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    base = OUT / stem
    kwargs = {"bbox_inches": "tight", "pad_inches": 0.035, "facecolor": "white"}
    fig.savefig(base.with_suffix(".svg"), **kwargs)
    fig.savefig(base.with_suffix(".pdf"), **kwargs)
    fig.savefig(base.with_suffix(".png"), dpi=450, **kwargs)
    fig.savefig(
        base.with_suffix(".tiff"),
        dpi=600,
        pil_kwargs={"compression": "tiff_lzw"},
        **kwargs,
    )
    plt.close(fig)


def arrow(ax, start, end, color=INK, lw=1.0, scale=10.0):
    item = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=scale,
        color=color,
        linewidth=lw,
        shrinkA=0,
        shrinkB=0,
        clip_on=False,
    )
    ax.add_patch(item)
    return item


def double_arrow(ax, start, end, text, offset=(0.0, 0.0)):
    item = FancyArrowPatch(
        start,
        end,
        arrowstyle="<->",
        mutation_scale=8,
        color=INK,
        linewidth=0.75,
        shrinkA=0,
        shrinkB=0,
        clip_on=False,
    )
    ax.add_patch(item)
    ax.text(
        0.5 * (start[0] + end[0]) + offset[0],
        0.5 * (start[1] + end[1]) + offset[1],
        text,
        ha="center",
        va="center",
        color=INK,
        fontsize=7.0,
    )


def info_box(ax, xy, size, title, lines, face, edge):
    x, y = xy
    w, h = size
    ax.add_patch(Rectangle((x, y), w, h, facecolor=face, edgecolor=edge, lw=0.9))
    ax.text(x + 0.025, y + h - 0.035, title, ha="left", va="top", weight="bold", color=INK)
    ax.text(
        x + 0.025,
        y + h - 0.095,
        "\n".join(lines),
        ha="left",
        va="top",
        fontsize=7.0,
        linespacing=1.42,
        color=INK,
    )


def plot_geometry() -> None:
    fig = plt.figure(figsize=(183 * MM, 98 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.22, 0.95], wspace=0.08,
                          left=0.035, right=0.985, bottom=0.065, top=0.92)
    ax = fig.add_subplot(gs[0, 0])
    notes = fig.add_subplot(gs[0, 1])

    ax.set_aspect("equal")
    ax.set_xlim(-13.4, 13.3)
    ax.set_ylim(-13.1, 14.2)
    ax.axis("off")
    ax.set_title("A2 三维挤出绝缘圆孔：俯视几何与边界位置", pad=8, weight="bold")

    ax.add_patch(Rectangle((-10, -10), 20, 20, facecolor="#FBFCFC", edgecolor=RED, lw=2.0))
    ax.add_patch(Circle((0, 0), A_UM, facecolor="white", edgecolor=BLUE, lw=2.0))

    arrow(ax, (-7.0, 12.1), (7.0, 12.1), color=ORANGE, lw=1.8, scale=12)
    ax.text(0, 12.65, r"远场方向  $\mathbf{E}_0=E_0\mathbf{e}_x$", ha="center", va="bottom",
            color=ORANGE, weight="bold")

    # Hot spots and stagnation points are analytical expectations, not fitted FEM data.
    for y in (-A_UM, A_UM):
        ax.add_patch(Circle((0, y), 0.34, facecolor=RED, edgecolor="white", lw=0.7, zorder=5))
    for x in (-A_UM, A_UM):
        ax.add_patch(Circle((x, 0), 0.26, facecolor="white", edgecolor=BLUE, lw=1.4, zorder=5))
    ax.text(0.65, A_UM, r"热点：$|E_t|/E_0=2$，$Q/Q_\infty=4$", va="center", color=RED)
    ax.text(0.65, -A_UM, r"热点：$\psi=270^\circ$", va="center", color=RED)
    ax.text(A_UM + 0.55, 0.55, "停滞点", color=BLUE, ha="left")
    ax.text(-A_UM - 0.55, 0.55, "停滞点", color=BLUE, ha="right")

    double_arrow(ax, (-10, -11.35), (10, -11.35), r"$20\ \mu\mathrm{m}$", offset=(0, -0.55))
    double_arrow(ax, (-11.35, -10), (-11.35, 10), r"$20\ \mu\mathrm{m}$", offset=(-0.80, 0))
    arrow(ax, (0, 0), (-A_UM * 0.71, A_UM * 0.71), color=INK, lw=0.8, scale=7)
    ax.text(-2.55, 1.42, r"$a=3\ \mu\mathrm{m}$", rotation=-43, color=INK)

    arrow(ax, (-9.0, -8.6), (-6.7, -8.6), color="#C84F45", lw=1.0, scale=9)
    arrow(ax, (-9.0, -8.6), (-9.0, -6.3), color="#3D8A77", lw=1.0, scale=9)
    ax.text(-6.35, -8.6, "$x$", va="center", color="#C84F45")
    ax.text(-9.0, -5.95, "$y$", ha="center", color="#3D8A77")
    ax.text(-9.65, 9.25, r"外周 $\Gamma_{\rm out}$", color=RED, weight="bold")
    ax.text(0, -0.15, "贯穿圆柱孔", ha="center", va="center", color=BLUE, weight="bold")

    notes.set_xlim(0, 1)
    notes.set_ylim(0, 1)
    notes.axis("off")
    notes.set_title("边界条件与三维挤出说明", pad=8, weight="bold")
    info_box(
        notes,
        (0.02, 0.69),
        (0.96, 0.25),
        "正方形四个外侧面",
        [
            r"$\phi=\phi_{\rm ex}(x,y)$，不是仅左右端线性电压",
            r"$\eta_n=\eta_p=0$",
            r"$\phi_{\rm ex}=-E_0x[1+a^2/(x^2+y^2)]$",
        ],
        RED_LIGHT,
        RED,
    )
    info_box(
        notes,
        (0.02, 0.39),
        (0.96, 0.25),
        "圆柱孔壁  $r=a$",
        [
            r"$\mathbf{D}\!\cdot\!\boldsymbol{\nu}=0$",
            r"$\mathbf{J}_n\!\cdot\!\boldsymbol{\nu}=\mathbf{J}_p\!\cdot\!\boldsymbol{\nu}=0$",
            "理想绝缘孔；边界通量在弱式中自然为零",
        ],
        BLUE_LIGHT,
        BLUE,
    )
    info_box(
        notes,
        (0.02, 0.09),
        (0.96, 0.25),
        "厚度方向与固定场",
        [
            r"$z\in[-0.5,0.5]\ \mu\mathrm{m}$，厚度 $t=1\ \mu\mathrm{m}$",
            r"上下表面：$D_z=J_{n,z}=J_{p,z}=0$",
            r"全部节点：$\mathbf{u}=0$，$\theta=T-T_0=0$",
        ],
        LIGHT,
        "#9AA7AE",
    )
    notes.text(
        0.02,
        0.01,
        "合成、常系数、各向同性、等温验证；不代表 ZnO 实验材料。",
        ha="left",
        va="bottom",
        fontsize=6.8,
        color=MUTED,
    )

    save_all(fig, "A2_geometry_boundary_conditions")


def plot_reference_sampling_geometry() -> None:
    """Draw the A2 model in the compact style of the user-provided reference."""
    nodes = np.genfromtxt(RUN / "A2_L2_nodes.csv", delimiter=",", names=True)
    wall = np.genfromtxt(RUN / "A2_L2_wall.csv", delimiter=",", names=True)

    center = nodes[(np.abs(nodes["y_m"]) < 1.0e-15) & (np.abs(nodes["z_m"]) < 1.0e-15)]
    center.sort(order="x_m")
    angle = np.degrees(wall["angle_rad"]) % 360.0
    order = np.argsort(angle)
    wall = wall[order]
    angle = angle[order]
    shown = np.arange(0, len(wall), 10, dtype=int)

    fig, ax = plt.subplots(figsize=(145 * MM, 128 * MM))
    fig.subplots_adjust(left=0.08, right=0.79, bottom=0.16, top=0.91)
    ax.set_aspect("equal")
    ax.set_xlim(-13.7, 16.5)
    ax.set_ylim(-14.0, 13.1)
    ax.axis("off")
    ax.set_title("A2 绝缘圆孔解析验证：几何、边界与折线取样点", fontsize=10.5, weight="bold", pad=8)

    # Reference-style grey domain with direct boundary annotations.
    ax.add_patch(Rectangle((-10, -10), 20, 20, facecolor="#E3E3E3",
                           edgecolor="#525BFF", lw=1.55))
    ax.add_patch(Circle((0, 0), A_UM, facecolor="white", edgecolor=INK, lw=1.15, zorder=4))

    # Coordinate convention and applied far-field direction.
    arrow(ax, (-8.4, 8.2), (-5.7, 8.2), color=INK, lw=0.8, scale=8)
    arrow(ax, (-8.4, 8.2), (-8.4, 10.2), color=INK, lw=0.8, scale=8)
    ax.text(-5.35, 8.2, "$x$", va="center", fontsize=7.0)
    ax.text(-8.4, 10.55, "$y$", ha="center", fontsize=7.0)
    ax.text(-5.15, 8.85, "$x,y$ 为面内；$z$ 为面外", fontsize=6.5, color=MUTED)
    arrow(ax, (-5.0, 6.4), (5.2, 6.4), color=ORANGE, lw=1.25, scale=10)
    ax.text(0.1, 6.85, r"$\mathbf{E}_0=E_0\mathbf{e}_x$", ha="center", color=ORANGE, fontsize=7.5)

    # Exact FEM points used in the visible line plots.
    ax.plot([-10, -A_UM], [0, 0], color=FEM, lw=0.75, ls="--", alpha=0.75, zorder=2)
    ax.plot([A_UM, 10], [0, 0], color=FEM, lw=0.75, ls="--", alpha=0.75, zorder=2)
    ax.plot(center["x_m"] * 1.0e6, np.zeros(len(center)), linestyle="none", marker="o",
            ms=4.1, mfc="white", mec=FEM, mew=1.25, zorder=9)
    ax.plot(wall["x_m"][shown] * 1.0e6, wall["y_m"][shown] * 1.0e6,
            linestyle="none", marker="o", ms=2.15, mfc=RED, mec="white", mew=0.25, zorder=8)

    # Radius and wall-angle convention.
    arrow(ax, (0, 0), (2.12, 2.12), color=INK, lw=0.8, scale=7)
    ax.text(0.85, 1.65, r"$a=3\ \mu\mathrm{m}$", rotation=42, fontsize=7.0)
    arc_t = np.linspace(0, np.pi / 2, 60)
    ax.plot(1.18 * np.cos(arc_t), 1.18 * np.sin(arc_t), color=MUTED, lw=0.7)
    ax.text(1.25, 0.55, r"$\psi$", fontsize=7.0, color=MUTED)
    ax.text(0, -0.35, "贯穿圆柱孔", ha="center", va="center", fontsize=7.2, weight="bold")

    # Dimension lines matching the reference-figure visual grammar.
    double_arrow(ax, (-10, -11.55), (10, -11.55), r"$L_x=20\ \mu\mathrm{m}$", offset=(0, -0.60))
    double_arrow(ax, (14.25, -10), (14.25, 10), r"$L_y=20\ \mu\mathrm{m}$", offset=(1.05, 0))
    ax.text(-9.8, -9.2, r"厚度 $t=1\ \mu\mathrm{m}$，$z\in[-0.5,0.5]\ \mu\mathrm{m}$",
            ha="left", va="center", fontsize=6.7, color=INK)

    # Boundary labels: same placement idea as the reference, but current A2 facts.
    ax.annotate(
        "理想绝缘孔壁（自然边界）\n"
        r"$\mathbf{D}\!\cdot\!\boldsymbol{\nu}=0$，"
        r"$\mathbf{J}_n\!\cdot\!\boldsymbol{\nu}=\mathbf{J}_p\!\cdot\!\boldsymbol{\nu}=0$",
        xy=(-1.75, 2.45), xytext=(-8.8, 4.2), fontsize=6.8, ha="left", va="bottom",
        arrowprops={"arrowstyle": "->", "lw": 0.75, "color": INK}, color=INK,
    )
    ax.text(
        10.45,
        2.8,
        "四个外侧面均施加\n"
        r"$\phi=\phi_{\rm ex}(x,y)$" "\n"
        r"$\eta_n=\eta_p=0$",
        ha="left", va="center", fontsize=7.0, color="#525BFF", linespacing=1.35,
    )
    ax.text(-9.75, -7.8, r"全域验证约束：$\mathbf{u}=0$，$\theta=0$",
            ha="left", va="center", fontsize=6.9, color=INK)

    # Sampling legend below the geometry.
    handles = [
        Line2D([], [], linestyle="none", marker="o", ms=5.0, mfc="white", mec=FEM, mew=1.25,
               label=r"中心线电势 $\phi(x,0,0)$：10 个 L2 节点"),
        Line2D([], [], linestyle="none", marker="o", ms=3.5, mfc=RED, mec=RED,
               label=r"孔壁 $E_t/Q_J$ 角向曲线：68 个显示点"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.42, -0.145), ncol=1,
              handletextpad=0.7, labelspacing=0.55, frameon=False, fontsize=6.8)
    ax.text(
        0.42,
        -0.205,
        r"误差范数由独立 $5\times5\times5$ Gauss 积分计算，不由上述显示点代替。",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=6.4,
        color=MUTED,
        clip_on=False,
    )

    save_all(fig, "A2_reference_style_geometry_sampling")

    with (OUT / "A2_plot_sampling_points.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["set", "point_id", "x_m", "y_m", "z_m", "angle_deg", "quantity", "source"])
        for i, row in enumerate(center, start=1):
            writer.writerow(["centerline", i, row["x_m"], row["y_m"], row["z_m"], "",
                             "phi(x,0,0)", "A2_L2_nodes.csv: y=z=0"])
        for j, idx in enumerate(shown, start=1):
            writer.writerow(["hole_wall_display", j, wall["x_m"][idx], wall["y_m"][idx], 0.0,
                             angle[idx], "Et/E0 and Q/Qinf", "A2_L2_wall.csv: sorted angle, stride 10"])

    sampling_contract = {
        "reference_image_role": "layout inspiration only; no numerical labels or boundary conditions copied",
        "current_model": "A2_F1_insulated_hole in solver_v2",
        "centerline_sampling": {
            "selection": "abs(y)<1e-15 and abs(z)<1e-15 in A2_L2_nodes.csv",
            "count": int(len(center)),
            "quantity": "potential phi on y=z=0, excluding the hole",
        },
        "wall_display_sampling": {
            "selection": "A2_L2_wall.csv sorted by angle modulo 360 degrees, every tenth record",
            "count": int(len(shown)),
            "full_wall_records": int(len(wall)),
            "quantities": ["abs(Et)/E0", "Q_J/Q_inf"],
        },
        "error_sampling": "independent 5x5x5 Gauss integration; not represented by the plotted-point markers",
        "scope_note": "The outer boundary is exact spatial Dirichlet on all four sides; it is not the 5 V one-sided boundary in the reference image.",
    }
    (OUT / "A2_reference_geometry_sampling_contract.json").write_text(
        json.dumps(sampling_contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def plot_results() -> dict[str, float]:
    wall = np.genfromtxt(RUN / "A2_L2_wall.csv", delimiter=",", names=True)
    nodes = np.genfromtxt(RUN / "A2_L2_nodes.csv", delimiter=",", names=True)
    values = read_metrics(RUN / "A2_export_stdout.txt")

    angle = np.degrees(wall["angle_rad"]) % 360.0
    order = np.argsort(angle)
    angle = angle[order]
    wall = wall[order]
    exact_angle = np.linspace(0.0, 360.0, 721)
    exact_et = 2.0 * np.abs(np.sin(np.radians(exact_angle)))
    exact_q = 4.0 * np.sin(np.radians(exact_angle)) ** 2

    center = nodes[(np.abs(nodes["y_m"]) < 1.0e-15) & (np.abs(nodes["z_m"]) < 1.0e-15)]
    center.sort(order="x_m")
    x_left = np.linspace(-L_UM, -A_UM, 201) * 1.0e-6
    x_right = np.linspace(A_UM, L_UM, 201) * 1.0e-6

    fig, axes = plt.subplots(2, 2, figsize=(183 * MM, 128 * MM))
    fig.subplots_adjust(left=0.085, right=0.975, bottom=0.100, top=0.855, wspace=0.30, hspace=0.48)

    ax = axes[0, 0]
    for x, label in ((x_left, "解析解"), (x_right, None)):
        ax.plot(x * 1.0e6, exact_phi(x), color=ANALYTIC, lw=1.3, marker="o",
                markevery=20, ms=2.4, mfc=ANALYTIC, mec=ANALYTIC, label=label)
    ax.plot(center["x_m"] * 1.0e6, center["phi_V"], linestyle="none", marker="o",
            ms=4.2, mfc="white", mec=FEM, mew=1.2, label="FEM（HEX20，L2）")
    ax.axvspan(-A_UM, A_UM, color="#E7ECEF", zorder=-2)
    ax.text(0, 0, "圆孔", ha="center", va="center", color=MUTED)
    ax.set(xlabel=r"中心线位置 $x$（$\mu$m），$y=z=0$", ylabel=r"电势 $\phi$（V）",
           title="中心线电势：数值解与解析解")
    ax.legend(loc="upper right")

    ax = axes[0, 1]
    ax.plot(exact_angle, exact_et, color=ANALYTIC, lw=1.3, marker="o", markevery=36,
            ms=2.4, mfc=ANALYTIC, mec=ANALYTIC, label="解析解")
    ax.plot(angle[::10], wall["Et_over_E0"][::10], linestyle="none", marker="o",
            ms=3.7, mfc="white", mec=FEM, mew=1.15, label="FEM（L2）")
    ax.set(xlabel=r"孔壁角度 $\psi$（°）", ylabel=r"切向电场 $|E_t|/E_0$",
           title="孔壁切向电场", xticks=[0, 90, 180, 270, 360])
    ax.legend(loc="center", bbox_to_anchor=(0.50, 0.70), ncol=2)

    ax = axes[1, 0]
    ax.plot(exact_angle, exact_q, color=ANALYTIC, lw=1.3, marker="o", markevery=36,
            ms=2.4, mfc=ANALYTIC, mec=ANALYTIC, label="解析解")
    ax.plot(angle[::10], wall["Q_over_Qinf"][::10], linestyle="none", marker="o",
            ms=3.7, mfc="white", mec=FEM, mew=1.15, label="FEM（L2）")
    ax.set(xlabel=r"孔壁角度 $\psi$（°）", ylabel=r"归一化焦耳热源 $Q_J/Q_\infty$",
           title="孔壁电场平方对应的焦耳热源", xticks=[0, 90, 180, 270, 360])
    ax.legend(loc="center", bbox_to_anchor=(0.50, 0.62), ncol=2)

    ax = axes[1, 1]
    levels = np.array([0, 1, 2])
    series = [
        ("phi_relative_L2", r"电势相对 $L^2$ 误差", FEM, "o", r"电势 $L^2$"),
        ("E_relative_L2", r"电场相对 $L^2$ 误差", ORANGE, "s", r"电场 $L^2$"),
        ("wall_normal_E_RMS_over_E0", r"孔壁法向场 RMS/$E_0$", ANALYTIC, "^", r"法向场 RMS/$E_0$"),
    ]
    for key, label, color, marker, short in series:
        y = np.array([values[f"a2_L{i}_{key}"] for i in levels])
        ax.semilogy(levels, y, color=color, marker=marker, ms=4.2, lw=1.15, label=label)
        ax.text(2.08, y[-1], short, color=color, fontsize=6.7, va="center", ha="left")
    ax.set(xlabel="网格层级", ylabel="归一化误差", title="独立积分误差随网格加密下降",
           xticks=[0, 1, 2], xticklabels=["L0", "L1", "L2"], xlim=(-0.10, 2.72))
    ax.text(
        0.03,
        0.05,
        "$L2$：$|E_t|_{\\max}/E_0=1.999857$\n"
        "$L2$：$Q_{\\max}/Q_\\infty=3.999428$\n"
        "判定：PASS",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.0,
        color=INK,
        bbox={"boxstyle": "round,pad=0.28", "facecolor": LIGHT, "edgecolor": "#AAB5BB", "lw": 0.7},
    )

    for ax in axes.ravel():
        ax.tick_params(direction="in", length=3.2, width=0.7)
        ax.grid(False)

    fig.suptitle("A2：三维挤出绝缘圆孔的恒电导电场集中解析验证",
                 fontsize=11, weight="bold", y=0.975)
    save_all(fig, "A2_analytic_fem_results")

    # Auditable source tables used by the figures.
    with (OUT / "A2_centerline_comparison.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x_m", "phi_fem_V", "phi_exact_at_fem_node_V"])
        for row in center:
            writer.writerow([row["x_m"], row["phi_V"], row["phi_exact_V"]])

    with (OUT / "A2_wall_comparison.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["angle_deg", "En_over_E0", "Et_fem_over_E0", "Et_exact_over_E0",
                         "Q_fem_over_Qinf", "Q_exact_over_Qinf"])
        for i in range(len(wall)):
            writer.writerow([angle[i], wall["En_over_E0"][i], wall["Et_over_E0"][i],
                             wall["Et_exact_over_E0"][i], wall["Q_over_Qinf"][i],
                             wall["Q_exact_over_Qinf"][i]])

    keys = ["phi_relative_L2", "E_relative_L2", "wall_normal_E_RMS_over_E0",
            "wall_Et_peak_over_E0", "wall_Q_peak_over_Qinf"]
    with (OUT / "A2_error_convergence.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["level", *keys, "nodes", "newton_iterations"])
        for i in levels:
            writer.writerow([f"L{i}", *[values[f"a2_L{i}_{key}"] for key in keys],
                             int(values[f"a2_L{i}_nodes"]), int(values[f"a2_L{i}_newton_iterations"])])
    return values


def plot_four_panel_composite() -> None:
    """Publication-focused four-panel A2 composite with an En boundary inset."""
    wall = np.genfromtxt(RUN / "A2_L2_wall.csv", delimiter=",", names=True)
    nodes = np.genfromtxt(RUN / "A2_L2_nodes.csv", delimiter=",", names=True)
    values = read_metrics(RUN / "A2_export_stdout.txt")

    angle = np.degrees(wall["angle_rad"]) % 360.0
    order = np.argsort(angle)
    angle = angle[order]
    wall = wall[order]
    shown = np.arange(0, len(wall), 10, dtype=int)
    center = nodes[(np.abs(nodes["y_m"]) < 1.0e-15) & (np.abs(nodes["z_m"]) < 1.0e-15)]
    center.sort(order="x_m")
    exact_angle = np.linspace(0.0, 360.0, 721)

    fig, axes = plt.subplots(2, 2, figsize=(183 * MM, 130 * MM))
    fig.subplots_adjust(left=0.065, right=0.98, bottom=0.125, top=0.89, wspace=0.25, hspace=0.36)
    geo, phi_ax, field_ax, q_ax = axes.ravel()

    # Geometry and actual plotted sampling sets.
    geo.set_aspect("equal")
    geo.set_xlim(-12.8, 13.5)
    geo.set_ylim(-12.5, 12.3)
    geo.axis("off")
    geo.set_title("几何、边界与折线取样", pad=5)
    geo.add_patch(Rectangle((-10, -10), 20, 20, facecolor="#E5E5E5",
                            edgecolor="#525BFF", lw=1.4))
    geo.add_patch(Circle((0, 0), A_UM, facecolor="white", edgecolor=INK, lw=1.0, zorder=4))
    geo.plot([-10, -A_UM], [0, 0], color=FEM, lw=0.65, ls="--", zorder=2)
    geo.plot([A_UM, 10], [0, 0], color=FEM, lw=0.65, ls="--", zorder=2)
    geo.plot(wall["x_m"][shown] * 1e6, wall["y_m"][shown] * 1e6, linestyle="none",
             marker="o", ms=1.9, mfc=RED, mec="white", mew=0.2, zorder=7)
    geo.plot(center["x_m"] * 1e6, np.zeros(len(center)), linestyle="none", marker="o",
             ms=3.7, mfc="white", mec=FEM, mew=1.1, zorder=8)
    arrow(geo, (-5.0, 7.2), (5.0, 7.2), color=ORANGE, lw=1.2, scale=9)
    geo.text(0, 7.7, r"$\mathbf{E}_0=E_0\mathbf{e}_x$", ha="center", color=ORANGE, fontsize=7.0)
    arrow(geo, (0, 0), (2.12, 2.12), color=INK, lw=0.7, scale=6)
    geo.text(1.05, 1.55, r"$a=3\ \mu\mathrm{m}$", rotation=41, fontsize=6.4)
    geo.text(-9.6, 9.1, r"$L_x=L_y=20\ \mu\mathrm{m}$；$t=1\ \mu\mathrm{m}$",
             fontsize=6.4, ha="left")
    geo.text(-9.6, 5.3, r"外周：$\phi=\phi_{\rm ex}$，$\eta_n=\eta_p=0$", fontsize=6.25)
    geo.annotate(r"孔壁：$D_\nu=J_{n\nu}=J_{p\nu}=0$", xy=(-1.9, 2.35), xytext=(-9.4, 3.7),
                 fontsize=6.25, arrowprops={"arrowstyle": "->", "lw": 0.65, "color": INK})
    geo.text(-9.6, -6.9, r"全域：$\mathbf{u}=0$，$\theta=0$", fontsize=6.25)
    geo.text(
        -9.6,
        -11.25,
        r"L2：$e_\phi=2.78\times10^{-4}$，$e_E=5.41\times10^{-3}$；PASS",
        fontsize=6.15,
        color=INK,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#AAB5BB", "lw": 0.6},
    )
    geo.legend(
        handles=[
            Line2D([], [], linestyle="none", marker="o", ms=4.2, mfc="white", mec=FEM,
                   label="中心线电势点"),
            Line2D([], [], linestyle="none", marker="o", ms=3.0, mfc=RED, mec=RED,
                   label="孔壁角向点"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.50, 0.12),
        ncol=2,
        columnspacing=1.1,
        handletextpad=0.4,
        fontsize=6.1,
    )

    # Centre-line potential comparison.
    x_left = np.linspace(-L_UM, -A_UM, 201) * 1e-6
    x_right = np.linspace(A_UM, L_UM, 201) * 1e-6
    for x, label in ((x_left, "解析解"), (x_right, None)):
        phi_ax.plot(x * 1e6, exact_phi(x), color=ANALYTIC, lw=1.2, marker="o",
                    markevery=20, ms=2.2, mfc=ANALYTIC, mec=ANALYTIC, label=label)
    phi_ax.plot(center["x_m"] * 1e6, center["phi_V"], linestyle="none", marker="o",
                ms=4.0, mfc="white", mec=FEM, mew=1.15, label="FEM（L2）")
    phi_ax.axvspan(-A_UM, A_UM, color="#E7ECEF", zorder=-2)
    phi_ax.text(0, 0, "圆孔", ha="center", va="center", color=MUTED)
    phi_ax.set(xlabel=r"中心线 $x$（$\mu$m），$y=z=0$", ylabel=r"电势 $\phi$（V）",
               title="中心线电势")
    phi_ax.legend(loc="upper right")

    # Tangential field plus a signed normal-field inset: the added comparison.
    et_exact = 2.0 * np.abs(np.sin(np.radians(exact_angle)))
    field_ax.plot(exact_angle, et_exact, color=ANALYTIC, lw=1.2, marker="o", markevery=36,
                  ms=2.2, mfc=ANALYTIC, mec=ANALYTIC, label="解析解")
    field_ax.plot(angle[shown], wall["Et_over_E0"][shown], linestyle="none", marker="o",
                  ms=3.4, mfc="white", mec=FEM, mew=1.05, label="FEM（L2）")
    field_ax.set(xlabel=r"孔壁角度 $\psi$（°）", ylabel=r"切向电场 $|E_t|/E_0$",
                 title="孔壁电场：切向集中与法向绝缘", xticks=[0, 90, 180, 270, 360])
    field_ax.legend(loc="upper center", bbox_to_anchor=(0.50, 0.98), ncol=2)

    inset = field_ax.inset_axes([0.57, 0.12, 0.39, 0.28])
    inset.plot(exact_angle, np.zeros_like(exact_angle), color=ANALYTIC, lw=0.8,
               marker="o", markevery=180, ms=1.7, mfc=ANALYTIC, mec=ANALYTIC)
    inset.plot(angle[shown], wall["En_over_E0"][shown], linestyle="none", marker="o",
               ms=2.3, mfc="white", mec=FEM, mew=0.8)
    en_lim = 1.20 * np.max(np.abs(wall["En_over_E0"]))
    inset.set_xlim(0, 360)
    inset.set_ylim(-en_lim, en_lim)
    inset.set_xticks([0, 180, 360])
    inset.tick_params(labelsize=5.3, length=2.0)
    inset.set_title(r"法向场 $E_n/E_0$；解析值 0", fontsize=5.8, pad=2)
    inset.text(0.03, 0.08, r"RMS=$8.76\times10^{-3}$", transform=inset.transAxes,
               fontsize=5.2, color=FEM)

    # Normalized field-squared Joule-source diagnostic.
    q_exact = 4.0 * np.sin(np.radians(exact_angle)) ** 2
    q_ax.plot(exact_angle, q_exact, color=ANALYTIC, lw=1.2, marker="o", markevery=36,
              ms=2.2, mfc=ANALYTIC, mec=ANALYTIC, label="解析解")
    q_ax.plot(angle[shown], wall["Q_over_Qinf"][shown], linestyle="none", marker="o",
              ms=3.4, mfc="white", mec=FEM, mew=1.05, label="FEM（L2）")
    q_ax.set(xlabel=r"孔壁角度 $\psi$（°）", ylabel=r"归一化焦耳热源 $Q_J/Q_\infty$",
             title="孔壁焦耳热源", xticks=[0, 90, 180, 270, 360])
    q_ax.legend(loc="upper center", bbox_to_anchor=(0.50, 0.66), ncol=2)

    for ax in (phi_ax, field_ax, q_ax):
        ax.tick_params(direction="in", length=3.0, width=0.7)
        ax.grid(False)

    fig.suptitle("A2：绝缘圆孔解析解—FEM 四联验证", fontsize=11, weight="bold", y=0.975)
    fig.text(
        0.50,
        0.020,
        "注：焦耳热曲线为恒电导条件下的电场平方诊断；温度固定，不表示温升。",
        ha="center",
        va="bottom",
        fontsize=6.6,
        color=MUTED,
    )
    save_all(fig, "A2_four_panel_validation")

    contract = {
        "title": "A2：绝缘圆孔解析解—FEM 四联验证",
        "panel_order": [
            "geometry, boundary conditions, and actual plotted sample sets",
            "centerline potential: analytical versus L2 FEM",
            "wall tangential field with an inset of signed wall-normal field versus exact zero",
            "normalized wall Joule-source diagnostic: analytical versus L2 FEM",
        ],
        "added_comparison": {
            "quantity": "signed hole-wall normal electric field En/E0",
            "analytical_value": 0.0,
            "reason": "directly tests the insulating-hole natural boundary and exposes angular discretization error",
            "L2_RMS_over_E0": values["a2_L2_wall_normal_E_RMS_over_E0"],
        },
        "evidence_source": "solver_v2/runs/F1_20260906_R0144/A2",
        "solver_rerun": False,
        "scope": "synthetic isothermal constant-conductivity A2 validation; not ZnO thermal prediction",
    }
    (OUT / "A2_four_panel_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def plot_reference_style_four_comparisons() -> None:
    """Match the user's 2x2 reference layout for the four requested A2 comparisons."""
    nodes = np.genfromtxt(RUN / "A2_L2_nodes.csv", delimiter=",", names=True)
    wall = np.genfromtxt(RUN / "A2_L2_wall.csv", delimiter=",", names=True)
    center_q, coordinate_difference = centerline_fem_joule_samples(nodes)

    center_nodes = nodes[
        (np.abs(nodes["y_m"]) < 1.0e-15) & (np.abs(nodes["z_m"]) < 1.0e-15)
    ]
    center_nodes.sort(order="x_m")
    center_xi = (center_nodes["x_m"] + L_UM * 1.0e-6) / (2.0 * L_UM * 1.0e-6)

    radius = np.hypot(nodes["x_m"], nodes["y_m"])
    wall_nodes = nodes[
        (np.abs(radius - A_UM * 1.0e-6) < 1.0e-13)
        & (np.abs(nodes["z_m"]) < 1.0e-15)
    ]
    wall_node_angle = np.mod(np.arctan2(wall_nodes["y_m"], wall_nodes["x_m"]), 2.0 * np.pi)
    wall_node_order = np.argsort(wall_node_angle)
    wall_nodes = wall_nodes[wall_node_order]
    wall_node_angle = wall_node_angle[wall_node_order]
    wall_node_show = np.arange(0, len(wall_nodes), 2, dtype=int)
    wall_phi_xi = np.r_[wall_node_angle[wall_node_show] / (2.0 * np.pi), 1.0]
    wall_phi_fem_mv = np.r_[wall_nodes["phi_V"][wall_node_show], wall_nodes["phi_V"][0]] * 1.0e3

    wall_angle = np.mod(wall["angle_rad"], 2.0 * np.pi)
    wall_order = np.argsort(wall_angle)
    wall = wall[wall_order]
    wall_xi_all = wall_angle[wall_order] / (2.0 * np.pi)
    wall_q_indices: list[int] = []
    for target in np.linspace(0.0, 1.0, 17):
        index = int(np.argmin(np.abs(wall_xi_all - target)))
        if index not in wall_q_indices:
            wall_q_indices.append(index)
    wall_q_indices_array = np.asarray(wall_q_indices, dtype=int)

    xi_hole_left = (L_UM - A_UM) / (2.0 * L_UM)
    xi_hole_right = (L_UM + A_UM) / (2.0 * L_UM)
    xi_left = np.linspace(0.0, xi_hole_left, 201)
    xi_right = np.linspace(xi_hole_right, 1.0, 201)
    x_left = (2.0 * xi_left - 1.0) * L_UM * 1.0e-6
    x_right = (2.0 * xi_right - 1.0) * L_UM * 1.0e-6
    wall_xi_exact = np.linspace(0.0, 1.0, 401)
    wall_angle_exact = 2.0 * np.pi * wall_xi_exact

    phi_left_mv = exact_phi(x_left) * 1.0e3
    phi_right_mv = exact_phi(x_right) * 1.0e3
    phi_wall_mv = -2.0 * E0 * A_UM * 1.0e-6 * np.cos(wall_angle_exact) * 1.0e3
    q_left = (1.0 - (A_UM * 1.0e-6 / x_left) ** 2) ** 2
    q_right = (1.0 - (A_UM * 1.0e-6 / x_right) ** 2) ** 2
    q_wall = 4.0 * np.sin(wall_angle_exact) ** 2

    fig, axes = plt.subplots(2, 2, figsize=(183 * MM, 134 * MM))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.095, top=0.975,
                        wspace=0.27, hspace=0.30)
    phi_center_ax, phi_wall_ax, q_center_ax, q_wall_ax = axes.ravel()

    analytical_style = {
        "color": REF_BLUE, "lw": 1.0, "marker": "o", "ms": 2.3,
        "mfc": REF_BLUE, "mec": REF_BLUE,
    }
    fem_style = {
        "linestyle": "none", "marker": "^", "ms": 3.8,
        "mfc": "white", "mec": REF_BLUE, "mew": 0.8,
    }

    phi_center_ax.plot(xi_left, phi_left_mv, markevery=20, label="解析解", **analytical_style)
    phi_center_ax.plot(xi_right, phi_right_mv, markevery=20, **analytical_style)
    phi_center_ax.plot(center_xi, center_nodes["phi_V"] * 1.0e3,
                       label="FEM（L2）", **fem_style)
    phi_center_ax.set_ylabel(r"中心线电势 $\phi$（mV）")
    phi_center_ax.set_xlabel(r"中心线位置 $x_1/L$")
    phi_center_ax.legend(loc="upper right")

    phi_wall_ax.plot(wall_xi_exact, phi_wall_mv, markevery=20,
                     label="解析解", **analytical_style)
    phi_wall_ax.plot(wall_phi_xi, wall_phi_fem_mv, label="FEM（L2）", **fem_style)
    phi_wall_ax.set_ylabel(r"孔壁电势 $\phi$（mV）")
    phi_wall_ax.set_xlabel(r"孔壁角度 $\psi/(2\pi)$")
    phi_wall_ax.legend(loc="upper left")

    q_center_ax.plot(xi_left, q_left, markevery=20, label="解析解", **analytical_style)
    q_center_ax.plot(xi_right, q_right, markevery=20, **analytical_style)
    q_center_ax.plot(center_q["x1_over_L"], center_q["Q_fem_over_Qinf"],
                     label="FEM（L2）", **fem_style)
    q_center_ax.set_ylabel(r"中心线焦耳热源 $Q_J/Q_\infty$")
    q_center_ax.set_xlabel(r"中心线位置 $x_1/L$")
    q_center_ax.legend(loc="upper center")

    q_wall_ax.plot(wall_xi_exact, q_wall, markevery=20,
                   label="解析解", **analytical_style)
    q_wall_ax.plot(wall_xi_all[wall_q_indices_array],
                   wall["Q_over_Qinf"][wall_q_indices_array],
                   label="FEM（L2）", **fem_style)
    q_wall_ax.set_ylabel(r"孔壁焦耳热源 $Q_J/Q_\infty$")
    q_wall_ax.set_xlabel(r"孔壁角度 $\psi/(2\pi)$")
    q_wall_ax.legend(loc="lower center")

    for ax in axes.ravel():
        ax.set_xlim(0.0, 1.0)
        ax.set_xticks(np.linspace(0.0, 1.0, 6))
        ax.grid(False)
        for side in ("top", "right", "bottom", "left"):
            ax.spines[side].set_visible(True)
            ax.spines[side].set_color("black")
            ax.spines[side].set_linewidth(0.75)
        ax.tick_params(axis="both", which="both", direction="in", top=True, right=True,
                       length=3.0, width=0.65, labelsize=7.0)
        ax.legend(frameon=False, fontsize=6.6, handlelength=2.2,
                  handletextpad=0.45, borderaxespad=0.35)
    phi_center_ax.set_ylim(-120.0, 120.0)
    phi_wall_ax.set_ylim(-70.0, 70.0)
    q_center_ax.set_ylim(-0.04, 0.92)
    q_wall_ax.set_ylim(-0.15, 4.25)
    save_all(fig, "A2_reference_style_four_comparisons")

    with (OUT / "A2_centerline_fem_joule_samples.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "x1_over_L", "x_m", "Q_fem_over_Qinf", "Q_exact_over_Qinf",
            "element_id", "r_natural", "side",
        ])
        for row in center_q:
            writer.writerow([
                row["x1_over_L"], row["x_m"], row["Q_fem_over_Qinf"],
                row["Q_exact_over_Qinf"], row["element_id"], row["r_natural"], row["side"],
            ])

    source_rows: list[list[object]] = []
    for xi, value in zip(np.r_[xi_left, xi_right], np.r_[phi_left_mv, phi_right_mv]):
        source_rows.append(["centerline_potential", xi, "analytic", value, "mV", True])
    for xi, value in zip(center_xi, center_nodes["phi_V"] * 1.0e3):
        source_rows.append(["centerline_potential", xi, "FEM_L2", value, "mV", True])
    for xi, value in zip(wall_xi_exact, phi_wall_mv):
        source_rows.append(["wall_potential", xi, "analytic", value, "mV", True])
    for xi, value in zip(wall_phi_xi, wall_phi_fem_mv):
        source_rows.append(["wall_potential", xi, "FEM_L2", value, "mV", True])
    for xi, value in zip(np.r_[xi_left, xi_right], np.r_[q_left, q_right]):
        source_rows.append(["centerline_joule", xi, "analytic", value, "QJ/Qinf", True])
    for row in center_q:
        source_rows.append([
            "centerline_joule", row["x1_over_L"], "FEM_L2",
            row["Q_fem_over_Qinf"], "QJ/Qinf", True,
        ])
    for xi, value in zip(wall_xi_exact, q_wall):
        source_rows.append(["wall_joule", xi, "analytic", value, "QJ/Qinf", True])
    displayed_wall_q = set(wall_q_indices)
    for index, (xi, value) in enumerate(zip(wall_xi_all, wall["Q_over_Qinf"])):
        source_rows.append([
            "wall_joule", xi, "FEM_L2", value, "QJ/Qinf", index in displayed_wall_q,
        ])
    with (OUT / "A2_reference_four_panel_source.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["panel", "normalized_coordinate", "series", "value", "unit", "displayed"])
        writer.writerows(source_rows)

    contract = {
        "title": "A2 analytical solution versus FEM — four requested comparisons",
        "reference_image": "codex-clipboard-0d52a753-f88f-4715-b4f4-e89705ad3a6a.png",
        "layout": {
            "figure_size_mm": [183, 134],
            "grid": "2x2",
            "panel_titles": False,
            "background": "white",
            "axes": "black four-sided box with inward ticks on all four sides",
            "analytical": "blue solid line with filled circular markers",
            "FEM": "unconnected blue open upward triangles",
            "normalized_horizontal_ranges": [0.0, 1.0],
        },
        "panels": [
            "centerline potential versus x1/L",
            "hole-wall potential versus psi/(2*pi)",
            "centerline normalized Joule source versus x1/L",
            "hole-wall normalized Joule source versus psi/(2*pi)",
        ],
        "analytical_expressions": {
            "centerline_potential": "phi=-E0*(x+a^2/x), abs(x)>=a",
            "wall_potential": "phi=-2*E0*a*cos(psi)",
            "centerline_joule": "QJ/Qinf=(1-a^2/x^2)^2, abs(x)>=a",
            "wall_joule": "QJ/Qinf=4*sin(psi)^2",
        },
        "centerline_FEM_joule_method": {
            "description": "evaluate E_h=-grad(phi_h) with the solver's HEX20 shape functions inside saved L2 elements; no nodal finite differencing",
            "natural_sampling": "r=+-1/sqrt(3), s=-1, t=0 on the theta=0 and theta=pi rays",
            "maximum_rebuilt_to_saved_node_coordinate_difference_m": coordinate_difference,
        },
        "displayed_comparison_maximum_absolute_differences": {
            "centerline_potential_mV": float(np.max(np.abs(
                center_nodes["phi_V"] - center_nodes["phi_exact_V"]
            )) * 1.0e3),
            "wall_potential_mV": float(np.max(np.abs(
                wall_nodes["phi_V"][wall_node_show] - wall_nodes["phi_exact_V"][wall_node_show]
            )) * 1.0e3),
            "centerline_QJ_over_Qinf": float(np.max(np.abs(
                center_q["Q_fem_over_Qinf"] - center_q["Q_exact_over_Qinf"]
            ))),
            "wall_QJ_over_Qinf": float(np.max(np.abs(
                wall["Q_over_Qinf"][wall_q_indices_array]
                - wall["Q_exact_over_Qinf"][wall_q_indices_array]
            ))),
        },
        "evidence": {
            "nodes": "solver_v2/runs/F1_20260906_R0144/A2/A2_L2_nodes.csv",
            "wall": "solver_v2/runs/F1_20260906_R0144/A2/A2_L2_wall.csv",
            "solver_rerun": False,
        },
        "scope": "synthetic isothermal constant-conductivity A2; QJ is a field-squared source diagnostic, not temperature rise",
    }
    (OUT / "A2_reference_four_panel_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def write_metadata(values: dict[str, float]) -> None:
    parameters = [
        ("outer_length_x", 20.0e-6, "m"),
        ("outer_length_y", 20.0e-6, "m"),
        ("thickness", 1.0e-6, "m"),
        ("hole_radius", 3.0e-6, "m"),
        ("far_field_E0", E0, "V/m"),
        ("reference_electron_density", 1.0e21, "1/m^3"),
        ("reference_hole_density", 1.0e21, "1/m^3"),
        ("permittivity", 1.0e-10, "F/m"),
        ("electron_mobility", 0.02, "m^2/(V s)"),
        ("hole_mobility", 0.005, "m^2/(V s)"),
    ]
    with (OUT / "A2_geometry_parameters.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["parameter", "value", "unit"])
        writer.writerows(parameters)

    boundaries = [
        ("all nodes", "u", "u_x=u_y=u_z=0", "essential"),
        ("all nodes", "theta", "theta=0", "essential"),
        ("four outer side faces", "phi", "phi=phi_exact(x,y)", "essential"),
        ("four outer side faces", "eta_n,eta_p", "eta_n=eta_p=0", "essential"),
        ("hole wall", "electric", "D dot nu=0", "natural"),
        ("hole wall", "carriers", "J_n dot nu=J_p dot nu=0", "natural"),
        ("top/bottom", "electric and carriers", "D_z=J_nz=J_pz=0", "natural"),
    ]
    with (OUT / "A2_boundary_conditions.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["region", "field", "condition", "type"])
        writer.writerows(boundaries)

    sources = [RUN / "A2_L2_nodes.csv", RUN / "A2_L2_wall.csv", RUN / "A2_export_stdout.txt"]
    contract = {
        "verification_id": "A2_F1_insulated_hole",
        "title_zh": "A2：三维挤出绝缘圆孔的恒电导电场集中解析验证",
        "title_en": "A2 — Analytical verification of constant-conductivity electric-field concentration around an extruded insulating circular hole",
        "program_baseline": "solver_v2 only",
        "solver_rerun": False,
        "evidence": [{"path": str(p.relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(p)} for p in sources],
        "analytical_expressions": {
            "potential": "phi=-E0*(r+a^2/r)*cos(psi)",
            "wall_tangential_field": "abs(E_t)/E0=2*abs(sin(psi))",
            "wall_normalized_joule_source": "Q_J/Q_inf=4*sin(psi)^2",
        },
        "geometry": {"outer_um": [20.0, 20.0, 1.0], "hole_radius_um": 3.0, "E0_V_per_m": E0},
        "finest_mesh_results": {
            "phi_relative_L2": values["a2_L2_phi_relative_L2"],
            "E_relative_L2": values["a2_L2_E_relative_L2"],
            "wall_normal_E_RMS_over_E0": values["a2_L2_wall_normal_E_RMS_over_E0"],
            "wall_Et_peak_over_E0": values["a2_L2_wall_Et_peak_over_E0"],
            "wall_Q_peak_over_Qinf": values["a2_L2_wall_Q_peak_over_Qinf"],
            "status": "PASS",
        },
        "figure_contract": {
            "backend": "Python/matplotlib",
            "analytical_role": "solid line plus filled circular markers; independently evaluated closed form",
            "fem_role": "unconnected open markers from saved solver output",
            "svg_text": "editable",
            "exports": ["SVG", "PDF", "PNG 450 dpi", "TIFF 600 dpi LZW"],
        },
        "scope_limits": [
            "synthetic isothermal material; not a ZnO calibration",
            "temperature and displacement are fixed",
            "normalized Q curve is the constant-conductivity field-squared diagnostic",
            "does not validate thermal feedback or real void dielectric/interface physics",
        ],
        "python_source": "solver_v2/scripts/plot_a2_skill_delivery.py",
    }
    (OUT / "A2_figure_contract.json").write_text(json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8")


def verify_exports() -> None:
    required = [
        "A2_geometry_boundary_conditions.svg", "A2_geometry_boundary_conditions.pdf",
        "A2_geometry_boundary_conditions.png", "A2_geometry_boundary_conditions.tiff",
        "A2_analytic_fem_results.svg", "A2_analytic_fem_results.pdf",
        "A2_analytic_fem_results.png", "A2_analytic_fem_results.tiff",
        "A2_centerline_comparison.csv", "A2_wall_comparison.csv",
        "A2_error_convergence.csv", "A2_geometry_parameters.csv",
        "A2_boundary_conditions.csv", "A2_figure_contract.json",
        "A2_reference_style_geometry_sampling.svg", "A2_reference_style_geometry_sampling.pdf",
        "A2_reference_style_geometry_sampling.png", "A2_reference_style_geometry_sampling.tiff",
        "A2_plot_sampling_points.csv", "A2_reference_geometry_sampling_contract.json",
        "A2_four_panel_validation.svg", "A2_four_panel_validation.pdf",
        "A2_four_panel_validation.png", "A2_four_panel_validation.tiff",
        "A2_four_panel_contract.json",
        "A2_reference_style_four_comparisons.svg", "A2_reference_style_four_comparisons.pdf",
        "A2_reference_style_four_comparisons.png", "A2_reference_style_four_comparisons.tiff",
        "A2_centerline_fem_joule_samples.csv", "A2_reference_four_panel_source.csv",
        "A2_reference_four_panel_contract.json",
    ]
    missing = [name for name in required if not (OUT / name).is_file()]
    if missing:
        raise RuntimeError(f"Missing exports: {missing}")
    for svg_name in ("A2_analytic_fem_results.svg", "A2_reference_style_four_comparisons.svg"):
        svg_text = (OUT / svg_name).read_text(encoding="utf-8")
        if "<text" not in svg_text or re.search(r"font-type:\s*path", svg_text):
            raise RuntimeError(f"SVG text editability check failed: {svg_name}")
    lines = [f"{sha256(OUT / name)}  {name}" for name in required]
    (OUT / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n", encoding="ascii")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    plot_geometry()
    plot_reference_sampling_geometry()
    values = plot_results()
    plot_four_panel_composite()
    plot_reference_style_four_comparisons()
    write_metadata(values)
    verify_exports()
    print(f"A2_SKILL_OUTPUT={OUT}")
    print("VERIFICATION_RESULT=PASS")


if __name__ == "__main__":
    main()
