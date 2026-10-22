"""Draw the editable A2 geometry, boundary-condition, mesh, and parameter figure."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, Polygon, Rectangle


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "figures" / "F1_20260906"
STEM = OUT / "A2_geometry_boundary_mesh_parameters"

MM = 1.0 / 25.4
Q = 1.602176634e-19
KB = 1.380649e-23
T0 = 300.0
HALF_WIDTH_UM = 10.0
RADIUS_UM = 3.0
THICKNESS_UM = 1.0
E0 = 1.0e4
N0 = 1.0e21
P0 = 1.0e21
MU_N = 0.02
MU_P = 0.005
D_N = KB * T0 / Q * MU_N
D_P = KB * T0 / Q * MU_P
SIGMA0 = Q * (N0 * MU_N + P0 * MU_P)
QINF = SIGMA0 * E0**2

COLORS = {
    "ink": "#26333D",
    "muted": "#66757F",
    "grid": "#A8B4BC",
    "light": "#EEF2F4",
    "blue": "#5B8DB8",
    "blue_light": "#DCEAF5",
    "red": "#C9675B",
    "red_light": "#F4DFDA",
    "teal": "#4E8E82",
    "teal_light": "#DDEDE9",
    "orange": "#D58A45",
    "orange_light": "#F7E9D8",
    "purple": "#8574A7",
}


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


def clean_axis(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("auto")
    ax.axis("off")


def panel_header(ax, letter, title):
    ax.text(0.005, 0.975, letter, va="top", ha="left", fontsize=9.0, weight="bold", color=COLORS["ink"])
    ax.text(0.040, 0.975, title, va="top", ha="left", fontsize=8.25, weight="bold", color=COLORS["ink"])


def arrow(ax, start, end, color=None, style="-|>", lw=0.9, mutation=8, zorder=5):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle=style,
        mutation_scale=mutation,
        linewidth=lw,
        color=color or COLORS["ink"],
        shrinkA=0,
        shrinkB=0,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def callout(ax, xy, width, height, title, lines, facecolor, edgecolor, fontsize=5.8):
    box = Rectangle(xy, width, height, facecolor=facecolor, edgecolor=edgecolor, linewidth=0.75, zorder=6)
    ax.add_patch(box)
    ax.text(
        xy[0] + 0.012,
        xy[1] + height - 0.026,
        title,
        ha="left",
        va="top",
        fontsize=6.45,
        weight="bold",
        color=COLORS["ink"],
        zorder=7,
    )
    ax.text(
        xy[0] + 0.012,
        xy[1] + height - 0.067,
        "\n".join(lines),
        ha="left",
        va="top",
        fontsize=fontsize,
        linespacing=1.28,
        color=COLORS["ink"],
        zorder=7,
    )
    return box


def project(u, v, z, origin=(0.36, 0.23), sx=0.31, sy=0.24, dx=0.075, dy=0.09):
    return (
        origin[0] + sx * u + dx * v,
        origin[1] + sy * z + dy * v,
    )


def draw_extruded_plate(ax):
    clean_axis(ax)
    panel_header(ax, "a", "A2 三维挤出绝缘圆孔与边界条件")

    top = [
        project(-1, -1, 0.55),
        project(1, -1, 0.55),
        project(1, 1, 0.55),
        project(-1, 1, 0.55),
    ]
    bottom = [
        project(-1, -1, -0.55),
        project(1, -1, -0.55),
        project(1, 1, -0.55),
        project(-1, 1, -0.55),
    ]

    ax.add_patch(Polygon([bottom[0], bottom[1], top[1], top[0]], closed=True,
                         facecolor="#E8EDF0", edgecolor=COLORS["ink"], lw=0.8, zorder=1))
    ax.add_patch(Polygon([bottom[1], bottom[2], top[2], top[1]], closed=True,
                         facecolor="#DCE5EA", edgecolor=COLORS["ink"], lw=0.8, zorder=1))
    ax.add_patch(Polygon(top, closed=True, facecolor="#F9FAFB", edgecolor=COLORS["ink"], lw=1.0, zorder=2))

    left_face = [top[0], top[3], bottom[3], bottom[0]]
    right_face = [top[1], top[2], bottom[2], bottom[1]]
    ax.add_patch(Polygon(left_face, closed=True, facecolor=COLORS["red_light"],
                         edgecolor=COLORS["red"], lw=1.0, alpha=0.92, zorder=3))
    ax.add_patch(Polygon(right_face, closed=True, facecolor=COLORS["blue_light"],
                         edgecolor=COLORS["blue"], lw=1.0, alpha=0.92, zorder=3))

    center_top = project(0, 0, 0.55)
    center_bottom = project(0, 0, -0.55)
    hole_w = 0.19
    hole_h = 0.065
    ax.add_patch(
        Polygon(
            [
                (center_top[0] - hole_w / 2, center_top[1]),
                (center_top[0] + hole_w / 2, center_top[1]),
                (center_bottom[0] + hole_w / 2, center_bottom[1]),
                (center_bottom[0] - hole_w / 2, center_bottom[1]),
            ],
            closed=True,
            facecolor="#CBD6DC",
            edgecolor="none",
            zorder=4,
        )
    )
    ax.add_patch(Ellipse(center_bottom, hole_w, hole_h, facecolor="#C4D0D6",
                         edgecolor=COLORS["ink"], lw=0.75, zorder=4))
    ax.add_patch(Ellipse(center_top, hole_w, hole_h, facecolor="white",
                         edgecolor=COLORS["ink"], lw=1.0, zorder=5))
    ax.plot(
        [center_top[0] - hole_w / 2, center_bottom[0] - hole_w / 2],
        [center_top[1], center_bottom[1]],
        color=COLORS["ink"],
        lw=0.75,
        zorder=5,
    )
    ax.plot(
        [center_top[0] + hole_w / 2, center_bottom[0] + hole_w / 2],
        [center_top[1], center_bottom[1]],
        color=COLORS["ink"],
        lw=0.75,
        zorder=5,
    )

    arrow(ax, (0.345, 0.735), (0.665, 0.735), COLORS["orange"], lw=1.5, mutation=10)
    ax.text(0.505, 0.775, r"远场 $\mathbf{E}_0=E_0\mathbf{e}_x$", ha="center", va="bottom",
            fontsize=7.0, weight="bold", color=COLORS["orange"])

    callout(
        ax,
        (0.012, 0.255),
        0.250,
        0.315,
        "四个外侧面",
        [
            r"$\phi=\phi_{\rm ex}(x,y)$",
            r"$\eta_n=\eta_p=0$",
            r"$\mathbf{u}=0,\ \theta=0$",
            "空间变化解析电势施加在完整外周",
        ],
        COLORS["red_light"],
        COLORS["red"],
    )
    callout(
        ax,
        (0.745, 0.270),
        0.245,
        0.300,
        "贯穿圆柱孔壁  r = a",
        [
            r"$\mathbf{D}\cdot\boldsymbol{\nu}=0$",
            r"$\mathbf{J}_n\cdot\boldsymbol{\nu}=0$",
            r"$\mathbf{J}_p\cdot\boldsymbol{\nu}=0$",
            r"$\mathbf{u}=0,\ \theta=0$",
        ],
        COLORS["blue_light"],
        COLORS["blue"],
    )
    callout(
        ax,
        (0.330, 0.035),
        0.395,
        0.155,
        "上下表面  z = ±t/2",
        [
            r"$\mathbf{D}\cdot\boldsymbol{\nu}=\mathbf{J}_n\cdot\boldsymbol{\nu}=\mathbf{J}_p\cdot\boldsymbol{\nu}=0$",
            r"$\mathbf{u}=0,\ \theta=0$；厚度方向仅保留三维离散",
        ],
        COLORS["light"],
        COLORS["grid"],
        fontsize=5.65,
    )

    arrow(ax, (0.262, 0.410), (0.330, 0.410), COLORS["red"], lw=0.85)
    arrow(ax, (0.745, 0.410), (center_top[0] + 0.085, center_top[1] + 0.01), COLORS["blue"], lw=0.85)
    arrow(ax, (0.530, 0.190), (0.530, 0.235), COLORS["muted"], lw=0.75)

    ax.text(0.515, 0.628, r"$20\times20\times1\ \mu{\rm m}$ 方板", ha="center", va="center",
            fontsize=7.0, weight="bold", color=COLORS["ink"])
    ax.text(0.515, 0.588, r"中心贯穿孔：$a=3\ \mu{\rm m}$", ha="center", va="center",
            fontsize=6.6, color=COLORS["ink"])
def draw_top_view(ax):
    clean_axis(ax)
    panel_header(ax, "b", "俯视解析场、孔壁热点与停滞点")

    x0, y0, size = 0.13, 0.13, 0.67
    ax.add_patch(Rectangle((x0, y0), size, size, facecolor="#FBFCFC", edgecolor=COLORS["ink"], lw=1.0))
    center = (x0 + size / 2, y0 + size / 2)
    radius = size * RADIUS_UM / (2 * HALF_WIDTH_UM)
    ax.add_patch(Circle(center, radius, facecolor="white", edgecolor=COLORS["ink"], lw=1.1, zorder=4))

    for yy in np.linspace(y0 + 0.10, y0 + size - 0.10, 5):
        arrow(ax, (x0 + 0.055, yy), (x0 + size - 0.055, yy), COLORS["orange"], lw=0.65, mutation=6, zorder=2)

    for sign in (-1, 1):
        hot = (center[0], center[1] + sign * radius)
        ax.add_patch(Circle(hot, 0.026, facecolor=COLORS["red"], edgecolor="white", lw=0.65, zorder=7))
        ax.text(hot[0] + 0.040, hot[1], "热点", va="center", ha="left", fontsize=5.7,
                color=COLORS["red"], weight="bold")
    for sign in (-1, 1):
        stop = (center[0] + sign * radius, center[1])
        ax.add_patch(Circle(stop, 0.015, facecolor="white", edgecolor=COLORS["blue"], lw=1.0, zorder=7))
    ax.text(center[0], center[1] + radius + 0.075, r"$|E_\psi|/E_0=2,\ Q_J/Q_\infty=4$",
            ha="center", va="bottom", fontsize=5.9, color=COLORS["red"], weight="bold")
    ax.text(center[0], center[1] - radius - 0.070, r"$\psi=270^\circ$", ha="center", va="top",
            fontsize=5.6, color=COLORS["muted"])
    ax.text(center[0] - radius - 0.025, center[1] - 0.050, "停滞点", ha="right", va="top",
            fontsize=5.5, color=COLORS["blue"])
    ax.text(center[0] + radius + 0.025, center[1] - 0.050, "停滞点", ha="left", va="top",
            fontsize=5.5, color=COLORS["blue"])

    ax.text(x0 - 0.035, center[1], r"$\phi>0$", rotation=90, ha="center", va="center",
            fontsize=6.2, color=COLORS["red"], weight="bold")
    ax.text(x0 + size + 0.035, center[1], r"$\phi<0$", rotation=90, ha="center", va="center",
            fontsize=6.2, color=COLORS["blue"], weight="bold")
    ax.text(center[0], y0 - 0.050, r"$x,y\in[-10,10]\ \mu{\rm m}$", ha="center", va="top",
            fontsize=6.0, color=COLORS["ink"])

    origin = (0.80, 0.085)
    arrow(ax, origin, (0.89, 0.085), COLORS["red"], lw=1.0)
    arrow(ax, origin, (0.80, 0.175), COLORS["teal"], lw=1.0)
    ax.text(0.902, 0.083, "x", fontsize=6.4, color=COLORS["red"], va="center")
    ax.text(0.797, 0.188, "y", fontsize=6.4, color=COLORS["teal"], ha="center")

    ax.text(
        0.50,
        0.895,
        r"$\phi_{\rm ex}=-E_0\left(r+a^2/r\right)\cos\psi$",
        ha="center",
        va="top",
        fontsize=6.8,
        color=COLORS["ink"],
    )
    ax.text(
        0.50,
        0.850,
        r"$E_r(a,\psi)=0,\quad |E_\psi(a,\psi)|/E_0=2|\sin\psi|$",
        ha="center",
        va="top",
        fontsize=6.3,
        color=COLORS["ink"],
    )


def map_radius(theta, fraction):
    outer = 1.0 / np.maximum(np.abs(np.cos(theta)), np.abs(np.sin(theta)))
    return RADIUS_UM / HALF_WIDTH_UM + fraction * (outer - RADIUS_UM / HALF_WIDTH_UM)


def draw_mesh_at(ax, center, width, nt, nr, label, nodes):
    cx, cy = center
    theta_fine = np.linspace(0, 2 * np.pi, 801)
    mesh_ax = ax.inset_axes([cx - width / 2, cy - width / 2, width, width])
    mesh_ax.set_xlim(-1.03, 1.03)
    mesh_ax.set_ylim(-1.03, 1.03)
    mesh_ax.set_aspect("equal", adjustable="box")
    mesh_ax.axis("off")
    mesh_ax.add_patch(Rectangle((-1, -1), 2, 2, facecolor="white",
                                edgecolor=COLORS["ink"], lw=0.75))
    for er in range(nr + 1):
        fraction = (er / nr) ** 2
        rr = map_radius(theta_fine, fraction)
        mesh_ax.plot(rr * np.cos(theta_fine), rr * np.sin(theta_fine),
                     color=COLORS["ink"] if er in (0, nr) else COLORS["grid"],
                     lw=0.72 if er in (0, nr) else 0.43)
    for it in range(nt):
        theta = 2 * np.pi * it / nt
        fractions = np.linspace(0, 1, 80)
        rr = map_radius(theta, fractions)
        mesh_ax.plot(rr * np.cos(theta), rr * np.sin(theta), color=COLORS["grid"], lw=0.34)
    ax.text(cx, cy - width / 2 - 0.055, label, ha="center", va="top", fontsize=5.9,
            weight="bold", color=COLORS["ink"])
    ax.text(cx, cy - width / 2 - 0.098, f"{nt * nr} HEX20；{nodes} 节点", ha="center", va="top",
            fontsize=5.35, color=COLORS["muted"])


def draw_mesh_panel(ax):
    clean_axis(ax)
    panel_header(ax, "c", "周向×径向三层曲边 HEX20 网格")

    draw_mesh_at(ax, (0.17, 0.56), 0.25, 16, 2, "L0：16×2", 304)
    draw_mesh_at(ax, (0.50, 0.56), 0.25, 24, 3, "L1：24×3", 624)
    draw_mesh_at(ax, (0.83, 0.56), 0.25, 32, 4, "L2：32×4", 1056)

    ax.add_patch(Rectangle((0.055, 0.080), 0.890, 0.205, facecolor=COLORS["light"],
                           edgecolor=COLORS["grid"], lw=0.7))
    ax.text(0.075, 0.245, "网格规则", ha="left", va="top", fontsize=6.2,
            weight="bold", color=COLORS["ink"])
    ax.text(
        0.075,
        0.205,
        "• 厚度方向保持 1 个二次单元层；解析场严格不随 z 变化\n"
        "• 径向采用平方分级，使孔边分辨率高于外边界区域\n"
        "• 孔壁由二次曲面近似，几何误差与场误差随网格共同下降",
        ha="left",
        va="top",
        fontsize=5.75,
        linespacing=1.35,
        color=COLORS["ink"],
    )


def draw_parameter_panel(ax):
    clean_axis(ax)
    panel_header(ax, "d", "验证专用参数、求解设置与解析预期")

    ax.add_patch(Rectangle((0.025, 0.820), 0.950, 0.080, facecolor=COLORS["orange_light"],
                           edgecolor=COLORS["orange"], lw=0.8))
    ax.text(0.500, 0.860, "合成、各向同性、等温参数；不代表 ZnO 实验材料",
            ha="center", va="center", fontsize=6.25, weight="bold", color=COLORS["ink"])

    rows = [
        (r"几何", r"$L_x=L_y=20\ \mu$m；$t=1\ \mu$m；$a=3\ \mu$m"),
        (r"远场", r"$E_0=1.0\times10^4$ V m$^{-1}$"),
        (r"介电", r"$\boldsymbol{\epsilon}=1.0\times10^{-10}\mathbf{I}$ F m$^{-1}$"),
        (r"载流子", r"$n_0=p_0=1.0\times10^{21}$ m$^{-3}$"),
        (r"迁移率", r"$\mu_n=0.020$；$\mu_p=0.005$ m$^2$ V$^{-1}$ s$^{-1}$"),
        (r"扩散", r"$D_n=5.1704\times10^{-4}$；$D_p=1.2926\times10^{-4}$ m$^2$ s$^{-1}$"),
        (r"常电导", r"$\sigma_0=4.00544$ S m$^{-1}$"),
        (r"远场热源", r"$Q_\infty=\sigma_0E_0^2=4.00544\times10^8$ W m$^{-3}$"),
        (r"关闭耦合", r"$\mathbf{e}=\boldsymbol{\alpha}=\mathbf{p}^\epsilon=0$；$R=0$"),
        (r"固定场", r"$\mathbf{u}=0$；$\theta=T-T_0=0$，$T_0=300$ K"),
        (r"Newton", r"最大 30 次；残差/更新容限均为 $10^{-11}$"),
    ]

    top = 0.785
    row_h = 0.058
    split = 0.285
    ax.plot([split, split], [top - row_h * len(rows), top], color=COLORS["light"], lw=0.7)
    for i, (name, value) in enumerate(rows):
        y_top = top - i * row_h
        y_mid = y_top - row_h / 2
        ax.plot([0.025, 0.975], [y_top, y_top], color=COLORS["grid"], lw=0.42)
        ax.text(0.040, y_mid, name, ha="left", va="center", fontsize=5.55,
                weight="bold", color=COLORS["ink"])
        ax.text(split + 0.015, y_mid, value, ha="left", va="center", fontsize=5.25,
                color=COLORS["ink"])
    ax.plot([0.025, 0.975], [top - row_h * len(rows), top - row_h * len(rows)],
            color=COLORS["grid"], lw=0.42)

    ax.text(0.025, 0.085, r"解析预测：孔壁上下侧 $|E_\psi|/E_0=2$，$Q_J/Q_\infty=4$；左右侧均为 0",
            ha="left", va="center", fontsize=5.55, color=COLORS["red"], weight="bold")
    ax.text(0.025, 0.035, r"焦耳热在本算例中用于后处理核验；温度自由度全部固定，不产生温升反馈。",
            ha="left", va="center", fontsize=5.45, color=COLORS["muted"])


def write_supporting_files():
    OUT.mkdir(parents=True, exist_ok=True)
    parameters = [
        ("plate_length_x", 20.0e-6, "m", "square outer boundary"),
        ("plate_length_y", 20.0e-6, "m", "square outer boundary"),
        ("plate_thickness", 1.0e-6, "m", "one HEX20 layer"),
        ("through_hole_radius", 3.0e-6, "m", "centered at x=y=0"),
        ("far_field_E0", E0, "V/m", "along +x"),
        ("reference_temperature", T0, "K", "theta fixed to zero"),
        ("permittivity_isotropic", 1.0e-10, "F/m", "diagonal tensor entries"),
        ("reference_electron_density", N0, "1/m^3", "synthetic neutral material"),
        ("reference_hole_density", P0, "1/m^3", "synthetic neutral material"),
        ("intrinsic_concentration", 1.0e21, "1/m^3", "n0 p0 = ni^2"),
        ("donor_concentration", 0.0, "1/m^3", "default"),
        ("acceptor_concentration", 0.0, "1/m^3", "default"),
        ("electron_mobility", MU_N, "m^2/(V s)", "isotropic"),
        ("hole_mobility", MU_P, "m^2/(V s)", "isotropic"),
        ("electron_diffusion_300K", D_N, "m^2/s", "Einstein relation"),
        ("hole_diffusion_300K", D_P, "m^2/s", "Einstein relation"),
        ("conductivity_sigma0", SIGMA0, "S/m", "q(n0 mu_n + p0 mu_p)"),
        ("far_field_joule_heat_Qinf", QINF, "W/m^3", "sigma0 E0^2"),
        ("dummy_elastic_diagonal", 1.0e11, "Pa", "mechanics fully constrained"),
        ("dummy_thermal_conductivity", 100.0, "W/(m K)", "temperature fully constrained"),
        ("piezoelectric", 0.0, "C/m^2", "disabled"),
        ("thermal_stress", 0.0, "Pa/K", "disabled"),
        ("pyroelectric", 0.0, "C/(m^2 K)", "disabled"),
        ("recombination", 0.0, "1/(m^3 s)", "none"),
        ("joule_factor", 1.0, "-", "post-processing/source evaluation active"),
        ("newton_max_iterations", 30, "-", "solver setting"),
        ("newton_residual_tolerance", 1.0e-11, "-", "solver setting"),
        ("newton_update_tolerance", 1.0e-11, "-", "solver setting"),
    ]
    with (OUT / "A2_geometry_model_parameters.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["parameter", "value", "unit", "note"])
        writer.writerows(parameters)

    mesh_rows = [
        ("L0", 16, 2, 1, 32, 304),
        ("L1", 24, 3, 1, 72, 624),
        ("L2", 32, 4, 1, 128, 1056),
    ]
    with (OUT / "A2_mesh_levels.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["level", "circumferential_elements", "radial_elements",
                         "through_thickness_elements", "total_HEX20_elements", "nodes"])
        writer.writerows(mesh_rows)

    boundary_rows = [
        ("all nodes", "mechanics", "u_x=u_y=u_z=0", "essential", "removes mechanical response from A2"),
        ("all nodes", "thermal", "theta=0", "essential", "isothermal; Joule heat cannot create temperature feedback"),
        ("four outer side faces", "electric", "phi=phi_ex(x,y)", "essential", "spatially varying exact potential on the full outer perimeter"),
        ("four outer side faces", "electron", "eta_n=0", "essential", "exact constant electron density"),
        ("four outer side faces", "hole", "eta_p=0", "essential", "exact constant hole density"),
        ("cylindrical hole wall", "electric", "D dot nu=0", "natural", "ideal insulating cavity"),
        ("cylindrical hole wall", "electron", "J_n dot nu=0", "natural", "no electron current through the hole"),
        ("cylindrical hole wall", "hole", "J_p dot nu=0", "natural", "no hole current through the hole"),
        ("top and bottom faces", "electric", "D dot nu=0", "natural", "z-independent exact field"),
        ("top and bottom faces", "electron", "J_n dot nu=0", "natural", "z-independent exact field"),
        ("top and bottom faces", "hole", "J_p dot nu=0", "natural", "z-independent exact field"),
    ]
    with (OUT / "A2_boundary_conditions.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["region", "field", "condition", "type", "purpose"])
        writer.writerows(boundary_rows)

    function_rows = [
        ("curved_geometry", "curved HEX20 approximation of a cylindrical through-hole", "covered"),
        ("spatial_dirichlet", "nonuniform exact potential on all four outer side faces", "covered"),
        ("natural_insulation", "zero normal electric displacement and carrier currents on the hole wall", "covered"),
        ("two_dimensional_field_in_three_dimensions", "z-invariant solution in an extruded 3D mesh", "covered"),
        ("electric_field_recovery", "E=-grad(phi) and relative volume error", "covered"),
        ("wall_field_decomposition", "normal and tangential electric-field components", "covered"),
        ("peak_location", "hot spots at 90 and 270 degrees; stagnation at 0 and 180 degrees", "covered"),
        ("normalized_joule_pattern", "Q_J/Q_inf=|E|^2/E0^2 for constant conductivity", "covered"),
        ("mesh_refinement", "L0/L1/L2 decrease of potential, field, and wall-normal errors", "covered"),
        ("temperature_feedback", "temperature-dependent mobility, diffusion, and self-heating", "not covered"),
        ("thermo_piezoelectric_coupling", "thermal expansion, piezoelectricity, and pyroelectricity", "not covered"),
        ("experimental_material_prediction", "quantitative ZnO specimen prediction", "not covered"),
    ]
    with (OUT / "A2_function_coverage.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["item", "meaning", "status"])
        writer.writerows(function_rows)

    contract = {
        "date": "2026-09-07",
        "core_conclusion": (
            "A2 isolates the geometric and electrical consequences of a through insulating circular hole, "
            "with an exact finite-domain boundary field and known 2x electric-field and 4x Joule-heating peaks."
        ),
        "archetype": "schematic-led composite",
        "backend": "Python/matplotlib only",
        "final_size_mm": [183, 136],
        "panel_map": {
            "a": "3D extruded square plate, through hole, coordinate convention, and boundary-condition classes",
            "b": "top-view analytical field, hot spots, stagnation points, and closed-form wall predictions",
            "c": "16x2, 24x3, and 32x4 curved HEX20 mesh sequence",
            "d": "synthetic material parameters, solver settings, and analytical expectations",
        },
        "source_data": [
            "frozen A2 implementation",
            "A2_geometry_model_parameters.csv",
            "A2_mesh_levels.csv",
            "A2_boundary_conditions.csv",
            "A2_function_coverage.csv",
        ],
        "statistics": "Not applicable; deterministic analytical verification",
        "reviewer_risks": [
            "A2 parameters are synthetic and must not be presented as ZnO data.",
            "The exact potential must be prescribed on all four outer side faces, not only left and right.",
            "Joule heat is evaluated while temperature is fixed; A2 does not validate thermal feedback.",
            "The 2x and 4x peaks are properties of the ideal constant-conductivity insulated-circle problem.",
        ],
        "editability": "SVG text remains editable; PDF uses TrueType fonts; schematic objects are vector paths.",
    }
    (OUT / "A2_geometry_model_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(183 * MM, 136 * MM))
    gs = fig.add_gridspec(
        2,
        12,
        height_ratios=[1.10, 1.0],
        hspace=0.085,
        wspace=0.16,
        left=0.025,
        right=0.985,
        top=0.975,
        bottom=0.035,
    )
    ax_a = fig.add_subplot(gs[0, :7])
    ax_b = fig.add_subplot(gs[0, 7:])
    ax_c = fig.add_subplot(gs[1, :7])
    ax_d = fig.add_subplot(gs[1, 7:])

    draw_extruded_plate(ax_a)
    draw_top_view(ax_b)
    draw_mesh_panel(ax_c)
    draw_parameter_panel(ax_d)

    fig.savefig(STEM.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.015)
    fig.savefig(STEM.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.015)
    fig.savefig(STEM.with_suffix(".png"), dpi=450, bbox_inches="tight", pad_inches=0.015)
    plt.close(fig)
    write_supporting_files()
    print(f"FIGURE_OK={STEM}")


if __name__ == "__main__":
    main()
