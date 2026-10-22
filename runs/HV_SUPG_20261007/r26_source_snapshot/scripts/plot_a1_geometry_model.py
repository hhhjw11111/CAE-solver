"""Draw the editable A1 geometry, boundary-condition, mesh, and parameter figure."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "figures" / "A1_100NM_010NM_NE20_DT10K_20260907"
STEM = OUT / "A1_geometry_boundary_mesh_parameters"

MM = 1.0 / 25.4
COLORS = {
    "ink": "#25313B",
    "muted": "#687681",
    "grid": "#AAB5BD",
    "light": "#EEF2F4",
    "blue": "#5B8DB8",
    "blue_light": "#DCEAF5",
    "red": "#C9675B",
    "red_light": "#F4DFDA",
    "teal": "#4E8E82",
    "orange": "#D58A45",
    "orange_light": "#F7E9D8",
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
    ax.text(0.040, 0.975, title, va="top", ha="left", fontsize=8.4, weight="bold", color=COLORS["ink"])


def callout(ax, xy, width, height, title, lines, facecolor, edgecolor):
    box = Rectangle(xy, width, height, facecolor=facecolor, edgecolor=edgecolor, linewidth=0.8, zorder=5)
    ax.add_patch(box)
    ax.text(xy[0] + 0.012, xy[1] + height - 0.028, title, ha="left", va="top", fontsize=6.7,
            weight="bold", color=COLORS["ink"], zorder=6)
    ax.text(xy[0] + 0.012, xy[1] + height - 0.072, "\n".join(lines), ha="left", va="top",
            fontsize=6.15, linespacing=1.32, color=COLORS["ink"], zorder=6)
    return box


def arrow(ax, start, end, color=None, style="-|>", lw=0.9, mutation=8, zorder=4):
    patch = FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=mutation, linewidth=lw,
                            color=color or COLORS["ink"], shrinkA=0, shrinkB=0, zorder=zorder)
    ax.add_patch(patch)
    return patch


def draw_prism(ax, left=0.245, bottom=0.235, length=0.47, height=0.31, dx=0.065, dy=0.095):
    p0 = (left, bottom)
    p1 = (left + length, bottom)
    p2 = (left + length, bottom + height)
    p3 = (left, bottom + height)
    q0 = (p0[0] + dx, p0[1] + dy)
    q1 = (p1[0] + dx, p1[1] + dy)
    q2 = (p2[0] + dx, p2[1] + dy)
    q3 = (p3[0] + dx, p3[1] + dy)

    ax.add_patch(Polygon([p3, p2, q2, q3], closed=True, facecolor="#F7F9FA", edgecolor=COLORS["ink"], lw=0.9))
    ax.add_patch(Polygon([p0, p1, p2, p3], closed=True, facecolor="white", edgecolor=COLORS["ink"], lw=1.0))
    ax.add_patch(Polygon([p0, p3, q3, q0], closed=True, facecolor=COLORS["blue_light"],
                         edgecolor=COLORS["blue"], lw=1.2, zorder=3))
    ax.add_patch(Polygon([p1, q1, q2, p2], closed=True, facecolor=COLORS["red_light"],
                         edgecolor=COLORS["red"], lw=1.2, zorder=3))
    for a, b in [(q0, q1), (q1, q2), (q2, q3), (q3, q0), (p0, q0), (p1, q1), (p2, q2), (p3, q3)]:
        ax.plot([a[0], b[0]], [a[1], b[1]], color=COLORS["ink"], lw=0.8, zorder=2)
    return {"p0": p0, "p1": p1, "p2": p2, "p3": p3, "q0": q0, "q1": q1, "q2": q2, "q3": q3}


def draw_geometry_panel(ax):
    clean_axis(ax)
    panel_header(ax, "a", "A1 三维棱柱与边界条件（解析场仅随 x 变化）")

    pts = draw_prism(ax)

    callout(
        ax, (0.012, 0.255), 0.205, 0.285, "左端  x = 0",
        [r"$u_x=0$", r"$\phi=0$", r"$\eta_n=\eta_p=0$", r"$\theta=0$"],
        COLORS["blue_light"], COLORS["blue"],
    )
    callout(
        ax, (0.790, 0.235), 0.198, 0.325, "右端  x = L",
        [r"$\sigma_{xx}=0$", r"$\phi=0$", r"$J_n\cdot\nu=J_p\cdot\nu=0$", r"$\theta=\Delta T=10$ K"],
        COLORS["red_light"], COLORS["red"],
    )
    callout(
        ax, (0.305, 0.705), 0.395, 0.195, "侧面与一维约束",
        [r"$D\cdot\nu=J_n\cdot\nu=J_p\cdot\nu=0$；热绝缘、零牵引",
         r"$u_y=u_z=0$；$\partial(\cdot)/\partial y=\partial(\cdot)/\partial z=0$"],
        COLORS["light"], COLORS["grid"],
    )

    arrow(ax, (0.217, 0.395), (0.246, 0.395), COLORS["blue"], lw=1.0)
    arrow(ax, (0.789, 0.395), (0.761, 0.395), COLORS["red"], lw=1.0)
    arrow(ax, (0.500, 0.704), (0.500, 0.640), COLORS["muted"], lw=0.8)

    ax.text(0.486, 0.438, r"$\Omega=[0,L]\times[0,b]\times[0,b]$", ha="center", va="center",
            fontsize=7.2, color=COLORS["ink"], weight="bold")
    ax.text(0.486, 0.370, r"$L=100$ nm, $b=10$ nm", ha="center", va="center", fontsize=7.0, color=COLORS["ink"])
    ax.text(0.486, 0.315, r"$\theta$: 0 → 10 K    $\phi$: 0 → 内部非零 → 0", ha="center", va="center",
            fontsize=6.6, color=COLORS["muted"])

    arrow(ax, (pts["p0"][0], 0.145), (pts["p1"][0], 0.145), COLORS["ink"], style="<->", lw=0.8)
    ax.text((pts["p0"][0] + pts["p1"][0]) / 2, 0.112, r"$L=100$ nm", ha="center", va="center",
            fontsize=6.4, color=COLORS["ink"])
    arrow(ax, (pts["p0"][0] - 0.027, pts["p0"][1]), (pts["p3"][0] - 0.027, pts["p3"][1]),
          COLORS["ink"], style="<->", lw=0.8)
    ax.text(pts["p0"][0] - 0.050, (pts["p0"][1] + pts["p3"][1]) / 2, "$b=10$ nm", ha="center", va="center",
            rotation=90, fontsize=6.4, color=COLORS["ink"])

    origin = (0.720, 0.145)
    arrow(ax, origin, (0.775, 0.145), COLORS["red"], lw=1.0)
    arrow(ax, origin, (0.720, 0.215), COLORS["teal"], lw=1.0)
    arrow(ax, origin, (0.748, 0.185), COLORS["blue"], lw=1.0)
    ax.text(0.782, 0.142, "x", fontsize=6.3, color=COLORS["red"], va="center")
    ax.text(0.716, 0.225, "y", fontsize=6.3, color=COLORS["teal"], ha="center")
    ax.text(0.754, 0.191, "z", fontsize=6.3, color=COLORS["blue"])

    ax.text(0.012, 0.055,
            r"验证对象：$u_x$、$\phi$、$\eta_n$、$\eta_p$、$\theta$；关闭焦耳热与复合（$Q_J=0$，$R=0$）",
            ha="left", va="center", fontsize=6.3, color=COLORS["muted"])


def prism_wireframe(ax, left, bottom, length, height, dx, dy, elements):
    p0 = (left, bottom)
    p1 = (left + length, bottom)
    p2 = (left + length, bottom + height)
    p3 = (left, bottom + height)
    q0 = (p0[0] + dx, p0[1] + dy)
    q1 = (p1[0] + dx, p1[1] + dy)
    q2 = (p2[0] + dx, p2[1] + dy)
    q3 = (p3[0] + dx, p3[1] + dy)
    for a, b in [(p0, p1), (p1, p2), (p2, p3), (p3, p0), (q0, q1), (q1, q2), (q2, q3),
                 (q3, q0), (p0, q0), (p1, q1), (p2, q2), (p3, q3)]:
        ax.plot([a[0], b[0]], [a[1], b[1]], color=COLORS["ink"], lw=0.75)
    for i in range(1, elements):
        t = i / elements
        a = (p0[0] + t * length, p0[1])
        b = (p3[0] + t * length, p3[1])
        c = (q0[0] + t * length, q0[1])
        d = (q3[0] + t * length, q3[1])
        ax.plot([a[0], b[0]], [a[1], b[1]], color=COLORS["grid"], lw=0.65)
        ax.plot([c[0], d[0]], [c[1], d[1]], color=COLORS["grid"], lw=0.65)
        ax.plot([b[0], d[0]], [b[1], d[1]], color=COLORS["grid"], lw=0.65)
        ax.plot([a[0], c[0]], [a[1], c[1]], color=COLORS["grid"], lw=0.55)
    return p0, p1, p2, p3, q0, q1, q2, q3


def draw_hex20_inset(ax):
    x0, y0, w, h, dx, dy = 0.745, 0.545, 0.170, 0.190, 0.035, 0.050
    pts = prism_wireframe(ax, x0, y0, w, h, dx, dy, 1)
    p0, p1, p2, p3, q0, q1, q2, q3 = pts
    corners = [p0, p1, p2, p3, q0, q1, q2, q3]
    edges = [(p0, p1), (p1, p2), (p2, p3), (p3, p0), (q0, q1), (q1, q2), (q2, q3), (q3, q0),
             (p0, q0), (p1, q1), (p2, q2), (p3, q3)]
    for x, y in corners:
        ax.add_patch(Circle((x, y), 0.007, facecolor=COLORS["ink"], edgecolor="white", lw=0.35, zorder=4))
    for a, b in edges:
        x, y = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
        ax.add_patch(Circle((x, y), 0.006, facecolor="white", edgecolor=COLORS["orange"], lw=0.8, zorder=4))
    ax.text(0.842, 0.785, "单个 HEX20", ha="center", va="bottom", fontsize=6.3, weight="bold", color=COLORS["ink"])
    ax.text(0.842, 0.505, "8 角点 + 12 棱中点", ha="center", va="top", fontsize=5.8, color=COLORS["muted"])


def segmented_bar(ax, x, y, width, height, n, color):
    ax.add_patch(Rectangle((x, y), width, height, facecolor="white", edgecolor=COLORS["ink"], lw=0.55))
    for i in range(1, n):
        xx = x + width * i / n
        ax.plot([xx, xx], [y, y + height], color=color, lw=0.42)


def draw_mesh_panel(ax):
    clean_axis(ax)
    panel_header(ax, "b", "HEX20 网格、收敛序列与 20 单元展示网格")

    prism_wireframe(ax, 0.070, 0.555, 0.555, 0.205, 0.070, 0.070, 20)
    ax.text(0.385, 0.835, "最终展示：Ne = 20，单元 5 nm × 10 nm × 10 nm", ha="center", va="bottom", fontsize=6.4,
            weight="bold", color=COLORS["ink"])
    ax.text(0.382, 0.505, "横截面 y、z 方向各 1 个单元；20 单元网格有 41 个轴向节点位置", ha="center", va="top",
            fontsize=6.1, color=COLORS["muted"])
    draw_hex20_inset(ax)

    ax.text(0.070, 0.420, "收敛序列", fontsize=6.5, weight="bold", color=COLORS["ink"], va="center")
    levels = [1, 2, 4, 8, 16]
    ys = [0.345, 0.285, 0.225, 0.165, 0.105]
    for n, y in zip(levels, ys):
        ax.text(0.073, y + 0.018, f"Ne={n:>2}", fontsize=5.8, color=COLORS["ink"], ha="left", va="center")
        segmented_bar(ax, 0.150, y, 0.700, 0.035, n, COLORS["blue"])
    ax.text(0.070, 0.035, "1/2/4/8/16 用于收敛阶；Ne=20 用于最终数值—解析曲线输出",
            fontsize=6.0, color=COLORS["muted"], ha="left", va="center")


def draw_parameter_panel(ax):
    clean_axis(ax)
    panel_header(ax, "c", "几何与合成验证参数")

    ax.add_patch(Rectangle((0.025, 0.855), 0.950, 0.085, facecolor=COLORS["orange_light"],
                           edgecolor=COLORS["orange"], lw=0.8))
    ax.text(0.500, 0.897, "验证专用参数，不代表 ZnO 实验材料", ha="center", va="center",
            fontsize=6.5, weight="bold", color=COLORS["ink"])

    left = [
        (r"$L$；$b$", r"$100$ nm；$10$ nm"),
        (r"$T_0$；$\Delta T$", "300 K；10 K"),
        (r"$C$", r"$1.0\times10^{11}$ Pa"),
        (r"$e$", r"$1.0$ C m$^{-2}$"),
        (r"$\epsilon$", r"$1.0\times10^{-8}$ F m$^{-1}$"),
        (r"$\alpha$", r"$1.0\times10^{-5}$ K$^{-1}$"),
        (r"$p^\epsilon$", r"$1.0\times10^{-8}$ C m$^{-2}$ K$^{-1}$"),
    ]
    right = [
        (r"$k$", r"$100$ W m$^{-1}$ K$^{-1}$"),
        (r"$n_0$；$p_0$", r"$1.0\times10^{21}$；$8.0\times10^{20}$ m$^{-3}$"),
        (r"$\mu_n$；$\mu_p$", r"$0.05$；$0.02$ m$^2$ V$^{-1}$ s$^{-1}$"),
        (r"$D_n$", r"$1.2926\times10^{-3}$ m$^2$ s$^{-1}$"),
        (r"$D_p$", r"$5.1704\times10^{-4}$ m$^2$ s$^{-1}$"),
        ("净掺杂", r"$n_0-p_0$（固定）"),
        ("热源；复合", r"$Q_J=0$；$R=0$"),
    ]

    def column(x0, width, rows):
        top = 0.800
        row_h = 0.094
        split = x0 + 0.30 * width
        ax.plot([split, split], [top - row_h * len(rows), top], color=COLORS["light"], lw=0.7)
        for i, (name, value) in enumerate(rows):
            y_top = top - i * row_h
            y_mid = y_top - row_h / 2
            ax.plot([x0, x0 + width], [y_top, y_top], color=COLORS["grid"], lw=0.45)
            ax.text(x0 + 0.012, y_mid, name, ha="left", va="center", fontsize=5.9,
                    weight="bold", color=COLORS["ink"])
            ax.text(split + 0.012, y_mid, value, ha="left", va="center", fontsize=5.55,
                    color=COLORS["ink"])
        y_bottom = top - row_h * len(rows)
        ax.plot([x0, x0 + width], [y_bottom, y_bottom], color=COLORS["grid"], lw=0.45)

    column(0.025, 0.455, left)
    column(0.520, 0.455, right)

    ax.text(0.025, 0.070, r"Einstein 关系：$D_s=k_BT_0\mu_s/q$；小扰动要求 $\max(|\eta_n|,|\eta_p|)\leq0.05$",
            ha="left", va="center", fontsize=5.75, color=COLORS["muted"])
    ax.text(0.025, 0.025, r"未知场：$u_x$、$\phi$、$\eta_n=\Delta n/n_0$、$\eta_p=\Delta p/p_0$、$\theta=T-T_0$",
            ha="left", va="center", fontsize=5.75, color=COLORS["muted"])


def write_supporting_files():
    parameters = [
        ("length_L", 1.0e-7, "m"),
        ("cross_section_width_b", 1.0e-8, "m"),
        ("presentation_axial_elements", 20, "HEX20"),
        ("reference_temperature_T0", 300.0, "K"),
        ("temperature_increment", 10.0, "K"),
        ("elastic_modulus_C", 1.0e11, "Pa"),
        ("piezoelectric_e", 1.0, "C/m^2"),
        ("permittivity", 1.0e-8, "F/m"),
        ("thermal_expansion", 1.0e-5, "1/K"),
        ("primary_pyroelectric", 1.0e-8, "C/(m^2 K)"),
        ("thermal_conductivity", 100.0, "W/(m K)"),
        ("reference_electron_density", 1.0e21, "1/m^3"),
        ("reference_hole_density", 8.0e20, "1/m^3"),
        ("electron_mobility", 0.05, "m^2/(V s)"),
        ("hole_mobility", 0.02, "m^2/(V s)"),
        ("electron_diffusion", 1.2926e-3, "m^2/s"),
        ("hole_diffusion", 5.1704e-4, "m^2/s"),
        ("joule_heating", 0.0, "disabled"),
        ("recombination", 0.0, "disabled"),
    ]
    with (OUT / "A1_geometry_model_parameters.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["parameter", "value", "unit"])
        writer.writerows(parameters)

    contract = {
        "date": "2026-09-07",
        "core_conclusion": "A1 is a three-dimensional HEX20 discretization of an analytically one-dimensional thermo-piezoelectric two-carrier verification problem.",
        "archetype": "schematic-led composite",
        "backend": "Python/matplotlib only",
        "final_size_mm": [183, 124],
        "panel_map": {
            "a": "3D prism, coordinate convention, one-dimensional restriction, and boundary conditions",
            "b": "HEX20 topology, the 1/2/4/8/16 convergence sequence, and the 20-element presentation mesh",
            "c": "synthetic geometry, material, carrier, and disabled-source parameters",
        },
        "source_data": "A1_geometry_model_parameters.csv and the frozen A1 implementation",
        "statistics": "Not applicable; deterministic discretization verification",
        "reviewer_risks": [
            "Synthetic verification constants must not be described as ZnO data.",
            "The computational coordinates use y,z in [0,b]; centered coordinates are an equivalent translation only.",
            "Only the x displacement is active and all exact fields depend on x alone.",
        ],
        "editability": "SVG text is retained as text; PDF uses TrueType fonts; all lines and patches are vector objects.",
    }
    (OUT / "A1_geometry_model_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(183 * MM, 124 * MM))
    gs = fig.add_gridspec(2, 12, height_ratios=[1.16, 1.0], hspace=0.075, wspace=0.12,
                          left=0.025, right=0.985, top=0.975, bottom=0.035)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :7])
    ax_c = fig.add_subplot(gs[1, 7:])

    draw_geometry_panel(ax_a)
    draw_mesh_panel(ax_b)
    draw_parameter_panel(ax_c)

    fig.savefig(STEM.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.015)
    fig.savefig(STEM.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.015)
    fig.savefig(STEM.with_suffix(".png"), dpi=450, bbox_inches="tight", pad_inches=0.015)
    fig.savefig(STEM.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pad_inches=0.015)
    plt.close(fig)
    write_supporting_files()
    print(f"FIGURE_OK={STEM}")


if __name__ == "__main__":
    main()
