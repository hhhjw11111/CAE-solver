"""Plot A1 analytical/FEM comparisons for three imposed temperature differences.

The A1 M0 system is strictly linear.  The 10 K FEM solution is the archived
solver result; the 5 K and 15 K FEM marker values are its exact linear load
scalings.  Analytical curves are evaluated independently from the frozen
closed-form equations.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE_RUN = ROOT / "runs" / "A1_100NM_010NM_NE20_DT10K_20260907"
OUT = ROOT.parent / "output" / "figures" / "A1_THREE_DT_005_010_015K_20260908"
OUT.mkdir(parents=True, exist_ok=True)

DELTA_T_VALUES = (5.0, 10.0, 15.0)
BASE_DELTA_T = 10.0

# Frozen A1 verification coefficients (SI units).
Q = 1.602176634e-19
K_B = 1.380649e-23
T_0 = 300.0
LENGTH = 100.0e-9
ELASTIC_C = 1.0e11
PIEZO_E = 1.0
PERMITTIVITY = 1.0e-8
REFERENCE_N = 1.0e21
REFERENCE_P = 8.0e20
MOBILITY_N = 5.0e-2
MOBILITY_P = 2.0e-2
DIFFUSION_N = K_B * T_0 * MOBILITY_N / Q
DIFFUSION_P = K_B * T_0 * MOBILITY_P / Q
THERMAL_EXPANSION = 1.0e-5
PYROELECTRIC = 1.0e-8

COLORS = {5.0: "#2676D2", 10.0: "#D94B40", 15.0: "#2C9A67"}
MARKERS = {5.0: "^", 10.0: "v", 15.0: "D"}
INK = "#161616"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial", "DejaVu Sans"],
        "font.size": 6.6,
        "axes.labelsize": 6.8,
        "xtick.labelsize": 6.2,
        "ytick.labelsize": 6.2,
        "legend.fontsize": 5.2,
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
    }
)


def read_base_fem() -> np.ndarray:
    raw = np.genfromtxt(BASE_RUN / "A1_ne20_nodes.csv", delimiter=",", names=True)
    order = np.argsort(raw["x_m"], kind="stable")
    raw = raw[order]
    rounded_x = np.round(raw["x_m"], 18)
    keep = np.r_[True, np.diff(rounded_x) != 0.0]
    return raw[keep]


def analytic_fields(x: np.ndarray, delta_t: float) -> dict[str, np.ndarray]:
    kappa_bar = PERMITTIVITY + PIEZO_E**2 / ELASTIC_C
    gamma_squared = Q * (
        REFERENCE_N * MOBILITY_N / DIFFUSION_N
        + REFERENCE_P * MOBILITY_P / DIFFUSION_P
    ) / kappa_bar
    gamma = np.sqrt(gamma_squared)
    source = (PYROELECTRIC + PIEZO_E * THERMAL_EXPANSION) * delta_t / (
        kappa_bar * LENGTH
    )
    phi = source / gamma_squared * (
        np.cosh(gamma * (x - 0.5 * LENGTH)) / np.cosh(0.5 * gamma * LENGTH) - 1.0
    )
    eta_n = MOBILITY_N / DIFFUSION_N * phi
    eta_p = -MOBILITY_P / DIFFUSION_P * phi
    displacement = (
        THERMAL_EXPANSION * delta_t * x**2 / (2.0 * LENGTH)
        - PIEZO_E / ELASTIC_C * phi
    )
    return {"ux": displacement, "phi": phi, "eta_n": eta_n, "eta_p": eta_p}


FIELDS = {
    "ux": {
        "base": "u_m",
        "scale": 1.0e12,
        "ylabel": r"轴向位移 $u_x$ (pm)",
    },
    "phi": {
        "base": "phi_V",
        "scale": 1.0e3,
        "ylabel": r"电势 $\varphi$ (mV)",
    },
    "eta_n": {
        "base": "eta_n",
        "scale": 100.0,
        "ylabel": r"电子相对浓度扰动 $\eta_n$ (%)",
    },
    "eta_p": {
        "base": "eta_p",
        "scale": 100.0,
        "ylabel": r"空穴相对浓度扰动 $\eta_p$ (%)",
    },
}


def style_axis(ax: plt.Axes, ylabel: str) -> None:
    ax.set_xlim(0.0, 1.0)
    ax.set_xticks(np.linspace(0.0, 1.0, 6))
    ax.set_xlabel(r"$x_1/L$")
    ax.set_ylabel(ylabel)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(INK)
        spine.set_linewidth(0.8)
    ax.tick_params(colors=INK, pad=2.0)


def paired_legend_handles() -> list[Line2D]:
    handles: list[Line2D] = []
    for delta_t in DELTA_T_VALUES:
        color = COLORS[delta_t]
        marker = MARKERS[delta_t]
        handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                lw=1.15,
                marker="o",
                markersize=2.8,
                markerfacecolor=color,
                markeredgecolor=color,
                label=rf"解析解，$\Delta T={delta_t:g}$ K",
            )
        )
        handles.append(
            Line2D(
                [0],
                [0],
                linestyle="none",
                marker=marker,
                markersize=4.0,
                markerfacecolor="white",
                markeredgecolor=color,
                markeredgewidth=0.9,
                label=rf"FEM，$\Delta T={delta_t:g}$ K",
            )
        )
    return handles


def draw_field(
    ax: plt.Axes,
    field: str,
    x_dense: np.ndarray,
    x_fem: np.ndarray,
    base: np.ndarray,
    show_legend: bool,
) -> None:
    spec = FIELDS[field]
    scale = spec["scale"]
    for delta_t in DELTA_T_VALUES:
        color = COLORS[delta_t]
        marker = MARKERS[delta_t]
        exact = analytic_fields(x_dense, delta_t)[field]
        fem = base[spec["base"]] * (delta_t / BASE_DELTA_T)
        ax.plot(
            x_dense / LENGTH,
            exact * scale,
            color=color,
            lw=1.15,
            marker="o",
            markersize=2.15,
            markerfacecolor=color,
            markeredgecolor=color,
            markeredgewidth=0.0,
            markevery=25,
            zorder=2,
        )
        ax.plot(
            x_fem / LENGTH,
            fem * scale,
            linestyle="none",
            marker=marker,
            markersize=3.6,
            markerfacecolor="white",
            markeredgecolor=color,
            markeredgewidth=0.85,
            markevery=4,
            zorder=3,
        )
    style_axis(ax, spec["ylabel"])
    ax.margins(y=0.08)
    if show_legend:
        legend_locations = {
            "ux": "upper left",
            "phi": "upper center",
            "eta_n": "upper center",
            "eta_p": "lower center",
        }
        ax.legend(
            handles=paired_legend_handles(),
            loc=legend_locations[field],
            frameon=False,
            handlelength=2.4,
            borderpad=0.15,
            labelspacing=0.22,
            handletextpad=0.45,
        )


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=450, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")


def write_source_data(base: np.ndarray) -> None:
    x = base["x_m"]
    rows = []
    for delta_t in DELTA_T_VALUES:
        exact = analytic_fields(x, delta_t)
        ratio = delta_t / BASE_DELTA_T
        for i, x_value in enumerate(x):
            rows.append(
                [
                    delta_t,
                    x_value,
                    x_value / LENGTH,
                    base["u_m"][i] * ratio,
                    exact["ux"][i],
                    base["phi_V"][i] * ratio,
                    exact["phi"][i],
                    base["eta_n"][i] * ratio,
                    exact["eta_n"][i],
                    base["eta_p"][i] * ratio,
                    exact["eta_p"][i],
                ]
            )
    np.savetxt(
        OUT / "A1_three_temperature_source_data.csv",
        np.asarray(rows),
        delimiter=",",
        header=(
            "delta_T_K,x_m,x_over_L,ux_fem_m,ux_analytic_m,phi_fem_V,phi_analytic_V,"
            "eta_n_fem,eta_n_analytic,eta_p_fem,eta_p_analytic"
        ),
        comments="",
    )


def verify_base_formula(base: np.ndarray) -> None:
    exact = analytic_fields(base["x_m"], BASE_DELTA_T)
    checks = {
        "ux": (base["u_exact_m"], exact["ux"]),
        "phi": (base["phi_exact_V"], exact["phi"]),
        "eta_n": (base["eta_n_exact"], exact["eta_n"]),
        "eta_p": (base["eta_p_exact"], exact["eta_p"]),
    }
    for name, (archived, recalculated) in checks.items():
        scale = max(float(np.max(np.abs(archived))), 1.0e-30)
        if not np.allclose(archived, recalculated, rtol=2.0e-11, atol=scale * 2.0e-13):
            raise RuntimeError(f"Independent analytical formula does not match archived A1 field: {name}")


def main() -> None:
    base = read_base_fem()
    verify_base_formula(base)
    write_source_data(base)
    x_fem = base["x_m"]
    x_dense = np.linspace(0.0, LENGTH, 501)

    for field in FIELDS:
        fig, ax = plt.subplots(figsize=(3.40, 2.75))
        draw_field(ax, field, x_dense, x_fem, base, show_legend=True)
        fig.tight_layout(pad=0.45)
        save_figure(fig, OUT / f"A1_three_DT_{field}")
        plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(7.20, 5.35))
    for ax, field in zip(axes.ravel(), FIELDS):
        draw_field(ax, field, x_dense, x_fem, base, show_legend=True)
    fig.subplots_adjust(left=0.095, right=0.985, bottom=0.095, top=0.975, wspace=0.26, hspace=0.30)
    save_figure(fig, OUT / "A1_three_DT_four_fields")
    plt.close(fig)

    max_eta_15 = float(np.max(np.abs(analytic_fields(x_dense, 15.0)["eta_n"])))
    contract = {
        "core_conclusion": "In the frozen linearized A1 model, increasing delta-T from 5 to 15 K scales all four fields linearly without changing their normalized spatial shapes.",
        "archetype": "quantitative grid plus four single-panel exports",
        "temperature_differences_K": list(DELTA_T_VALUES),
        "fields": ["u_x", "phi", "eta_n", "eta_p"],
        "analytic_curves": "Independently evaluated from the frozen A1 closed-form equations on 501 x locations; 21 uniformly spaced filled markers are displayed on each curve.",
        "fem_markers": "10 K archived Ne=20 FEM solution; 5 K and 15 K markers are exact linear load scalings of that solution.",
        "maximum_absolute_eta_at_15K": max_eta_15,
        "small_perturbation_limit": 0.05,
        "source_data": "A1_three_temperature_source_data.csv",
        "exports": ["SVG", "PDF", "PNG 450 dpi", "TIFF 600 dpi"],
    }
    (OUT / "A1_three_temperature_figure_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"A1_THREE_TEMPERATURE_FIGURES_OK={OUT}")


if __name__ == "__main__":
    main()
