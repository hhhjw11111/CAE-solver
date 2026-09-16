"""Plot A1 analytical comparisons in the user's Program/Analytic reference style."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "runs" / "A1_100NM_010NM_NE20_DT10K_20260907"
OUT = ROOT.parent / "output" / "figures" / "A1_100NM_010NM_NE20_DT10K_20260907"
OUT.mkdir(parents=True, exist_ok=True)

PROGRAM = "#111111"
ANALYTIC = "#E21A1A"
GRID = "#D9D9D9"
TEXT = "#111111"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial", "DejaVu Sans"],
        "font.size": 7.0,
        "axes.titlesize": 7.8,
        "axes.labelsize": 7.0,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.7,
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 1.0,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.major.size": 5.0,
        "ytick.major.size": 5.0,
    }
)


def read_nodes():
    raw = np.genfromtxt(DATA / "A1_ne20_nodes.csv", delimiter=",", names=True)
    # The exact solution is one-dimensional. Select one representative line through y=z=0.
    order = np.argsort(raw["x_m"])
    raw = raw[order]
    keep = np.r_[True, np.diff(np.round(raw["x_m"], 18)) != 0.0]
    return raw[keep]


def read_metrics():
    values = {}
    for line in (DATA / "A1_stdout.txt").read_text(encoding="utf-8-sig").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        try:
            values[key] = float(value)
        except ValueError:
            pass
    return values


def style_axis(ax, title, ylabel):
    ax.set_title(title, loc="left", fontweight="bold", color=TEXT, pad=7.0)
    ax.set_xlabel("x/L", color=TEXT)
    ax.set_ylabel(ylabel, color=TEXT)
    ax.set_xlim(0.0, 1.0)
    ax.set_xticks(np.linspace(0.0, 1.0, 6))
    ax.grid(axis="y", color=GRID, linewidth=0.8, alpha=0.85)
    ax.set_axisbelow(True)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_color(PROGRAM)
        ax.spines[spine].set_linewidth(1.0)
    ax.tick_params(colors=TEXT)


def style_convergence_axis(ax):
    ax.set_title("A1 convergence", loc="left", fontweight="bold", color=TEXT, pad=7.0)
    ax.set_xlabel("Elements along x", color=TEXT)
    ax.set_ylabel("Relative error", color=TEXT)
    ax.grid(axis="y", color=GRID, linewidth=0.8, alpha=0.85, which="both")
    ax.set_axisbelow(True)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_color(PROGRAM)
        ax.spines[spine].set_linewidth(1.0)
    ax.tick_params(colors=TEXT)


def plot_field(ax, x, program, analytic, title, ylabel, scale=1.0, legend=False):
    ax.plot(x, program * scale, color=PROGRAM, linewidth=1.35, label="Program", zorder=2)
    ax.plot(x, analytic * scale, color=ANALYTIC, linewidth=1.20, marker="o",
            markersize=3.2, markerfacecolor="white", markeredgecolor=ANALYTIC,
            markeredgewidth=0.95, markevery=max(1, len(x) // 17), label="Analytic", zorder=3)
    style_axis(ax, title, ylabel)
    if legend:
        ax.legend(loc="upper left", frameon=False, handlelength=2.2, borderpad=0.2,
                  labelspacing=0.35, handletextpad=0.55)


def save_figure(fig, path):
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".png"), dpi=450, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")


def main():
    d = read_nodes()
    x = d["x_m"] / d["x_m"][-1]
    entries = [
        ("u_m", "u_exact_m", "x/L - u_x", "Axial displacement, u_x (pm)", 1.0e12, "u_x"),
        ("phi_V", "phi_exact_V", "x/L - phi", "Potential, phi (mV)", 1.0e3, "phi"),
        ("eta_n", "eta_n_exact", "x/L - eta_n", "Electron perturbation, eta_n (%)", 100.0, "eta_n"),
        ("eta_p", "eta_p_exact", "x/L - eta_p", "Hole perturbation, eta_p (%)", 100.0, "eta_p"),
        ("theta_K", "theta_exact_K", "x/L - theta", "Temperature rise, theta (K)", 1.0, "theta"),
    ]

    axis_columns = [d["x_m"] * 1.0e9, x]
    axis_header = ["x_nm", "x_over_L"]
    for program_key, exact_key, _, _, _, stem in entries:
        axis_columns.extend([d[program_key], d[exact_key], d[program_key] - d[exact_key]])
        axis_header.extend([f"{stem}_program", f"{stem}_analytic", f"{stem}_error"])
    np.savetxt(
        DATA / "A1_ne20_axis_comparison.csv",
        np.column_stack(axis_columns),
        delimiter=",",
        header=",".join(axis_header),
        comments="",
    )

    # Individual panels: closest to the user's reference image.
    for program_key, exact_key, title, ylabel, scale, stem in entries:
        fig, ax = plt.subplots(figsize=(5.6, 4.0))
        plot_field(ax, x, d[program_key], d[exact_key], title, ylabel, scale, legend=True)
        fig.tight_layout(pad=0.6)
        save_figure(fig, OUT / f"A1_{stem}_program_vs_analytic")
        plt.close(fig)

    # Composite figure for the manuscript: five fields plus convergence.
    fig, axes = plt.subplots(2, 3, figsize=(7.20, 4.90), constrained_layout=True)
    axes = axes.ravel()
    for i, (program_key, exact_key, title, ylabel, scale, _) in enumerate(entries):
        plot_field(axes[i], x, d[program_key], d[exact_key], title, ylabel, scale, legend=(i == 0))

    metrics = read_metrics()
    ne = np.array([1, 2, 4, 8, 16], dtype=float)
    l2 = np.array([metrics[f"a1_ne{int(n)}_phi_l2"] for n in ne])
    h1 = np.array([metrics[f"a1_ne{int(n)}_phi_h1"] for n in ne])
    ax = axes[5]
    ax.loglog(ne, l2, color=PROGRAM, marker="o", markersize=4.0, linewidth=1.35, label="L2 error")
    ax.loglog(ne, h1, color=ANALYTIC, marker="o", markerfacecolor="white",
              markeredgecolor=ANALYTIC, markersize=4.0, linewidth=1.2, label="H1 seminorm")
    style_convergence_axis(ax)
    ax.set_yscale("log")
    ax.set_xticks(ne)
    ax.set_xticklabels([str(int(n)) for n in ne])
    ax.legend(loc="upper right", frameon=False, handlelength=2.2, borderpad=0.2, labelspacing=0.35)

    convergence_fig, convergence_ax = plt.subplots(figsize=(5.6, 4.0))
    convergence_ax.loglog(ne, l2, color=PROGRAM, marker="o", markersize=5.0,
                          linewidth=1.35, label="L2 error")
    convergence_ax.loglog(ne, h1, color=ANALYTIC, marker="o", markerfacecolor="white",
                          markeredgecolor=ANALYTIC, markersize=5.0, linewidth=1.2,
                          label="H1 seminorm")
    style_convergence_axis(convergence_ax)
    convergence_ax.set_xticks(ne)
    convergence_ax.set_xticklabels([str(int(n)) for n in ne])
    convergence_ax.legend(loc="upper right", frameon=False, handlelength=2.2,
                          borderpad=0.2, labelspacing=0.35)
    convergence_fig.tight_layout(pad=0.6)
    save_figure(convergence_fig, OUT / "A1_error_convergence")
    plt.close(convergence_fig)

    fig.suptitle("A1  100 nm × 10 nm × 10 nm，ΔT = 10 K：程序解—解析解验证",
                 x=0.01, ha="left", fontsize=9.5, fontweight="bold", color=TEXT)
    save_figure(fig, OUT / "A1_linearized_multifield_program_vs_analytic")
    plt.close(fig)

    contract = {
        "figure_name": "A1 一维热—压电—载流子多场耦合线性化解析解验证",
        "style_reference": "Program black solid line; Analytic red line with open-circle markers; white background; light-gray horizontal grid.",
        "source_data": "A1 100 nm × 10 nm × 10 nm, Ne=20, delta-T=10 K nodal output and verification log",
        "panels": ["u_x", "phi", "eta_n", "eta_p", "theta", "L2/H1 convergence"],
        "units": {"u_x": "pm", "phi": "mV", "eta_n": "%", "eta_p": "%", "theta": "K"},
        "statistics": "Deterministic numerical-versus-analytic verification; no replicates or inferential statistics.",
        "editability": "SVG text remains text; PDF uses TrueType text; curves, markers and axes are vector objects.",
        "interpretation": "Agreement validates the linearized coupled field implementation and HEX20 convergence, not the full temperature-dependent nonlinear ZnO production model.",
    }
    (OUT / "A1_reference_style_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("A1_REFERENCE_STYLE_FIGURES_OK")


if __name__ == "__main__":
    main()
