"""Build source-backed H1 LQ-versus-NT voltage-sweep figures.

LQ values are exact voltage scalings of the accepted 3 V H1-LQ FEM solution.
NT values are independently solved FEM states at every nonzero voltage; 0 V is
the exact zero state.  High-voltage states that converge but fail only the frozen
independent heat-flux audit are retained and visibly flagged.  No nonlinear
values are interpolated before the line is drawn.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
F1_RUNS = ROOT / "runs" / "F1_FAST_HOLE_L0_20260907"
F2_RUNS = ROOT / "runs" / "F2_H1_VOLTAGE_SWEEP_L0_20260908"
INPUTS = ROOT / "cases" / "F2_H1_VOLTAGE_SWEEP_20260908"
OUT = ROOT.parent / "output" / "F2_孔洞电压扫描_20260908"
TOOLS = F1_RUNS / "evidence" / "delivery_scripts"
sys.path.insert(0, str(TOOLS))
import plot_f1_main_results as fe  # noqa: E402
from audit_hex20_boundary_geometry import read_mesh  # noqa: E402


VOLTAGES = np.arange(0.0, 5.0 + 0.25, 0.5)
LQ_ANCHOR_V = 3.0
LQ_CASE = F1_RUNS / "FAST_H1_LQ_L0"
NT_3V_CASE = F1_RUNS / "FAST_H1_NT_L0"
METRIC_KEYS = (
    "joule_W", "terminal_power_W", "theta_max_K", "electron_right_A",
    "ring_E_mean_V_m", "ring_Q_mean_W_m3", "ring_Q_p95_W_m3",
    "ring_theta_mean_K", "maximum_abs_n_relative", "minimum_mu_n_ratio",
    "heat_balance_relative", "electron_balance_relative",
    "joule_quadrature_4_5_relative", "maximum_accepted_residual",
    "maximum_accepted_update", "maximum_stage_balance",
)


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "Arial", "DejaVu Sans"],
    "font.size": 8.2,
    "axes.titlesize": 9.0,
    "axes.labelsize": 8.6,
    "legend.fontsize": 7.6,
    "axes.linewidth": 0.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": False,
    "ytick.right": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def case_tag(voltage: float) -> str:
    return f"V{round(voltage * 10):03d}"


def read_metrics(path: Path, *, allow_heat_flux_fail: bool = False) -> dict[str, float | str]:
    result: dict[str, float | str] = {}
    for line in (path / "metrics.txt").read_text(encoding="utf-8").splitlines():
        if "=" not in line or line.startswith(("continuation_columns", "balance_detail", "stage")):
            continue
        key, value = line.split("=", 1)
        value = value.strip().rstrip(",")
        try:
            result[key.strip()] = float(value)
        except ValueError:
            result[key.strip()] = value
    acceptance = result.get("case_acceptance")
    if acceptance == "PASS_NOT_MESH_ACCEPTANCE":
        return result
    only_heat_flux_failed = (
        allow_heat_flux_fail and acceptance == "FAIL"
        and "f1_solve=COMPLETE" in (path / "stdout.txt").read_text(encoding="utf-8")
        and (path / "nodes.csv").is_file()
        and float(result["maximum_stage_balance"]) < 1.0e-8
        and float(result["maximum_accepted_residual"]) < 1.0e-9
        and float(result["maximum_accepted_update"]) < 1.0e-9
        and float(result["electron_balance_relative"]) < 0.01
        and float(result["joule_quadrature_4_5_relative"]) < 0.001
        and float(result["heat_balance_relative"]) >= 0.01
    )
    if not only_heat_flux_failed:
        raise ValueError(f"Unaccepted or incomplete case: {path}")
    return result


def exact_zero_nt() -> dict[str, float | str]:
    result: dict[str, float | str] = {key: 0.0 for key in METRIC_KEYS}
    result.update({
        "minimum_mu_n_ratio": 1.0,
        "case_acceptance": "EXACT_ZERO_STATE",
        "source": "exact_zero_state",
    })
    return result


def nt_case(voltage: float) -> Path | None:
    if voltage == 0:
        return None
    if voltage == 3.0:
        return NT_3V_CASE
    return F2_RUNS / f"F2_H1_NT_L0_{case_tag(voltage)}"


def lq_scaled(anchor: dict[str, float | str], voltage: float) -> dict[str, float | str]:
    a = voltage / LQ_ANCHOR_V
    out = dict(anchor)
    for key in ("joule_W", "terminal_power_W", "theta_max_K", "ring_Q_mean_W_m3",
                "ring_Q_p95_W_m3", "ring_theta_mean_K"):
        out[key] = float(anchor[key]) * a * a
    for key in ("electron_right_A", "ring_E_mean_V_m"):
        out[key] = float(anchor[key]) * a
    out["maximum_abs_n_relative"] = 0.0
    out["minimum_mu_n_ratio"] = 1.0
    out["voltage_V"] = voltage
    out["source"] = "exact_scaling_of_accepted_H1_LQ_3V_FEM"
    return out


def save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.png", dpi=450, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def node_state(path: Path) -> np.ndarray:
    data = np.genfromtxt(path / "nodes.csv", delimiter=",", names=True, usecols=range(13))
    return np.column_stack([data[c] for c in (
        "ux_m", "uy_m", "uz_m", "phi_V", "eta_n", "eta_p", "theta_K"
    )])


def add_hole(ax: plt.Axes) -> None:
    ax.add_patch(Circle((0, 0), 3.0, fill=False, ec="#222222", lw=0.65, zorder=20))


def load_rows() -> tuple[list[dict], list[dict]]:
    anchor = read_metrics(LQ_CASE)
    lq_rows, nt_rows = [], []
    for voltage in VOLTAGES:
        lq = lq_scaled(anchor, float(voltage))
        nt_path = nt_case(float(voltage))
        nt = exact_zero_nt() if nt_path is None else read_metrics(nt_path, allow_heat_flux_fail=True)
        nt["voltage_V"] = float(voltage)
        nt["source"] = "exact_zero_state" if nt_path is None else str(nt_path.relative_to(ROOT))
        lq_rows.append(lq)
        nt_rows.append(nt)
    return lq_rows, nt_rows


def export_csv(lq_rows: list[dict], nt_rows: list[dict]) -> list[dict]:
    records = []
    for voltage, lq, nt in zip(VOLTAGES, lq_rows, nt_rows):
        row = {"voltage_V": voltage, "lq_source": lq["source"], "nt_source": nt["source"],
               "nt_acceptance": nt["case_acceptance"]}
        for key in METRIC_KEYS:
            row[f"lq_{key}"] = lq.get(key, "")
            row[f"nt_{key}"] = nt.get(key, "")
            lv, nv = lq.get(key), nt.get(key)
            row[f"relative_gap_{key}"] = "" if not isinstance(lv, float) or lv == 0 else (float(nv)-lv)/abs(lv)
        records.append(row)
    with (OUT / "F2_voltage_sweep_source.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)
    return records


def plot_load_response(lq_rows: list[dict], nt_rows: list[dict]) -> None:
    specs = [
        ("theta_max_K", 1.0, "最大温升 (K)"),
        ("joule_W", 1e3, "总焦耳功率 (mW)"),
        ("electron_right_A", 1e3, "右电极电子电流 (mA)"),
        ("ring_E_mean_V_m", 1e-5, "孔周平均电场 (kV/cm)"),
        ("ring_Q_mean_W_m3", 1e-12, r"孔周平均焦耳热 ($10^{12}$ W m$^{-3}$)"),
        ("minimum_mu_n_ratio", 1.0, r"最小迁移率比 $\mu_n/\mu_{n0}$"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(7.20, 4.55), layout="constrained")
    for ax, (key, scale, ylabel) in zip(axes.flat, specs):
        lq = np.array([float(r[key]) for r in lq_rows]) * scale
        nt = np.array([float(r[key]) for r in nt_rows]) * scale
        ax.scatter(VOLTAGES, lq, s=25, marker="o", color="#2676D2", zorder=4,
                   label="线性 LQ（3 V FEM精确缩放）")
        first_fail = next((i for i, r in enumerate(nt_rows)
                           if r["case_acceptance"] == "FAIL"), len(nt_rows))
        ax.plot(VOLTAGES[:first_fail], nt[:first_fail], color="#D94B40", lw=1.75,
                label="非线性 NT（通过验收）")
        if first_fail < len(nt_rows):
            start = max(0, first_fail-1)
            ax.plot(VOLTAGES[start:], nt[start:], color="#D94B40", lw=1.75, ls="--",
                    label="NT（收敛；热流审计未通过）")
        ax.set(xlabel="电压幅值 (V)", ylabel=ylabel, xlim=(-0.08, 5.08), xticks=np.arange(0, 6))
        ax.ticklabel_format(axis="y", useOffset=False)
        ax.margins(y=0.12)
    axes[0, 0].legend(loc="upper left")
    fig.suptitle("中心圆孔ZnO薄板：线性与非线性电压响应", fontsize=10.2)
    save(fig, "F2_线性点_非线性线_六指标")


def plot_gap(records: list[dict]) -> None:
    active = records[1:]
    x = np.array([float(r["lq_theta_max_K"]) for r in active])
    specs = [
        ("theta_max_K", "最大温升差"),
        ("joule_W", "总焦耳功率差"),
        ("ring_E_mean_V_m", "孔周平均电场差"),
        ("ring_Q_mean_W_m3", "孔周平均焦耳热差"),
    ]
    fig, ax = plt.subplots(figsize=(4.70, 3.35), layout="constrained")
    colors = ["#D94B40", "#6A51A3", "#238B45", "#CB6D2D"]
    for (key, label), color in zip(specs, colors):
        gap = 100*np.array([float(r[f"relative_gap_{key}"]) for r in active])
        first_fail = next((i for i, r in enumerate(active)
                           if r["nt_acceptance"] == "FAIL"), len(active))
        ax.plot(x[:first_fail], gap[:first_fail], lw=1.55, color=color, label=label)
        if first_fail < len(active):
            start = max(0, first_fail-1)
            ax.plot(x[start:], gap[start:], lw=1.55, ls="--", color=color)
    ax.axhline(0, color="#444444", lw=0.65)
    ax.set(xlabel="线性LQ预测的最大温升 (K)", ylabel="NT相对LQ偏差 (%)")
    ax.legend(loc="best")
    ax.text(0.98, 0.04, "虚线段：L0独立热流审计未通过", transform=ax.transAxes,
            ha="right", va="bottom", color="#555555", fontsize=7.2)
    ax.margins(x=.03, y=.12)
    save(fig, "F2_非线性偏离随温升")


def plot_clouds_5v() -> None:
    xyz, conn, _ = read_mesh(INPUTS / "H1_L0.mesh")
    lq3 = node_state(LQ_CASE)
    ratio = 5.0 / 3.0
    lq5 = lq3.copy()
    lq5[:, 0:3] *= ratio  # not displayed; retained only for a complete state array
    lq5[:, 3] *= ratio
    lq5[:, 4:6] = 0.0
    lq5[:, 6] *= ratio**2
    nt5_case = F2_RUNS / "F2_H1_NT_L0_V050"
    nt5_metrics = read_metrics(nt5_case, allow_heat_flux_fail=True)
    nt5 = node_state(nt5_case)
    lq = fe.plane(xyz, conn, lq5, "LQ", 0.0)
    nt = fe.plane(xyz, conn, nt5, "NT", 0.0)
    fields = [
        ("phi", 1.0, "电势", "V", "viridis"),
        ("E", 1e-5, "电场强度", "kV/cm", "viridis"),
        ("Q", 1e-12, "焦耳热", r"$10^{12}$ W m$^{-3}$", "magma"),
        ("theta", 1.0, "温升", "K", "magma"),
    ]
    fig, axes = plt.subplots(4, 2, figsize=(6.70, 10.1), layout="constrained")
    for row, (key, scale, name, unit, cmap) in enumerate(fields):
        arrays = [lq[key]*scale, nt[key]*scale]
        bounds = (min(a.min() for a in arrays), max(a.max() for a in arrays))
        for col, (model, data, values) in enumerate(zip(("线性 LQ", "非线性 NT"), (lq, nt), arrays)):
            suffix = "" if model == "线性 LQ" else "（L0热流审计未通过）"
            fe.panel(fig, axes[row, col], data, values,
                     f"{model}：{name}，5 V，z=0{suffix}", unit, cmap, bounds)
            add_hole(axes[row, col])
    fig.suptitle("中心圆孔模型5 V场分布对照（各行共用色标范围）", fontsize=10.2)
    save(fig, "F2_5V_线性非线性云图对照")

    diffs = [
        (nt["E"]-lq["E"], 1e-5, "电场差 NT−LQ", "kV/cm"),
        (nt["Q"]-lq["Q"], 1e-12, "焦耳热差 NT−LQ", r"$10^{12}$ W m$^{-3}$"),
        (nt["theta"]-lq["theta"], 1.0, "温升差 NT−LQ", "K"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(7.20, 2.65), layout="constrained")
    for ax, (values, scale, title, unit) in zip(axes, diffs):
        scaled = values*scale
        limit = max(abs(float(scaled.min())), abs(float(scaled.max())))
        fe.panel(fig, ax, nt, scaled, title, unit, "RdBu_r", (-limit, limit))
        add_hole(ax)
    fig.suptitle("5 V非线性修正的空间位置（NT：L0独立热流审计未通过）", fontsize=10.2)
    save(fig, "F2_5V_NT减LQ差值云图")
    np.savez_compressed(OUT / "F2_5V_midplane_cloud_source.npz",
                        xyz=lq["xyz"], triangles=lq["triangles"],
                        lq_phi=lq["phi"], nt_phi=nt["phi"],
                        lq_E=lq["E"], nt_E=nt["E"],
                        lq_Q=lq["Q"], nt_Q=nt["Q"],
                        lq_theta=lq["theta"], nt_theta=nt["theta"])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=False)
    lq_rows, nt_rows = load_rows()
    records = export_csv(lq_rows, nt_rows)
    plot_load_response(lq_rows, nt_rows)
    plot_gap(records)
    plot_clouds_5v()
    max_nt = nt_rows[-1]
    contract = {
        "conclusion": "The H1 nonlinear solution departs progressively from the frozen LQ voltage scaling as self-heating grows; high-voltage L0 points are exploratory where the independent heat-flux audit fails.",
        "geometry": "20 um x 20 um x 1 um ZnO plate, centered insulated cylindrical hole r=3 um, c axis +z.",
        "boundary_conditions": "x=-10/+10 um: phi=0/-V, eta_n=eta_p=0, theta=0; remaining external and hole faces natural; six minimal rigid-body constraints.",
        "line_encoding": "LQ blue unconnected points; NT red solid line for accepted points and red dashed line for converged high-voltage points that fail only the independent heat-flux audit. 0 V is the exact zero state.",
        "lq_provenance": "Exact scaling of accepted H1-LQ 3 V FEM fields: phi,E,J proportional to V; Q and theta proportional to V^2. No scaled point is described as an independent solve.",
        "nt_provenance": "Accepted H1-NT L0 solves at 0.5 V increments; the accepted F1 3 V case is reused without modification.",
        "acceptance_scope": "Fast L0 only. PASS means solve/flux acceptance; FAIL high-voltage points are included only when the nonlinear solve and every frozen criterion except the 1% independent heat-flux balance pass. No point is mesh-independent or experimentally validated.",
        "nt_acceptance_by_voltage": {str(float(v)): r["case_acceptance"] for v, r in zip(VOLTAGES, nt_rows)},
        "five_volt_summary": {key: max_nt.get(key) for key in METRIC_KEYS},
    }
    (OUT / "figure_contract.json").write_text(json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"F2_FIGURE_OUTPUT={OUT}")


if __name__ == "__main__":
    main()
