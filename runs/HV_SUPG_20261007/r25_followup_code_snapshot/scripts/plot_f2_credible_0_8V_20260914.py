"""Credibility-aware 0--8 V F2 summary.

This figure deliberately separates (i) the algebraic LQ reference, (ii) NT
states that pass the frozen solve/flux gates but remain unresolved by the
high-Peclet diagnostic, (iii) retained exploratory states that fail a gate,
and (iv) voltages for which no complete NT state exists.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import shutil

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "output/F2_孔洞电压扩展_0_10V_20260909/F2_voltage_sweep_0_10V_source.csv"
OUT = ROOT / "output/F2_0_8V_可信性审计_20260914"
OUT.mkdir(parents=True, exist_ok=True)

with SRC.open(encoding="utf-8-sig", newline="") as stream:
    source_rows = [row for row in csv.DictReader(stream)
                   if float(row["voltage_V"]) <= 8.0]
by_voltage = {round(float(row["voltage_V"]), 6): row for row in source_rows}
anchor = by_voltage[3.0]
voltages = np.arange(0.0, 8.01, 0.5)

specs = [
    ("theta_max_K", 1.0, 2, r"最大温升 $\theta_{\max}$ (K)", "max temperature"),
    ("joule_W", 1.0e3, 2, r"总焦耳功率 $P_J$ (mW)", "total Joule power"),
    ("ring_E_mean_V_m", 1.0e-5, 1, r"孔周平均电场 $\langle |E|\rangle$ (kV/cm)", "annulus electric field"),
    ("ring_Q_mean_W_m3", 1.0e-12, 2,
     r"孔周平均焦耳热密度 $\langle Q_J\rangle$ ($10^{12}$ W m$^{-3}$)",
     "annulus Joule density"),
]

def number(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.nan
    return result if math.isfinite(result) else math.nan

def lq_value(key, exponent, voltage):
    return number(anchor["lq_" + key]) * (voltage / 3.0) ** exponent

def nt_gate(row):
    if row is None:
        return "NO_COMPLETE_STATE", "NO_STATE", "NO_STATE", "NO_COMPLETE_NT_STATE"
    complete = row["nt_solve_complete"].strip().lower() == "true"
    status = row["nt_status"].strip()
    if not complete:
        return "NO_COMPLETE_STATE", "NO_STATE", "NO_STATE", "NO_COMPLETE_NT_STATE"
    heat = number(row["nt_heat_balance_relative"])
    quad = number(row["nt_joule_quadrature_4_5_relative"])
    solve_gate = "PASS" if status == "PASS_NOT_MESH_ACCEPTANCE" else "FAIL"
    heat_gate = "PASS" if (math.isfinite(heat) and heat < 1.0e-2 and
                           math.isfinite(quad) and quad < 1.0e-3) else "FAIL"
    voltage = number(row["voltage_V"])
    if solve_gate == "PASS" and heat_gate == "PASS":
        if voltage >= 3.5:
            cls = "SOLVE_BALANCE_PASS_HIGH_PE_UNRESOLVED"
        else:
            cls = "SOLVE_BALANCE_PASS_PE_NOT_AUDITED"
    else:
        cls = "EXPLORATORY_STATE_GATE_FAIL"
    return cls, solve_gate, heat_gate, cls

# A compact credibility table is part of the deliverable, not an inference
# from the plotted curves.
credibility_rows = []
for voltage in voltages:
    key_voltage = round(float(voltage), 6)
    row = by_voltage.get(key_voltage)
    if abs(voltage) < 1.0e-12:
        nt_class, solve_gate, heat_gate, pe_gate = (
            "EXACT_ZERO_STATE", "PASS", "PASS", "NOT_APPLICABLE")
    else:
        nt_class, solve_gate, heat_gate, _ = nt_gate(row)
        if nt_class == "SOLVE_BALANCE_PASS_HIGH_PE_UNRESOLVED":
            pe_gate = "FAIL_HIGH_PE_DIAGNOSTIC"
        elif nt_class == "SOLVE_BALANCE_PASS_PE_NOT_AUDITED":
            pe_gate = "NOT_AUDITED"
        elif nt_class == "EXPLORATORY_STATE_GATE_FAIL":
            pe_gate = "NOT_CLEARED"
        else:
            pe_gate = "NO_STATE"
    nt_value = number(row["nt_theta_max_K"]) if row is not None and row["nt_theta_max_K"].strip() else math.nan
    credibility_rows.append({
        "voltage_V": f"{voltage:.1f}",
        "lq_theta_max_K": f"{lq_value('theta_max_K', 2, voltage):.16e}",
        "nt_theta_max_K": "" if math.isnan(nt_value) else f"{nt_value:.16e}",
        "nt_source": "" if row is None else row["nt_source"],
        "nt_status": "" if row is None else row["nt_status"],
        "solve_gate": solve_gate,
        "heat_flux_gate": heat_gate,
        "peclet_gate": pe_gate,
        "credibility_class": nt_class,
        "paper_quantitative_use": "NO" if nt_class != "EXACT_ZERO_STATE" else "YES_ZERO_CHECK",
    })

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "Arial", "DejaVu Sans"],
    "font.size": 10.5,
    "axes.labelsize": 11.5,
    "legend.fontsize": 8.5,
    "axes.linewidth": 1.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "axes.unicode_minus": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.0))
fig.subplots_adjust(left=0.095, right=0.975, bottom=0.11, top=0.965,
                    wspace=0.30, hspace=0.36)
plot_records = []

for ax, (key, scale, exponent, ylabel, sampling) in zip(axes.flat, specs):
    # Per the frozen figure convention, the linearized reference is shown as
    # points (algebraic scaling at the requested voltage values), while the
    # available nonlinear states are connected by a line.  This avoids
    # visually implying that the LQ curve consists of independent 8 V solves.
    lq_points = np.array([lq_value(key, exponent, value) * scale for value in voltages])
    ax.scatter(voltages, lq_points, s=24, color="#2676D2",
               label="线性化 LQ 点（3 V 场精确缩放）", zorder=3)

    accepted_rows = [r for r in source_rows
                     if r["nt_solve_complete"].strip().lower() == "true"
                     and r["nt_status"].strip() == "PASS_NOT_MESH_ACCEPTANCE"
                     and number(r["nt_heat_balance_relative"]) < 1.0e-2
                     and number(r["nt_joule_quadrature_4_5_relative"]) < 1.0e-3]
    accepted_x = np.array([number(r["voltage_V"]) for r in accepted_rows])
    accepted_y = np.array([number(r["nt_" + key]) * scale for r in accepted_rows])
    accepted_mask = np.isfinite(accepted_x) & np.isfinite(accepted_y)
    ax.plot(accepted_x[accepted_mask], accepted_y[accepted_mask], "-",
            lw=1.6, color="#D94B40",
            label="NT：守恒通过，高 Pe 未解除", zorder=2)
    ax.scatter(accepted_x[accepted_mask], accepted_y[accepted_mask],
               s=24, facecolors="white", edgecolors="#D94B40", linewidths=1.1,
               zorder=4)

    failed_rows = [r for r in source_rows
                   if r["nt_solve_complete"].strip().lower() == "true"
                   and (r["nt_status"].strip() != "PASS_NOT_MESH_ACCEPTANCE"
                        or number(r["nt_heat_balance_relative"]) >= 1.0e-2
                        or number(r["nt_joule_quadrature_4_5_relative"]) >= 1.0e-3)]
    failed_x = np.array([number(r["voltage_V"]) for r in failed_rows])
    failed_y = np.array([number(r["nt_" + key]) * scale for r in failed_rows])
    failed_mask = np.isfinite(failed_x) & np.isfinite(failed_y)
    ax.scatter(failed_x[failed_mask], failed_y[failed_mask], marker="x", s=34,
               color="#D98C27", linewidths=1.3,
               label="NT：保留状态但门禁失败", zorder=4)

    ax.set(xlim=(0, 8), xlabel=r"电压幅值 $V$ (V)", ylabel=ylabel)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    upper = max(float(np.nanmax(lq_points)),
                float(np.nanmax(failed_y)) if failed_y.size else 0.0,
                float(np.nanmax(accepted_y)) if accepted_y.size else 0.0)
    ax.set_ylim(-0.04 * upper, 1.12 * upper)
    ax.tick_params(length=4, width=1)
    ax.legend(loc="upper left", frameon=False, handlelength=2.3)
    ax.text(0.99, 0.02, "NT 无完整状态：5.5、6.5、7.0、7.5、8.0 V",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
            color="#666666")

    for voltage in voltages:
        plot_records.append({
            "model": "LQ", "voltage_V": f"{voltage:.1f}", "metric": key,
            "value_SI": f"{lq_value(key, exponent, voltage):.16e}",
            "display_value": f"{lq_value(key, exponent, voltage) * scale:.16e}",
            "status": "LINEAR_REFERENCE_BASELINE",
            "sampling": sampling,
        })
    for r in source_rows:
        value = number(r["nt_" + key])
        if math.isfinite(value):
            plot_records.append({
                "model": "NT", "voltage_V": r["voltage_V"], "metric": key,
                "value_SI": f"{value:.16e}",
                "display_value": f"{value * scale:.16e}",
                "status": nt_gate(r)[0], "sampling": sampling,
            })

stem = OUT / "F2_0_8V_可信性分层_四指标"
for extension in ("png", "pdf", "svg", "tiff"):
    options = {"pil_kwargs": {"compression": "tiff_lzw"}} if extension == "tiff" else {}
    fig.savefig(stem.with_suffix("." + extension),
                dpi=600 if extension == "tiff" else 450,
                bbox_inches="tight", **options)
plt.close(fig)

with (OUT / "plot_data_credible.csv").open("w", encoding="utf-8-sig", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(plot_records[0]))
    writer.writeheader()
    writer.writerows(plot_records)

with (OUT / "credibility_status_0_8V.csv").open("w", encoding="utf-8-sig", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(credibility_rows[0]))
    writer.writeheader()
    writer.writerows(credibility_rows)

contract = {
    "source": str(SRC.relative_to(ROOT)),
    "source_sha256": hashlib.sha256(SRC.read_bytes()).hexdigest(),
    "geometry_um": [20, 20, 1], "hole_radius_um": 3,
    "voltage_axis_V": [0, 8], "voltage_step_V": 0.5,
    "LQ": "Blue points at 0.5 V increments, exact algebraic scaling of the accepted 3 V LQ FEM state; not independent solves.",
    "NT_solve_balance_pass": "Open red circles connected by a solid red segment; high-Pe transport remains unresolved.",
    "NT_gate_fail": "Orange x markers; retained exploratory states are not quantitative results.",
    "NT_no_complete_state_V": [5.5, 6.5, 7.0, 7.5, 8.0],
    "credibility_rule": "No NT point is paper-quantitative until an exponential-fitting/stabilized transport audit passes.",
    "annulus_um": {"radial_range": [3, 3.5], "z_range": [-0.5, 0.5]},
    "sampling": "5x5x5 global field samples for theta_max; 5x5x3 physical-radius weighted annulus for E and Q.",
}
(OUT / "figure_contract.json").write_text(json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8")
(OUT / "source_status_0_8V.csv").write_text(SRC.read_text(encoding="utf-8-sig"), encoding="utf-8-sig")
shutil.copyfile(__file__, OUT / Path(__file__).name)

print(f"CREDIBLE_F2_0_8V_DONE={OUT}")
