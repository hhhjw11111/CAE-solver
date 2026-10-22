"""Postprocess the requested 6--10 V extension without inventing failed states.

The linear LQ reference is evaluated at every requested voltage by exact
scaling of the accepted 3 V FEM anchor. An NT value is plotted only when the
solver wrote a complete state. Failed high-voltage continuation attempts are
reported separately by their last accepted Joule-load factor.
"""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import postprocess_f2_voltage_sweep as base  # noqa: E402


LOW_NT_CASE = base.nt_case
HIGH_RUNS = ROOT / "runs" / "F2_H1_HIGH_VOLTAGE_L0_20260909"
HIGH_INPUTS = ROOT / "cases" / "F2_H1_HIGH_VOLTAGE_20260909"
OUT = ROOT.parent / "output" / "F2_孔洞电压扩展_0_10V_20260909"
VOLTAGES = np.r_[np.arange(0.0, 5.0 + 0.25, 0.5), np.arange(6.0, 10.0 + 0.5, 1.0)]
EXTRA_KEYS = ("E_above_1e6_volume_fraction", "maximum_abs_p_relative",
              "displacement_max_m", "stress_von_mises_max_Pa")
ATTEMPT_RE = re.compile(
    r"newton_iteration=(\d+), voltage_factor=\s*([0-9.E+\-]+), "
    r"joule_factor=\s*([0-9.E+\-]+), residual=\s*([0-9.E+\-]+)"
)


def high_case_tag(voltage: float) -> str:
    return f"V{round(voltage * 10):03d}"


def nt_case(voltage: float) -> Path | None:
    if voltage <= 5.0:
        return LOW_NT_CASE(voltage)
    return HIGH_RUNS / f"F2_H1_NT_L0_{high_case_tag(voltage)}"


def continuation_evidence(path: Path) -> tuple[float, float, str]:
    stdout = (path / "stdout.txt").read_text(encoding="utf-8", errors="replace")
    groups: list[tuple[float, list[float]]] = []
    for line in stdout.splitlines():
        match = ATTEMPT_RE.search(line)
        if not match:
            continue
        iteration, voltage_factor, joule_factor, residual = match.groups()
        if float(voltage_factor) < 0.999999 or float(joule_factor) <= 0.0:
            continue
        # Solver stdout prints only four decimals; recover the frozen 1/64
        # continuation grid before reporting the last accepted factor.
        factor = round(float(joule_factor) * 64.0) / 64.0
        if int(iteration) == 0 or not groups:
            groups.append((factor, [float(residual)]))
        else:
            groups[-1][1].append(float(residual))
    accepted: list[float] = []
    for current, following in zip(groups, groups[1:]):
        if following[0] > current[0] + 1.0e-12:
            accepted.append(current[0])
    last_accepted = max(accepted, default=0.0)
    last_residual = groups[-1][1][-1] if groups else np.nan
    stderr = (path / "stderr.txt").read_text(encoding="utf-8", errors="replace").strip()
    message = next((line.split("=", 1)[1].strip() for line in stderr.splitlines()
                    if line.startswith("status_message=")), stderr)
    return last_accepted, last_residual, message


def failed_nt(voltage: float, path: Path) -> dict:
    last_factor, last_residual, message = continuation_evidence(path)
    row = {key: np.nan for key in base.METRIC_KEYS + EXTRA_KEYS}
    row.update({
        "voltage_V": voltage,
        "source": str(path.relative_to(ROOT)),
        "case_acceptance": "NEWTON_FAILED_NO_STATE",
        "solve_complete": False,
        "last_accepted_joule_factor": last_factor,
        "last_attempt_residual": last_residual,
        "failure_message": message,
    })
    return row


def load_rows() -> tuple[list[dict], list[dict]]:
    anchor = base.read_metrics(base.LQ_CASE)
    lq_rows, nt_rows = [], []
    for voltage in VOLTAGES:
        voltage = float(voltage)
        lq = base.lq_scaled(anchor, voltage)
        path = nt_case(voltage)
        if path is None:
            nt = base.exact_zero_nt()
            nt.update({"solve_complete": True, "last_accepted_joule_factor": 1.0,
                       "failure_message": ""})
        else:
            stdout = (path / "stdout.txt").read_text(encoding="utf-8", errors="replace")
            if "f1_solve=COMPLETE" in stdout and (path / "nodes.csv").is_file():
                nt = base.read_metrics(path, allow_heat_flux_fail=True)
                nt.update({"solve_complete": True, "last_accepted_joule_factor": 1.0,
                           "failure_message": ""})
            else:
                nt = failed_nt(voltage, path)
        nt["voltage_V"] = voltage
        nt["source"] = "exact_zero_state" if path is None else str(path.relative_to(ROOT))
        lq_rows.append(lq)
        nt_rows.append(nt)
    return lq_rows, nt_rows


def export_csv(lq_rows: list[dict], nt_rows: list[dict]) -> list[dict]:
    records: list[dict] = []
    metric_keys = base.METRIC_KEYS + EXTRA_KEYS
    for voltage, lq, nt in zip(VOLTAGES, lq_rows, nt_rows):
        row = {
            "voltage_V": float(voltage), "lq_source": lq["source"],
            "nt_source": nt["source"], "nt_status": nt["case_acceptance"],
            "nt_solve_complete": nt["solve_complete"],
            "last_accepted_joule_factor": nt["last_accepted_joule_factor"],
            "last_attempt_residual": nt.get("last_attempt_residual", ""),
            "failure_message": nt.get("failure_message", ""),
        }
        for key in metric_keys:
            lv, nv = lq.get(key, ""), nt.get(key, "")
            row[f"lq_{key}"] = lv
            row[f"nt_{key}"] = "" if isinstance(nv, float) and not np.isfinite(nv) else nv
            row[f"relative_gap_{key}"] = (
                (float(nv) - float(lv)) / abs(float(lv))
                if isinstance(lv, float) and lv != 0.0 and isinstance(nv, float) and np.isfinite(nv)
                else ""
            )
        records.append(row)
    with (OUT / "F2_voltage_sweep_0_10V_source.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    return records


def finite_nt(nt_rows: list[dict], key: str) -> tuple[np.ndarray, np.ndarray]:
    x, y = [], []
    for voltage, row in zip(VOLTAGES, nt_rows):
        value = row.get(key, np.nan)
        if isinstance(value, float) and np.isfinite(value):
            x.append(float(voltage)); y.append(value)
    return np.asarray(x), np.asarray(y)


def plot_load_response(lq_rows: list[dict], nt_rows: list[dict]) -> None:
    specs = [
        ("theta_max_K", 1.0, "最大温升 (K)"),
        ("joule_W", 1e3, "总焦耳功率 (mW)"),
        ("electron_right_A", 1e3, "右电极电子电流 (mA)"),
        ("ring_E_mean_V_m", 1e-5, "孔周平均电场 (kV/cm)"),
        ("ring_Q_mean_W_m3", 1e-12, r"孔周平均焦耳热 ($10^{12}$ W m$^{-3}$)"),
        ("minimum_mu_n_ratio", 1.0, r"最小迁移率比 $\mu_n/\mu_{n0}$"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(7.35, 4.65), layout="constrained")
    for ax, (key, scale, ylabel) in zip(axes.flat, specs):
        lq = np.array([float(r[key]) for r in lq_rows]) * scale
        ax.scatter(VOLTAGES, lq, s=22, color="#2676D2", zorder=4,
                   label="线性 LQ（精确比例缩放点）")
        x, y = finite_nt(nt_rows, key); y *= scale
        accepted = x <= 3.5 + 1e-12
        ax.plot(x[accepted], y[accepted], color="#D94B40", lw=1.7,
                label="非线性 NT（完整解）")
        flagged = x >= 3.5 - 1e-12
        ax.plot(x[flagged], y[flagged], color="#D94B40", lw=1.7, ls="--",
                label="NT完整解；热流审计未通过")
        ax.axvline(6.0, color="#777777", lw=.7, ls=":")
        ax.set(xlabel="电压幅值 (V)", ylabel=ylabel, xlim=(-.15, 10.15),
               xticks=np.arange(0, 11, 2))
        ax.margins(y=.12)
    axes[0, 0].legend(loc="upper left")
    fig.suptitle("中心圆孔ZnO薄板：0–10 V请求范围与实际可得非线性解", fontsize=10.2)
    fig.text(.5, .004, "蓝点延伸至10 V；红线仅连接真实完整场解（当前止于6 V）",
             ha="center", fontsize=7.4, color="#555555")
    base.save(fig, "F2_0_10V_线性点_非线性实际解线_六指标")


def plot_gap(records: list[dict]) -> None:
    finite = [r for r in records[1:] if r["nt_solve_complete"]]
    x = np.array([float(r["lq_theta_max_K"]) for r in finite])
    specs = [("theta_max_K", "最大温升差"), ("joule_W", "总焦耳功率差"),
             ("ring_E_mean_V_m", "孔周平均电场差"),
             ("ring_Q_mean_W_m3", "孔周平均焦耳热差")]
    fig, ax = plt.subplots(figsize=(4.75, 3.4), layout="constrained")
    for (key, label), color in zip(specs, ["#D94B40", "#6A51A3", "#238B45", "#CB6D2D"]):
        gap = 100 * np.array([float(r[f"relative_gap_{key}"]) for r in finite])
        accepted = np.array([float(r["voltage_V"]) <= 3.5 for r in finite])
        ax.plot(x[accepted], gap[accepted], lw=1.5, color=color, label=label)
        flagged = np.array([float(r["voltage_V"]) >= 3.5 for r in finite])
        ax.plot(x[flagged], gap[flagged], lw=1.5, ls="--", color=color)
    ax.axhline(0, color="#444444", lw=.65)
    ax.set(xlabel="线性LQ预测的最大温升 (K)", ylabel="NT相对LQ偏差 (%)")
    ax.legend(loc="best")
    ax.text(.98, .04, "仅统计完整场解；虚线未通过L0热流审计", transform=ax.transAxes,
            ha="right", va="bottom", color="#555555", fontsize=7.2)
    base.save(fig, "F2_非线性偏离随温升_实际完整点")


def plot_high_voltage_diagnostics(nt_rows: list[dict]) -> None:
    high = [(float(v), r) for v, r in zip(VOLTAGES, nt_rows) if v >= 6.0]
    fig, axes = plt.subplots(2, 2, figsize=(6.35, 4.9), layout="constrained")
    specs = [("theta_max_K", 1.0, "完整解最大温升 (K)"),
             ("heat_balance_relative", 100.0, "独立热流误差 (%)"),
             ("E_above_1e6_volume_fraction", 100.0, r"$|E|>10$ kV/cm体积分数 (%)")]
    for ax, (key, scale, ylabel) in zip(axes.flat[:3], specs):
        x = [v for v, r in high if np.isfinite(float(r.get(key, np.nan)))]
        y = [float(r[key]) * scale for v, r in high if np.isfinite(float(r.get(key, np.nan)))]
        ax.plot(x, y, color="#D94B40", marker="o", lw=1.45)
        ax.set(xlabel="电压幅值 (V)", ylabel=ylabel, xlim=(5.7, 10.3), xticks=range(6, 11))
        ax.margins(y=.15)
    axes[0, 1].axhline(1.0, color="#333333", ls="--", lw=.8, label="冻结门限 1%")
    axes[0, 1].legend(loc="best")
    ax = axes[1, 1]
    x = [v for v, _ in high]
    y = [100 * float(r["last_accepted_joule_factor"]) for _, r in high]
    colors = ["#D94B40" if r["solve_complete"] else "#6B6B6B" for _, r in high]
    ax.bar(x, y, color=colors, width=.58)
    ax.axhline(100, color="#333333", lw=.75, ls="--")
    ax.set(xlabel="电压幅值 (V)", ylabel="最后接受的焦耳热载荷 (%)",
           xlim=(5.5, 10.5), xticks=range(6, 11), ylim=(0, 108))
    for voltage, value, (_, row) in zip(x, y, high):
        label = "完整" if row["solve_complete"] else "失败"
        ax.text(voltage, value + 2, label, ha="center", va="bottom", fontsize=7)
    fig.suptitle("6–10 V高电压计算：完整解指标与续接可达性", fontsize=10.2)
    base.save(fig, "F2_6_10V_高电压求解状态")


def plot_clouds_6v() -> dict:
    xyz, conn, _ = base.read_mesh(HIGH_INPUTS / "H1_L0.mesh")
    lq3 = base.node_state(base.LQ_CASE)
    ratio = 6.0 / base.LQ_ANCHOR_V
    lq6 = lq3.copy()
    lq6[:, 0:3] *= ratio
    lq6[:, 3] *= ratio
    lq6[:, 4:6] = 0.0
    lq6[:, 6] *= ratio**2
    nt6_case = HIGH_RUNS / "F2_H1_NT_L0_V060"
    nt6_metrics = base.read_metrics(nt6_case, allow_heat_flux_fail=True)
    nt6 = base.node_state(nt6_case)
    lq = base.fe.plane(xyz, conn, lq6, "LQ", 0.0)
    nt = base.fe.plane(xyz, conn, nt6, "NT", 0.0)
    fields = [("phi", 1.0, "电势", "V", "viridis"),
              ("E", 1e-5, "电场强度", "kV/cm", "viridis"),
              ("Q", 1e-12, "焦耳热", r"$10^{12}$ W m$^{-3}$", "magma"),
              ("theta", 1.0, "温升", "K", "magma")]
    fig, axes = plt.subplots(4, 2, figsize=(6.7, 10.1), layout="constrained")
    for row, (key, scale, name, unit, cmap) in enumerate(fields):
        arrays = [lq[key] * scale, nt[key] * scale]
        bounds = (min(float(a.min()) for a in arrays), max(float(a.max()) for a in arrays))
        for col, (model, data, values) in enumerate(zip(("线性 LQ", "非线性 NT"), (lq, nt), arrays)):
            suffix = "" if col == 0 else "（L0热流审计未通过）"
            base.fe.panel(fig, axes[row, col], data, values,
                          f"{model}：{name}，6 V，z=0{suffix}", unit, cmap, bounds)
            base.add_hole(axes[row, col])
    fig.suptitle("最高完整点6 V场分布对照（每行共用色标范围）", fontsize=10.2)
    base.save(fig, "F2_6V_线性非线性云图对照")
    diffs = [(nt["E"] - lq["E"], 1e-5, "电场差 NT−LQ", "kV/cm"),
             (nt["Q"] - lq["Q"], 1e-12, "焦耳热差 NT−LQ", r"$10^{12}$ W m$^{-3}$"),
             (nt["theta"] - lq["theta"], 1.0, "温升差 NT−LQ", "K")]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.65), layout="constrained")
    for ax, (values, scale, title, unit) in zip(axes, diffs):
        scaled = values * scale
        limit = max(abs(float(scaled.min())), abs(float(scaled.max())))
        base.fe.panel(fig, ax, nt, scaled, title, unit, "RdBu_r", (-limit, limit))
        base.add_hole(ax)
    fig.suptitle("6 V非线性修正（最高完整点；L0热流审计未通过）", fontsize=10.2)
    base.save(fig, "F2_6V_NT减LQ差值云图")
    np.savez_compressed(OUT / "F2_6V_midplane_cloud_source.npz",
                        xyz=lq["xyz"], triangles=lq["triangles"],
                        lq_phi=lq["phi"], nt_phi=nt["phi"], lq_E=lq["E"], nt_E=nt["E"],
                        lq_Q=lq["Q"], nt_Q=nt["Q"], lq_theta=lq["theta"], nt_theta=nt["theta"])
    return nt6_metrics


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=False)
    base.OUT = OUT
    lq_rows, nt_rows = load_rows()
    records = export_csv(lq_rows, nt_rows)
    plot_load_response(lq_rows, nt_rows)
    plot_gap(records)
    plot_high_voltage_diagnostics(nt_rows)
    nt6 = plot_clouds_6v()
    contract = {
        "geometry": "20 um x 20 um x 1 um ZnO plate; centered through-hole r=3 um; c axis +z.",
        "boundary_conditions": "x faces: phi=0/-V, eta_n=eta_p=0, theta=0; all other faces natural; minimal rigid-body constraints.",
        "voltage_points_requested_V": [float(v) for v in VOLTAGES],
        "lq_provenance": "Exact scaling of the accepted 3 V H1-LQ FEM anchor; blue points are not independent solves.",
        "nt_provenance": "Red curves contain actual complete saved FEM states only. Failed continuation attempts are omitted from field curves and reported in the status panel.",
        "encoding": "LQ blue unconnected points through 10 V; NT red line only through the highest complete point; dashed NT segments fail the independent L0 heat-flux audit.",
        "acceptance_scope": "All high-voltage L0 states are exploratory. No failed point is interpolated, extrapolated, or shown as a solved field.",
        "nt_status_by_voltage": {str(float(v)): r["case_acceptance"] for v, r in zip(VOLTAGES, nt_rows)},
        "last_accepted_joule_factor_by_voltage": {str(float(v)): r["last_accepted_joule_factor"] for v, r in zip(VOLTAGES, nt_rows) if v >= 6.0},
        "six_volt_complete_state": {key: nt6.get(key) for key in base.METRIC_KEYS + EXTRA_KEYS},
    }
    (OUT / "figure_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"F2_HIGH_VOLTAGE_OUTPUT={OUT}")


if __name__ == "__main__":
    main()
