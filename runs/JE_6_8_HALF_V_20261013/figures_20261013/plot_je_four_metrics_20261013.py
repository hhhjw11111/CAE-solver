"""Plot four existing voltage-response metrics without extrapolating NT cases.

The same-mesh figure compares original drift heat and current work only.
The overview includes a separately identified LQ reference on its original mesh.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import MultipleLocator, MaxNLocator

SOLVER = Path(__file__).resolve().parents[1]
ROOT = SOLVER.parent
RUN = SOLVER / "runs/JE_6_8_HALF_V_20261013"
DEFAULT_OUTPUT = RUN / "figures_20261013"
ORIGINAL = SOLVER / "runs/HV_SUPG_20261007/r13_XY0_Z2_V060"
LQ = SOLVER / "runs/HV_SUPG_20261007/production_r25_clean/candidate"
MESH = SOLVER / "cases/HV_SURFACE_20261008/XY0_Z2/H1_XY0_Z2.mesh"
BLUE = "#1675d1"
RED = "#d63c32"
GREY = "#4d5660"
SPECS = [
    ("theta_max_K", 1.0, 2, r"最大温升 $\theta_{\max}$ (K)", "maximum_temperature"),
    ("joule_W", 1e3, 2, r"总热源功率 $P_Q$ (mW)", "total_heat_power"),
    ("ring_E_mean_V_m", 1e-5, 1, r"孔周平均电场 $\langle |E|\rangle$ (kV/cm)", "annulus_electric_field"),
    ("ring_Q_mean_W_m3", 1e-12, 2, "孔周平均热源密度 " + r"$\langle Q\rangle$" + "\n" + r"($10^{12}$ W m$^{-3}$)", "annulus_heat_density"),
]


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def metrics(directory):
    result = {}
    for line in (directory / "metrics.txt").read_text(encoding="utf-8-sig").splitlines():
        if "=" not in line:
            continue
        key, val = line.split("=", 1)
        try:
            result[key] = float(val.strip())
        except ValueError:
            result[key] = val.strip()
    return result


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def setup_style():
    for name in ["msyh.ttc", "msyhbd.ttc"]:
        font = Path("C:/Windows/Fonts") / name
        if font.exists():
            font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Microsoft YaHei", "Arial", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8.5, "legend.fontsize": 7.1,
        "axes.linewidth": 0.8, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "axes.unicode_minus": False,
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "mathtext.fontset": "dejavusans",
    })


def save_figure(fig, stem, tiff=True):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes:
        for label in [ax.xaxis.label, ax.yaxis.label]:
            box = label.get_window_extent(renderer)
            assert box.x0 >= -1 and box.y0 >= -1, (stem, label.get_text())
    formats = ["png", "svg", "pdf"] + (["tiff"] if tiff else [])
    for ext in formats:
        kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}
        fig.savefig(stem.with_suffix("." + ext), dpi=600 if ext == "tiff" else 360,
                    bbox_inches="tight", pad_inches=.07, **kw)
    assert "<text" in stem.with_suffix(".svg").read_text(encoding="utf-8")
    plt.close(fig)


def panel(ax, spec, original, lq, current, overview):
    key, scale, exponent, ylabel, _ = spec
    vx = np.array([v["voltage_V"] for v in current])
    vy = np.array([v[key] for v in current]) * scale
    control_color = GREY if overview else BLUE
    if overview:
        dense = np.linspace(0, 8, 321)
        reference = lq[key] * (dense / 3) ** exponent * scale
        ax.plot(dense, reference, "--", color=BLUE, lw=1.35,
                label="原热源 LQ（L0参考）", zorder=2)
    ax.plot(vx, vy, "-o", color=RED, lw=1.4, ms=3.7,
            mfc="white", mew=1.0, label="J·E 候选 NT（XY0_Z2）", zorder=4)
    ax.plot([6.0], [original[key] * scale], "D", color=control_color,
            ms=4.5, mfc="white", mew=1.1,
            label="原热源 NT（同网格，6 V）", zorder=5)
    values = list(vy) + [original[key] * scale]
    if overview:
        values.append(lq[key] * (8 / 3) ** exponent * scale)
    ax.set(xlim=(-.02, 8.13) if overview else (5.9, 8.13),
           ylim=(0, max(values) * 1.13), xlabel=r"电压幅值 $V$ (V)", ylabel=ylabel)
    ax.xaxis.set_major_locator(MultipleLocator(1.0 if overview else .5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=4))
    ax.tick_params(length=3.2, width=.8, labelsize=8)
    ax.legend(loc="upper left", frameon=False, handlelength=2.3,
              borderaxespad=.35, labelspacing=.45)


def make_grid(out, original, lq, current, overview):
    fig, axes = plt.subplots(2, 2, figsize=(183 / 25.4, 140 / 25.4))
    fig.subplots_adjust(left=.115, right=.985, bottom=.15, top=.865,
                        wspace=.34, hspace=.39)
    title = "现有数据四指标：0–8 V 参考视图" if overview else "同网格四指标：6–8 V 实际计算点"
    subtitle = ("LQ为原热源、不同网格参考；J·E仅有6–8 V数据"
                if overview else "两组均为NT；蓝色为原热源6 V对照，红色为J·E候选分支")
    fig.suptitle(title, y=.973, fontsize=11)
    fig.text(.5, .925, subtitle, ha="center", fontsize=8, color=GREY)
    for idx, (ax, spec) in enumerate(zip(axes.flat, SPECS)):
        panel(ax, spec, original, lq, current, overview)
        ax.text(0, 1.06, "abcd"[idx], transform=ax.transAxes,
                fontsize=9, fontweight="bold")
    note = ("LQ按3 V结果缩放；NT只连接实际计算点，未补造低压或高压解。"
            if overview else "圆点为6.0、6.5、7.0、7.5、8.0 V；原热源高压缺失点未外推。")
    fig.text(.5, .057, note, ha="center", fontsize=7.5, color=GREY)
    fig.text(.5, .025, "J·E为诊断分支；空间精度和电功闭合仍待验证。",
             ha="center", fontsize=7.5, color=GREY)
    stem = "four_metrics_reference_0_8V" if overview else "four_metrics_same_mesh_6_8V"
    save_figure(fig, out / stem)


def sample_check(directory, m):
    # Archived independent 5x5x5 samples, read without changing solver results.
    data = np.loadtxt(directory / "field_samples.csv", delimiter=",", skiprows=1,
                      usecols=(0, 1, 2, 3, 7, 8))
    assert data.shape == (3136 * 125, 6)
    assert np.isfinite(data).all() and (data[:, 3] > 0).all()
    integrated = np.sum(data[:, 3] * data[:, 4])
    assert abs(integrated / m["joule_W"] - 1) < 1e-8
    peak = int(np.argmax(data[:, 5]))
    assert abs(data[peak, 5] - m["theta_max_K"]) < 1e-7
    return data[:, :4], {
        "sample_count": len(data), "integrated_heat_W": float(integrated),
        "sampled_theta_max_K": float(data[peak, 5]),
        "theta_peak_location_um": (data[peak, :3] * 1e6).tolist(),
        "heat_density_min_W_m3": float(data[:, 4].min()),
        "heat_density_max_W_m3": float(data[:, 4].max()),
        "field_samples_sha256": sha(directory / "field_samples.csv"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    setup_style()
    original, lq = metrics(ORIGINAL), metrics(LQ)
    current = []
    for v in [6, 6.5, 7, 7.5, 8]:
        directory = RUN / f"V{round(v * 10):03d}"
        m = metrics(directory)
        assert m["model"] == "NT" and m["voltage_V"] == v
        assert m["case_acceptance"] == "PASS_NOT_MESH_ACCEPTANCE"
        assert m["final_voltage_factor"] == m["final_joule_factor"] == 1
        assert "DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE" in m["heat_source"]
        assert all(np.isfinite(m[key]) for key, *_ in SPECS)
        current.append(m)
    assert original["model"] == "NT" and original["voltage_V"] == 6
    assert original["case_acceptance"] == "PASS_NOT_MESH_ACCEPTANCE"
    assert lq["model"] == "LQ" and lq["voltage_V"] == 3
    assert lq["solution_mode"] == "THREE_STAGE_LQ_NO_RETURN_FEEDBACK"
    assert lq["case_acceptance"] == "PASS_NOT_MESH_ACCEPTANCE"
    records = []
    source_files = []
    for directory, m, branch in [(ORIGINAL, original, "ORIGINAL_NT_CONTROL")] + [
        (RUN / f"V{round(m['voltage_V'] * 10):03d}", m, "CURRENT_WORK_NT") for m in current
    ]:
        src = directory / "metrics.txt"
        records.append({
            "branch": branch, "voltage_V": m["voltage_V"],
            "theta_max_K": m["theta_max_K"], "heat_power_mW": m["joule_W"] * 1e3,
            "annulus_E_kV_cm": m["ring_E_mean_V_m"] * 1e-5,
            "annulus_Q_1e12_W_m3": m["ring_Q_mean_W_m3"] * 1e-12,
            "status": m["case_acceptance"], "mesh": "XY0_Z2",
            "heat_source": m["heat_source"], "source_metrics": str(src.relative_to(ROOT)),
            "metrics_sha256": sha(src),
        })
        source_files.extend([src, directory / "nodes.csv"])
    for v in np.arange(0, 8.01, .5):
        records.append({
            "branch": "ORIGINAL_LQ_SCALING_REFERENCE", "voltage_V": float(v),
            "theta_max_K": lq["theta_max_K"] * (v / 3) ** 2,
            "heat_power_mW": lq["joule_W"] * (v / 3) ** 2 * 1e3,
            "annulus_E_kV_cm": lq["ring_E_mean_V_m"] * v / 3 * 1e-5,
            "annulus_Q_1e12_W_m3": lq["ring_Q_mean_W_m3"] * (v / 3) ** 2 * 1e-12,
            "status": "EXACT_MODEL_SCALING_NOT_INDEPENDENT_SOLVE", "mesh": "L0_REFERENCE_ONLY",
            "heat_source": lq["heat_source"],
            "source_metrics": str((LQ / "metrics.txt").relative_to(ROOT)),
            "metrics_sha256": sha(LQ / "metrics.txt"),
        })
    write_csv(out / "plot_data.csv", records)
    old_points, old_check = sample_check(ORIGINAL, original)
    new_points, new_check = sample_check(RUN / "V060", current[0])
    assert np.array_equal(old_points, new_points), "The 6 V samples are not the same mesh/weights."
    q0, q1 = original["joule_W"], current[0]["joule_W"]
    t0, t1 = original["theta_max_K"], current[0]["theta_max_K"]
    explanation = {
        "same_mesh_and_sampling": True, "original_6V": old_check, "current_work_6V": new_check,
        "heat_power_change_percent": 100 * (q1 / q0 - 1),
        "peak_temperature_change_percent": 100 * (t1 / t0 - 1),
        "descriptive_peak_to_source_ratio_K_per_mW": {
            "original": t0 / (q0 * 1e3), "current_work": t1 / (q1 * 1e3),
            "change_percent": 100 * ((t1 / q1) / (t0 / q0) - 1),
            "interpretation": "Describes the solved source distribution and peak location; not a material constant.",
        },
        "fixed_current_work_6V_state_electrical_work_decomposition_W": {
            "drift_work": q1 - current[0]["diffusion_work_W"],
            "diffusion_work": current[0]["diffusion_work_W"], "total_current_work": q1,
            "interpretation": "Same-state diagnostic identity, not a separate coupled solve or causal isolation.",
        },
        "electrical_work_gap_percent": [
            {"voltage_V": m["voltage_V"], "gap": 100 * abs(m["total_current_work_W"] - m["terminal_power_W"]) /
             max(abs(m["total_current_work_W"]), abs(m["terminal_power_W"]))} for m in current
        ],
    }
    (out / "temperature_diagnostics.json").write_text(
        json.dumps(explanation, indent=2, ensure_ascii=False), encoding="utf-8")
    provenance = []
    for file in source_files + [LQ / "metrics.txt", MESH]:
        provenance.append({"file": str(file.relative_to(ROOT)), "bytes": file.stat().st_size, "sha256": sha(file)})
    write_csv(out / "source_sha256.csv", provenance)
    make_grid(out, original, lq, current, overview=True)
    make_grid(out, original, lq, current, overview=False)
    for spec in SPECS:
        fig, ax = plt.subplots(figsize=(105 / 25.4, 90 / 25.4))
        fig.subplots_adjust(left=.20, right=.97, bottom=.24, top=.91)
        panel(ax, spec, original, lq, current, overview=True)
        fig.text(.5, .105, "LQ为原热源、不同网格参考；J·E仅有6–8 V实际点。",
                 ha="center", fontsize=6.6, color=GREY)
        fig.text(.5, .055, "诊断分支；未作空间精度或物理认证。", ha="center", fontsize=6.6, color=GREY)
        save_figure(fig, out / spec[4], tiff=True)
    shutil.copyfile(__file__, out / Path(__file__).name)
    qa = {"all_plotted_NT_full_load_and_accepted": True, "no_NT_extrapolation": True,
          "LQ_scaling_explicitly_labelled": True, "same_mesh_6V_samples_equal": True,
          "sampled_power_and_peak_match_metrics": True, "editable_svg_text": True,
          "label_bounds_check": "PASS", "visual_inspection": "PENDING"}
    (out / "qa_checks.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(out), "current_work_8V": records[5],
                      "temperature_diagnostics": explanation}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
