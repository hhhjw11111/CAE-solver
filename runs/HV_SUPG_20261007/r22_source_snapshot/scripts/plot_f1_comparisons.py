"""Six-case comparisons and qualified interpretation from actual F1 tables."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from summarize_f1_release import MAIN


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run-root", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()
    root, out = args.run_root.resolve(), args.output_dir.resolve()
    data = json.loads((root/"main_results.json").read_text(encoding="utf-8"))
    mesh = json.loads((root/"mesh_comparison.json").read_text(encoding="utf-8-sig"))
    if not mesh["passed"]:
        raise ValueError("Mesh audit is not accepted")
    lookup = {r["case"]: r for r in data}
    if set(lookup) != set(MAIN):
        raise ValueError("Expected exactly the frozen six main cases")
    families = ["LQ", "N0", "NT"]
    cases = {("H0", m): lookup[f"MAIN_H0_{m}_L2"] for m in families}
    cases.update({("H1", m): lookup[f"MAIN_H1_{m}_L2"] for m in ["LQ", "N0"]})
    cases["H1", "NT"] = lookup["MESH_H1_NT_L2"]
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 9, "legend.fontsize": 8,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .7,
        "svg.fonttype": "none", "pdf.fonttype": 42,
    })
    specs = [
        ("joule_W", 1e3, "Total Joule power", "mW"),
        ("theta_max_K", 1, "Peak temperature rise", "K"),
        ("ring_Q_mean_W_m3", 1e-12, "Annulus mean Joule heat", r"$10^{12}$ W m$^{-3}$"),
        ("ring_E_p95_V_m", 1e-5, "Annulus electric-field P95", "kV/cm"),
        ("ring_abs_n_relative_mean", 100, "Annulus mean electron deviation", "%"),
        ("right_mean_ux_m", 1e9, "Right-face mean displacement", "nm"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.6), layout="constrained")
    for ax, (key, scale, title, label), letter in zip(axes.flat, specs, "abcdef"):
        for geometry, color, marker in [("H0", "#42658A", "o"), ("H1", "#BB713D", "s")]:
            y = [cases[geometry, m][key]*scale for m in families]
            ax.plot(range(3), y, color=color, marker=marker, lw=1, ms=4, label=geometry)
        ax.set(xticks=range(3), xticklabels=families, title=title, ylabel=label)
        ax.text(-.17, 1.10, letter, transform=ax.transAxes, fontweight="bold")
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        ax.margins(x=.16, y=.18)
    axes[0, 0].legend(frameon=False)
    out.mkdir(parents=True, exist_ok=True)
    for ext in ["png", "svg", "pdf"]:
        fig.savefig(out/f"six_case_comparison.{ext}", dpi=360, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    contrasts = []
    for key, scale, title, unit in specs:
        pairs = [(f"H1-H0 / {m}", cases["H1", m], cases["H0", m]) for m in families]
        pairs += [(f"NT-N0 / {g}", cases[g, "NT"], cases[g, "N0"]) for g in ["H0", "H1"]]
        pairs += [(f"N0-LQ / {g}", cases[g, "N0"], cases[g, "LQ"]) for g in ["H0", "H1"]]
        for name, a, b in pairs:
            change = (a[key]-b[key])/max(abs(a[key]), abs(b[key]), 1e-300)
            uncertainty = mesh["metric_changes"].get(key, {}).get("L1_L2")
            contrasts.append({
                "contrast": name, "metric": key, "a": a[key], "b": b[key],
                "difference": a[key]-b[key], "signed_relative_difference": change,
                "H1_NT_mesh_change_reference": uncertainty,
                "resolution_screen": (
                    "NO_MATCHED_MESH_METRIC" if uncertainty is None else
                    "UNRESOLVED_AT_CURRENT_PRECISION" if abs(change) < 2*uncertainty else
                    "EXCEEDS_TWO_TIMES_H1_NT_MESH_CHANGE"
                ),
                "qualification": "H1-NT mesh change is a working screen, not a rigorous error bound for other cases.",
            })
    (root/"contrasts.json").write_text(json.dumps(contrasts, indent=2, allow_nan=False), encoding="utf-8")
    with (root/"contrasts.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(contrasts[0]))
        writer.writeheader()
        writer.writerows(contrasts)
    h0, h1 = cases["H0", "NT"], cases["H1", "NT"]
    paragraphs = [
        f"完整 NT 模型中，无孔与有孔总焦耳功率分别为 {h0['joule_W']*1000:.6g} 和 "
        f"{h1['joule_W']*1000:.6g} mW；最高温升分别为 {h0['theta_max_K']:.6g} 和 "
        f"{h1['theta_max_K']:.6g} K。环带热源均值的有孔/无孔比为 "
        f"{h1['ring_Q_mean_W_m3']/h0['ring_Q_mean_W_m3']:.6g}。这些量分别描述整体功率、"
        "峰值温度和局部热源，不能相互替代。",
    ]
    for g in ["H0", "H1"]:
        n0, nt = cases[g, "N0"], cases[g, "NT"]
        paragraphs.append(
            f"{g} 中 NT 相对 N0 的总功率变化为 {(nt['joule_W']/n0['joule_W']-1)*100:.5g}%，"
            f"最高温升变化为 {nt['theta_max_K']-n0['theta_max_K']:.6g} K。"
            "这一差值同时包含迁移率与扩散系数的温度相关性及场反馈，不能解释为纯迁移率效应。")
    unresolved = [r for r in contrasts if r["resolution_screen"] == "UNRESOLVED_AT_CURRENT_PRECISION"]
    if unresolved:
        paragraphs.append("以下差值小于两倍对应 H1-NT 网格变化量，当前数值精度下不可分辨：" +
                          "；".join(f"{r['contrast']} 的 {r['metric']}" for r in unresolved) + "。")
    paragraphs += [
        "LQ 等温场的载流子扰动接近舍入误差，不构成孔洞必然导致一阶载流子集中的证据。"
        "完整模型中的热反馈、空间电场与屏蔽共同决定载流子重分布。",
        "所有数值均属于冻结几何、边界、候选材料与漂移电导型焦耳热闭合下的条件性计算；"
        "不对应某个已实验标定 ZnO 样品的定量预测。",
    ]
    (root/"interpretation.json").write_text(json.dumps({"paragraphs": paragraphs}, ensure_ascii=False, indent=2), encoding="utf-8")
    (root/"interpretation.md").write_text("# F1 主结果解释\n\n"+"\n\n".join(paragraphs)+"\n", encoding="utf-8")
    contract = {
        "claim": "Compare geometry and model responses using the frozen six accepted deterministic FE solutions.",
        "archetype": "quantitative grid", "backend": "Python/matplotlib",
        "source": str(root/"main_results.csv"), "source_contrasts": str(root/"contrasts.csv"),
        "panels": [s[0] for s in specs], "width_inches": 7.2, "height_inches": 4.6,
        "error_bars": "None. These are deterministic solutions, not statistical estimates.",
        "integrity": "No additional parameter sweep, fabricated data, or statistical confidence claim.",
    }
    (out/"comparison_figure_contract.json").write_text(json.dumps(contract, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
