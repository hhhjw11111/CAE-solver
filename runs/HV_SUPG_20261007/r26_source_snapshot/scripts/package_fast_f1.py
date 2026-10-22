"""Package accepted fast results, source records and a browsable image index."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import html
import json
import shutil


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--delivery-root", type=Path, required=True)
    args = parser.parse_args()
    root, delivery = args.run_root.resolve(), args.delivery_root.resolve()
    project = Path(__file__).resolve().parents[2]
    audit = read(root/"fast_audit.json")
    if not audit["case_numerical_checks_passed"] or not audit["repeat"]["representative_repeat_executed"]:
        raise ValueError("Cannot package incomplete cases or an unexecuted repeat")
    rows = read(root/"main_results.json")
    source, evidence = delivery/"04_数据表与图源", root/"evidence"
    if not (source/"visual_review.md").exists():
        raise ValueError("Visual review record required")
    audit["visual_review"] = "REVIEWED_REPORT_9_PAGES_AND_100_CURRENT_FIGURES"
    audit["known_limitations"] = [
        "Formal three-level L2 mesh acceptance not completed.",
        "Strict byte-identical repeat failed; no replacement acceptance threshold introduced.",
        "H0 annulus indicator quadrature underestimates exact annulus volume by about 28.48%.",
        "Unsmoothed stress maps show element-gradient seams; not physical discontinuities.",
    ]
    (root/"fast_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report_qa = read(root/"report_qa.json")
    report_qa["visual_review"] = "REVIEWED_9_PAGE_RENDER_AND_DETAIL_PAGES"
    (root/"report_qa.json").write_text(json.dumps(report_qa, indent=2), encoding="utf-8")
    gallery_qa = read(source/"gallery_qa.json")
    gallery_qa["visual_review"] = "REVIEWED_ALL_CURRENT_IMAGES_IN_7_CONTACT_SHEETS_AND_H1_NT_DETAIL"
    (source/"gallery_qa.json").write_text(json.dumps(gallery_qa, indent=2), encoding="utf-8")
    for name in ("fast_audit.json", "main_results.json", "main_results.csv", "report_qa.json"):
        shutil.copy2(root/name, source/name)
    for name in ("A1_joint_analytic", "A2_hole_analytic"):
        for ext in ("png", "svg", "pdf"):
            shutil.copy2(project/"output/figures/F1_20260906"/(name+"."+ext),
                         delivery/"03_组合对比图"/(name+"."+ext))
    for name in ("A1_error_source.csv", "A2_error_source.csv"):
        shutil.copy2(project/"output/figures/F1_20260906"/name, source/name)
    scripts = evidence/"delivery_scripts"
    scripts.mkdir(exist_ok=True)
    for name in [
        "finalize_fast_f1.py", "plot_fast_f1_gallery.py", "audit_fast_f1_repeat.py",
        "summarize_f1_release.py", "plot_f1_main_results.py", "audit_f1_results.py",
        "audit_hex20_boundary_geometry.py", "test_fast_f1_profiles.py",
        "run_fast_case.ps1", "run_fast_hole_matrix.ps1", "generate_f1_fast_cases.py",
        "package_fast_f1.py", "qa_fast_f1_artifacts.py",
    ]:
        shutil.copy2(Path(__file__).parent/name, scripts/name)
    for name in ("build_fast_f1_report.py", "build_case_freeze_pdf.py"):
        shutil.copy2(project/"scripts"/name, scripts/name)
    shutil.copy2(source/"visual_review.md", evidence/"visual_review.md")
    shutil.copy2(project/"tmp/pdfs/F1_fast_20260907/artifact_checks.json", source/"artifact_checks.json")
    shutil.copy2(project/"solver_v2/VERSION", evidence/"VERSION_corrected_after_build")
    shutil.copytree(project/"solver_v2/material_audit", evidence/"material_audit", dirs_exist_ok=True)
    manifest = read(evidence/"build_manifest.json")
    correction = {
        "record_type": "RETROSPECTIVE_METADATA_CORRECTION_NOT_REBUILD",
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "original_manifest_sha256": sha(evidence/"build_manifest.json"),
        "original_manifest_solver_version": manifest["solver_version"],
        "actual_embedded_version": "0.14.4-f1-fast-threaded",
        "actual_executable_sha256": sha(evidence/"TPSolverV2.exe"),
        "explanation": "VERSION metadata was corrected after the threaded build; original build manifest retained unmodified.",
        "early_H0_environment": {
            "record_type": "RETROSPECTIVE_FROM_EXECUTION_HISTORY",
            "cases": [r["case"] for r in rows if r["case"].startswith("FAST_H0")],
            "MKL_NUM_THREADS": 8, "OMP_NUM_THREADS": 8, "TPS_PARDISO_STORAGE": "IC",
            "MKL_DYNAMIC": "not explicitly set by the original launcher",
            "warning": "No contemporaneous per-case provenance.json was written for these three runs.",
        },
        "H1_N0_continuation_change": "Initial/maximum step changed to 1.0; minimum step/cutback/tolerances unchanged. Prior input retained in input_history.",
        "plot_sample_check_change": "Near-zero electric components normalized by the electric-vector norm; threshold unchanged. No PDE acceptance changed.",
        "controller_note": "The sequential controller stopped at already-existing FAST_H1_NT_L0 after five passes; the independently prioritized H1-NT passed and is included by case-level audit.",
    }
    (evidence/"metadata_correction.json").write_text(json.dumps(correction, indent=2), encoding="utf-8")
    shutil.copy2(evidence/"metadata_correction.json", source/"metadata_correction.json")
    for row in rows:
        case = row["case"]
        if sha(root/case/"input.nml") != sha(evidence/"inputs"/(case+".nml")):
            raise ValueError("Input snapshot mismatch: "+case)
    dirs = ["01_折线对比图", "02_单幅场图", "03_组合对比图"]
    counts = {d: len(list((delivery/d).glob("*.png"))) for d in dirs}
    h0 = next(r for r in rows if r["label"] == "H0-NT")
    h1 = next(r for r in rows if r["label"] == "H1-NT")
    table = "\n".join(
        f"| {r['label']} | {r['joule_W']*1e3:.6f} | {r['theta_max_K']:.6f} | "
        f"{r['ring_Q_mean_W_m3']/1e12:.6f} | {r['elapsed_seconds']/60:.2f} |"
        for r in rows)
    text = f"""# F1 六算例汇报成果

本目录按 Asia/Shanghai 运行日期 2026-09-07 归档。
六个 L0 算例均通过收敛与守恒验收，H1-NT 代表性复算已执行。
复算三个核心文件不逐字节一致，原严格复现门槛未通过，不能将本交付说成原全部门禁通过。
复算最大电势差 {audit['repeat']['maximum_absolute_field_difference'][7]:.4g} V，
最大温升差 {audit['repeat']['maximum_absolute_field_difference'][10]:.4g} K；
原场尺度下最大状态差 {audit['repeat']['maximum_solver_scaled_state_difference']:.4g}。
原三层 L2 网格验收未完成；候选 ZnO 未经实验定标。本目录不替代最终精度验收。

## 快速入口

- `F1_六算例汇报.pdf`：中文汇报报告。
- `图像索引.html`：打开后浏览全部图片，点击图片查看原尺寸。
- `01_折线对比图`：{counts[dirs[0]]} 张独立对比图，含六个整体指标和两条空间剖面。
- `02_单幅场图`：{counts[dirs[1]]} 张单幅场图，六个算例分别输出。
- `03_组合对比图`：{counts[dirs[2]]} 张组合图，其中 A1/A2 为保留的解析验证。
- `04_数据表与图源`：结果 CSV/JSON、两条剖面 CSV、18 个表面/中面 NPZ、审计记录。
- 每张图片均提供 PNG、SVG、PDF。SVG 中曲线与文字可编辑，场图着色层为栅格。

## 主结果

| 算例 | 总功率 mW | 最高温升 K | 环带平均 Q（10^12 W/m3） | 用时 min |
|---|---:|---:|---:|---:|
{table}

NT 中，有孔/无孔总功率比 {h1['joule_W']/h0['joule_W']:.6f}，
环带平均热源比 {h1['ring_Q_mean_W_m3']/h0['ring_Q_mean_W_m3']:.6f}；
最高温升差 {h1['theta_max_K']-h0['theta_max_K']:.6f} K。
整体功率下降与孔周局部热源集中可以同时出现，不应混为同一指标。
但当前 H0 环带积分点筛选的体积约低估 28.5%，局部平均值与 P95 存在粗网格采样偏差；
上述环带热源比只作当前离散趋势，不能作为已收敛定量结论。

## 图像与单位

- H0 为无孔，H1 为有孔；LQ 为弱耦合基准，N0 冻结迁移率和扩散系数，NT 恢复二者温度相关。
- 同名场图共用六算例色标。电子/空穴图为相对各自平衡浓度的百分比，并非绝对载流子数。
- 空穴参考浓度极低，相对变化不代表显著导电贡献。LQ 的微小载流子起伏不得解释成物理集中。
- `x_y4_z0`：沿 x，y=4 μm、z=0，整条直线避开孔洞。
- `y_x0_ztop`：沿 y，x=0、z=0.5 μm；有孔曲线在 -3<y<3 μm 处明确断开。
- 剖面 CSV 使用 SI 单位，`s_um` 除外；`n_change/p_change` 为无量纲相对变化，绘图乘 100。
- `in_material=0` 的 CSV 行是孔内空缺，不是零值；不做跨孔插值。
- 由原 HEX20 形函数评价，不人为平滑；导数在单元交界处取确定的单侧值。
- 应力图可见单元梯度接缝，不能把这些接缝当作物理间断；该图主要用于定位后续检查区域。
- 所有曲线是确定性数值结果，没有虚构误差条。模型之间连线不是参数扫描。
- NT−N0 为 μ(T) 与 D(T) 的共同反馈；N0−LQ 也不能归因于单一效应。

## 计算时间与优化边界

原 L1 正式计算约 57 分钟，包含多次非线性迭代与稀疏直接分解。
本轮采用 L0（H0:1960、H1:2240 个 HEX20 单元）、线程化直接求解、
内存内存储及可用提交内存检查。H1-N0 减少加载分段，仍保留自动回退和原容差。
前三个 H0 使用 8 线程，后续与复算使用 4 线程。部分算例并发运行；
表中单例耗时不能直接相加当作总墙钟时间，也不能用于宣称纯线程加速比。
未通过性能剖析证明磁盘读写是唯一或主要瓶颈，不作此结论。

## 证据与限制

- 原始输入、节点场、积分场、日志、编译信息、源码与二进制：
  `../../solver_v2/runs/F1_FAST_HOLE_L0_20260907/`。
- 原 L1 成功、L2 异常退出及更早验证：
  `../../solver_v2/runs/F1_20260906_R0144/`，全部保留。
- 原快速构建清单存在版本元数据滞后，原件未改；更正说明见 `metadata_correction.json`。
- 本次六算例全部完成是用户授权的缩减网格范围，不等于原 L2 冻结方案全部通过。
- 原收敛/守恒阈值未放宽。小于网格变化量的差异不可作精确机制结论。
"""
    (delivery/"README_先读我.md").write_text(text, encoding="utf-8")
    panels = []
    for folder in dirs:
        items = []
        for path in sorted((delivery/folder).glob("*.png")):
            relative = path.relative_to(delivery).as_posix()
            safe = html.escape(relative, quote=True)
            caption = html.escape(path.stem)
            items.append(f'<figure><a href="{safe}"><img loading="lazy" src="{safe}" alt="{caption}"></a>'
                         f'<figcaption>{caption}</figcaption><a href="{safe[:-3]}svg">SVG</a> · '
                         f'<a href="{safe[:-3]}pdf">PDF</a></figure>')
        panels.append(f'<section id="{folder}"><h2>{folder} · {len(items)} 张</h2>'
                      '<div class="grid">'+"".join(items)+"</div></section>")
    page = """<!doctype html><html lang="zh-CN"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>F1 六算例图像索引</title>
<style>body{font:15px Arial,"Microsoft YaHei",sans-serif;margin:0;color:#26323a;background:#fff}
header,main{max-width:1240px;margin:auto;padding:24px}header{border-bottom:1px solid #d9e0e3}
h1{font-size:25px;margin:0 0 12px}h2{font-size:19px;margin:24px 0 16px}p{line-height:1.65}
nav{display:flex;flex-wrap:wrap;gap:20px}a{color:#18736b}section{padding-bottom:32px}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,280px),1fr));gap:22px}
figure{margin:0;padding-bottom:14px;border-bottom:1px solid #e1e5e7;min-width:0}
img{display:block;width:100%;height:245px;object-fit:contain;background:white}
figcaption{font-size:13px;overflow-wrap:anywhere;margin:10px 0;line-height:1.5}</style>
<header><h1>F1 六算例图像索引</h1><p>L0 快速结果 · 六算例通过 · 复算有浮点末位差异，非逐字节一致 · 原 L2 网格验收未完成</p>
<nav><a href="F1_六算例汇报.pdf">中文汇报</a><a href="01_折线对比图">折线图文件夹</a>
<a href="#01_折线对比图">折线对比</a><a href="#02_单幅场图">单幅场图</a>
<a href="#03_组合对比图">组合对比</a><a href="04_数据表与图源/main_results.csv">结果 CSV</a></nav></header><main>
"""
    (delivery/"图像索引.html").write_text(page+"".join(panels)+"</main></html>", encoding="utf-8")
    files = {}
    for path in sorted(delivery.rglob("*")):
        if path.is_file() and path.name != "delivery_manifest.json":
            files[path.relative_to(delivery).as_posix()] = {"bytes": path.stat().st_size, "sha256": sha(path)}
    summary = {"status": audit["delivery_scope_status"], "images_by_folder": counts,
               "files": files, "recorded_utc": datetime.now(timezone.utc).isoformat()}
    (delivery/"delivery_manifest.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    record = {p.relative_to(evidence).as_posix(): sha(p)
              for p in sorted(evidence.rglob("*")) if p.is_file() and p.name != "evidence_manifest.json"}
    (evidence/"evidence_manifest.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps({"status": summary["status"], "images": counts, "files": len(files)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
