"""Audit completed corrected-geometry runs and archive their actual scope.

This script never reruns a solver, changes a mesh, or promotes historical plots.
All three input directories must already contain complete deterministic logs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from audit_g12_radius_output import METRICS, RADII, audit_values, parse
from audit_hex20_boundary_geometry import audit, read_mesh

CORE_GATES = (
    "V01", "G02", "G03", "G04", "L01", "G05", "G06", "G07", "G08",
    "G12C", "A1J", "G09", "G10", "G11",
)
EXPECTED_SOURCE_CHANGES = {
    "core/build_info_m.f90",
    "verification/radius_matrix_check_m.f90",
    "verification/verification_driver_m.f90",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_repeat(directory: Path, prefix: str) -> dict:
    one, two = (directory / f"{prefix}_run_{n}.txt" for n in (1, 2))
    if one.read_bytes() != two.read_bytes():
        raise ValueError(f"Repeat mismatch: {directory.name}/{prefix}")
    values = parse(one.read_text(encoding="utf-8-sig"))
    if values.get("verification_result") != "PASS":
        raise ValueError(f"Gate did not pass: {directory.name}/{prefix}")
    return values


def check_build(directory: Path) -> dict:
    manifest = json.loads((directory / "build_manifest.json").read_text(encoding="utf-8-sig"))
    checked = {}
    for entry in manifest["files"]:
        relative = entry["path"].replace("\\", "/")
        if relative.startswith("src/"):
            archived = directory / "source_snapshot" / relative[4:]
            checked[relative[4:]] = sha(archived)
        elif relative.endswith("/TPSolverV2.exe"):
            archived = directory / "TPSolverV2.exe"
        else:
            continue
        if sha(archived) != entry["sha256"]:
            raise ValueError(f"Build provenance mismatch: {archived}")
    if not checked or not (directory / "TPSolverV2.exe").is_file():
        raise ValueError(f"Build provenance incomplete: {directory}")
    return {
        "version": manifest["solver_version"],
        "executable_sha256": sha(directory / "TPSolverV2.exe"),
        "source_hashes": checked,
    }


def check_geometry(directory: Path, case_directories: list[Path]) -> list[dict]:
    results = []
    for case in case_directories:
        mesh_path = case / "mesh.dat"
        coords, connectivity, _ = read_mesh(mesh_path)
        checks = audit(coords, connectivity)
        if any(item["status"] != "PASS_SAMPLED" for item in checks.values()):
            raise ValueError(f"Independent geometry sampling failed: {mesh_path}")
        results.append({
            "case_directory": case.relative_to(directory).as_posix(),
            "case_sha256": sha(case / "case.nml"),
            "mesh_sha256": sha(mesh_path),
            "sampling_not_global_positivity_proof": True,
            "checks": checks,
        })
    return results


def hashes(directory: Path) -> None:
    records = [
        f"{sha(path)}  {path.relative_to(directory).as_posix()}"
        for path in sorted(directory.rglob("*"))
        if path.is_file() and path.name != "SHA256SUMS.txt"
    ]
    (directory / "SHA256SUMS.txt").write_text("\n".join(records) + "\n", encoding="ascii")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core", type=Path, required=True)
    parser.add_argument("--g12", type=Path, required=True)
    parser.add_argument("--g13", type=Path, required=True)
    args = parser.parse_args()
    core, g12, g13 = (path.resolve(strict=True) for path in (args.core, args.g12, args.g13))
    core_values = {gate: check_repeat(core, gate) for gate in CORE_GATES}
    radius_values = check_repeat(g12, "G12")
    radius_audit = audit_values(radius_values)
    if radius_audit["status"] != "PASS":
        raise ValueError(radius_audit)
    endpoint_values = {radius: check_repeat(g13, radius) for radius in ("R060", "R120")}
    builds = {label: check_build(directory) for label, directory in
              (("core", core), ("g12", g12), ("g13", g13))}
    if builds["core"] != builds["g12"]:
        raise ValueError("Core and G12 must use the exact same build")
    old_sources, new_sources = builds["g13"]["source_hashes"], builds["core"]["source_hashes"]
    if old_sources.keys() != new_sources.keys():
        raise ValueError("Unreviewed source-tree additions/removals between G13 and core")
    source_changes = {name for name in old_sources if old_sources[name] != new_sources[name]}
    if source_changes != EXPECTED_SOURCE_CHANGES:
        raise ValueError(f"Unreviewed source changes between G13 and core: {source_changes}")
    old_negative = (core / "old_mesh_negative_control.txt").read_text(encoding="utf-8-sig")
    if "Invalid HEX20 boundary-trace geometry" not in old_negative or "inspect_result=PASS" in old_negative:
        raise ValueError("Missing old-geometry rejection evidence")

    core_cases = sorted(path.parent for path in (core / "input_snapshot").rglob("case.nml"))
    g12_cases = [g12 / "input_snapshot" / radius for radius in ("H00", "R060", "R090", "R120")]
    g13_cases = [g13 / f"{radius}_{level}" for radius in ("R060", "R120") for level in ("L0", "L1", "L2")]
    if len(core_cases) != 6:
        raise ValueError("Expected six archived core cases")
    # G13 endpoints and G12 must represent exactly the same fine meshes.
    for radius in ("R060", "R120"):
        if sha(g12 / "input_snapshot" / radius / "mesh.dat") != sha(g13 / f"{radius}_L2" / "mesh.dat"):
            raise ValueError(f"G12/G13 fine-mesh mismatch: {radius}")
    if sha(g12 / "input_snapshot/R090/mesh.dat") != sha(core / "input_snapshot/G10_mesh_convergence/L2/mesh.dat"):
        raise ValueError("G12/G10 R090 fine-mesh mismatch")
    if sha(g12 / "input_snapshot/H00/mesh.dat") != sha(core / "input_snapshot/G11_matched_control/H00_L2/mesh.dat"):
        raise ValueError("G12/G11 H00 mesh mismatch")

    geometries = {
        "core": check_geometry(core, core_cases),
        "g12": check_geometry(g12, g12_cases),
        "g13": check_geometry(g13, g13_cases),
    }
    summary = {
        "verification_status": "VERIFIED_FOR_DECLARED_CLI_GATES",
        "whole_project_publication_package": "NOT_CURRENTLY_VERIFIED",
        "directories": {name: str(path) for name, path in (("core", core), ("g12", g12), ("g13", g13))},
        "builds": builds,
        "source_changes_between_g13_and_core": sorted(source_changes),
        "core_exact_repeat_gates": list(CORE_GATES),
        "g12_output_audit": radius_audit,
        "g12_exact_repeat": True,
        "g13_exact_repeat": {"R060": True, "R120": True},
        "old_mesh_rejected_before_newton": True,
        "geometry": geometries,
        "new_visual_studio_cross_build_comparison": "NOT_RUN",
        "paper_figures_rebound_to_repaired_evidence": False,
        "independent_hole_boundary_flux_acceptance": "NOT_ESTABLISHED_BY_THESE_GATES",
    }
    (core / "repaired_evidence_audit.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    labels = (
        "总焦耳功率", "全域平均焦耳热", "环带焦耳热均值", "环带焦耳热 P95",
        "环带电场均值", "电子重分布指标", "空穴重分布指标", "环带平均温升", "环带平均位移模",
    )
    metric_lines = ["| 对比量 | 60 nm | 90 nm | 120 nm |", "|---|---:|---:|---:|"]
    for label, metric in zip(labels, METRICS):
        percents = [100 * float(radius_values[f"g12_r{r}_contrast_{metric}"]) for r in RADII]
        metric_lines.append("| " + label + " | " + " | ".join(f"{v:+.4f}%" for v in percents) + " |")
    mesh_keys = (
        "total_joule_power", "annulus_joule_mean", "annulus_joule_p95",
        "annulus_electric_mean", "annulus_n_relative_mean", "annulus_temperature_mean",
    )
    mesh_labels = ("总功率", "环带焦耳均值", "焦耳 P95", "电场均值", "电子重分布", "温升均值")
    mesh_lines = ["| L1→L2 相对变化 | 60 nm | 90 nm | 120 nm |", "|---|---:|---:|---:|"]
    mesh_values = (endpoint_values["R060"], core_values["G10"], endpoint_values["R120"])
    for label, key in zip(mesh_labels, mesh_keys):
        percents = [100 * float(values[f"g10_{key}_change_l1_l2"]) for values in mesh_values]
        mesh_lines.append("| " + label + " | " + " | ".join(f"{v:.4f}%" for v in percents) + " |")
    report = f"""# 修复后孔洞数值重验证报告 — 2026-09-06

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: deterministic reproducibility validation
- Origin Date: 2026-09-06（UTC+08:00；带 Z 的运行目录名使用 UTC）
- Verification Status: VERIFIED，仅限下述 CLI 门禁
- Version Label: core/G12 `{builds['core']['version']}`；G13 `{builds['g13']['version']}`

## 结论与边界

修复径向中节点后，G09/G10/G11/G12/G13 以及基础门禁、A1J 均通过各自冻结判据，
并完成字节完全一致的双运行。旧翻折网格在求解前被拒绝。
G12 的逐半径参考值保存和验收遗漏已修复；全部 27 个对比量通过独立算术核对。

**这不等于整套论文结果已最终冻结。** 新版本 Visual Studio 跨构建复现未重做；
旧 `final_results` 场图和数据仍然过期；这些门禁没有给出实际孔洞独立边界通量的
预冻结验收阈值，也不是材料实验验证。G10/G13 验证既定积分方式下六项指标的网格稳定性，
不验证未经积分的点值极大值，也不单独排除环带后处理积分误差。

## 版本与来源

- 核心双运行：`{core.name}`，14 项门禁。
- G12 双运行：`{g12.name}`，第二次读取归档输入。
- G13 双运行：`{g13.name}`，两个端点各两次三级网格运行。
- G13 与新版源码只差版本信息、G12 审计模块、验证入口三处；
  控制方程、本构、单元、Newton、网格检查和 G10 路径源码逐文件一致。
- 已核对两个版本的归档源码/可执行文件与构建清单 SHA-256；不把两种构建冒充同一二进制。
- G12 与 G10/G13 的对应最细网格哈希一致。
- 共 16 份归档算例通过独立 3 阶/5 阶体积分点、5 阶边界迹及 9×9×9 闭合采样。
  这是几何采样证据，不是参数域连续 Jacobian 正性的数学证明。

## 网格敏感性

{chr(10).join(mesh_lines)}

沿用既定阈值：总功率 <1%，焦耳 P95 <10%，其余四项 <5%；
每个半径至少四项随加密改善。未根据新结果放宽阈值。

## 修复后半径矩阵

共同输入仍为 `600 nm × 600 nm × 20 nm`、固定 `−5 V`、
`TPS_G08_REFERENCE`、M2、孔外 30 nm 环带。表格为
`100×[H02(r)−H00(相同虚拟环带 r)]/H00(相同虚拟环带 r)`。
它是参考材料数值算例，不是实验标定材料参数下的定量预测。

{chr(10).join(metric_lines)}

本表支持的现象是：在这些固定条件下，局部焦耳热和载流子重分布指标增大，
而总功率与环带平均温升可以下降；两者并不矛盾。但表格本身不能把各耦合路径的因果贡献分离，
也不能把半径趋势推广为普适规律。电子/空穴指标是相对浓度偏移的汇总指标，不是净载流子总数。

## 本轮新增防错检查

- G12C：90 项无效输入均被拒绝，双运行一致。
- 独立输出审计：27 个对比量自洽；
  最大归一化算术误差 `{radius_audit['maximum_scaled_contrast_discrepancy']:.3e}`。
- 审计器自测：119 项负控制，含原始“各半径共用 R120 打印基准”的错误。
- 旧网格负控制：`old_mesh_negative_control.txt`，在 Newton 前拒绝。
- 机器可读证据：`repaired_evidence_audit.json`；完整文件校验值：各目录 `SHA256SUMS.txt`。

## 尚需闭环

1. 在当前源码上重新做 Visual Studio/CLI 跨构建比较。
2. 对实际孔洞生产输出的独立边界通量进行专门验收，不能把组装守恒残差替代这一项。
3. 新数据与新图源绑定、可视化检查后，才能替换历史 `final_results`。
"""
    (core / "repaired_validation_report.md").write_text(report, encoding="utf-8")
    for directory in (core, g12, g13):
        hashes(directory)
    print("REPAIRED_CLI_GATES=VERIFIED")
    print("WHOLE_PROJECT_PACKAGE=NOT_CURRENTLY_VERIFIED")
    print("REPORT=" + str(core / "repaired_validation_report.md"))


if __name__ == "__main__":
    main()
