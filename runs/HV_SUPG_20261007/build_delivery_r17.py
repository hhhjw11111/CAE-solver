"""Assemble a traceable, explicitly incomplete campaign report from real evidence."""
from pathlib import Path
import csv
import json
import sys

RUN = Path(__file__).absolute().parent
ROOT = RUN.parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from run_hv_sweep import classify, read_metrics, sha
from audit_hv_stabilization import METRICS


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def main():
    out = RUN / 'delivery_r17'
    out.mkdir(exist_ok=False)
    snapshot = read(RUN / 'r17_source_snapshot/snapshot_manifest.json')['files']
    assert all(sha(ROOT / r['path']) == r['sha256'] for r in snapshot)
    py = read(RUN / 'python_verification_r17.json')
    assert py['exit_code'] == 0
    assert all(sha(ROOT / p) == h for p, h in py['files'].items())
    cli = read(RUN / 'regression_r17/report.json')
    assert len(cli) == 16 and all(r['passed'] and r['byte_identical'] for r in cli)
    prod = read(RUN / 'production_r17_serial/numerical_comparison.json')
    assert all(r['exit_code'] == 0 for r in read(RUN / 'production_r17_serial/executions.json'))
    assert prod['same_topology'] and all(s == 'PASS_NOT_MESH_ACCEPTANCE' for s in prod['old_frozen_acceptance'])

    # Two completed same-mesh states can disprove the tolerance, but cannot
    # complete the requested three-scale study or establish a spatial certificate.
    pair_dirs = [RUN / n for n in ['r17_fine_scale05_V060', 'r17_XY1_Z3_V060']]
    pair_metrics, pair_manifests = [], []
    for directory in pair_dirs:
        state, m = classify(directory, read(directory / 'execution.json')['exit_code'])
        assert state == 'SOLVE_FLUX_PASS_MESH_PENDING'
        assert read(directory / 'field_audit.json')['passed'] is True
        manifest = read(directory / 'manifest.json')
        for p, h in [('mesh', 'mesh_sha256'), ('executable', 'executable_sha256')]:
            assert sha(Path(manifest[p])) == manifest[h]
        pair_metrics.append(m)
        pair_manifests.append(manifest)
    for k in ['mesh_sha256', 'executable_sha256', 'model']:
        assert pair_manifests[0][k] == pair_manifests[1][k]
    for k in ['voltage_V', 'model', 'source', 'heat_source', 'supg_directional_tau']:
        assert pair_metrics[0][k] == pair_metrics[1][k]
    assert [m['supg_scale'] for m in pair_metrics] == [.5, 1.]
    comparisons = {}
    for k in METRICS:
        a, b = [m[k] for m in pair_metrics]
        span = abs(a-b)/max(abs(a), abs(b), 1e-300)
        comparisons[k] = dict(values=[a, b], relative_span=span, limit=.01, passed=span < .01)
    assert comparisons['ring_abs_p_relative_mean']['passed'] is False
    pair = dict(scope='TWO_COMPLETED_SCALES_DISPROVE_TOLERANCE_NOT_A_COMPLETE_THREE_SCALE_STUDY',
                scales=[.5, 1.], comparisons=comparisons, sensitivity_requirement_satisfied=False,
                scale2_complete=False, complete_quantitative_0_8V=False,
                source_runs=[str(p) for p in pair_dirs])
    (out / 'stabilization_two_scale_failure.json').write_text(json.dumps(pair, indent=2, allow_nan=False))

    with (RUN / 'trajectory_r13_L0_serial/status.csv').open(newline='') as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 17 and [float(r['voltage']) for r in rows] == [i/2 for i in range(17)]
    # Keep the original coarse table unchanged; this index separately identifies
    # the better 6 V evidence and does not claim a common-mesh curve.
    for r in rows:
        r['evidence'] = str(RUN / 'trajectory_r13_L0_serial/status.csv')
        r['mesh_scope'] = 'L0_DIAGNOSTIC_NOT_SPATIAL_ACCEPTANCE'
        r['stabilization_status'] = 'PENDING'
        r['fully_validated'] = False
        if float(r['voltage']) == 6:
            directory = RUN / 'r17_XY1_Z3_V060'
            state, m = classify(directory, read(directory / 'execution.json')['exit_code'])
            for k in list(r):
                if k in m:
                    r[k] = m[k]
            r.update(state=state, exit_code=0, evidence=str(directory), mesh_scope='XY1_Z3_MESH_PENDING',
                     stabilization_status='FAILED_COMPLETED_0.5_VS_1_COMPARISON')
    with (out / 'status_0_8V.csv').open('x', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    branch = read(RUN / 'fixed_voltage_branch_replay_r17.json')
    assert all(r['node_bytes_identical'] and r['state'] == 'SOLVE_FLUX_PASS_MESH_PENDING' for r in branch)
    for r in branch:
        assert sha(RUN / f"r17_inverse_XY1_Z3/accepted_{r['point']:05d}_nodes.csv") == r['source_sha256']
        assert sha(RUN / f"r17_fine_replay_branch_{r['point']:02d}/nodes.csv") == r['replay_sha256']
    summary = dict(version=(ROOT / 'VERSION').read_text().strip(),
                   complete_quantitative_0_8V=False, preserved_legacy_cli_checks=16,
                   legacy_cli_byte_identical=True, python_test_exit=py['exit_code'],
                   production_strict_byte_identical=prod['strict_byte_identical'],
                   production_maximum_scaled_node_difference=prod['maximum_scaled_node_difference'],
                   validated_source_snapshot_records=len(snapshot),
                   highest_fixed_voltage_solve_and_flux_pass_V=max(r['voltage_V'] for r in branch),
                   exact_fold_voltage_verified=False, dynamical_stability_verified=False,
                   material_experimental_calibration_verified=False,
                   missing_grid_voltages_V=[float(r['voltage']) for r in rows if r['state']=='NO_ACCEPTED_CHECKPOINT'],
                   limitation='Mixed-mesh evidence index; no common-mesh full-interval curve or missing-voltage interpolation.')
    (out / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False))

    table = '\n'.join(f"| {float(r['voltage']):g} | {r['state']} | "
                      f"{100*float(r['heat_balance_relative']):.5f}% |" if r['heat_balance_relative'] != ''
                      else f"| {float(r['voltage']):g} | {r['state']} | 无完整结果 |" for r in rows)
    report = f'''# 0–8 V 程序完善：当前交付与未通过项

程序已加入并验证高压数值功能；**0–8 V 全区间可信结果尚未完成，7–8 V 不能作为已计算的 NT 数据使用。**
当前版本为 `{summary['version']}`，本报告是实际运行证据汇总，不是全区间验收证书。

## 已落地的程序改动

- 默认关闭的一致残差 SUPG，包含曲面 HEX20 二阶导数和一致 Jacobian；高 Pe 双载流子解析检验及切线检查已有记录。
- 可选全方向对数步长控制、停滞回滚、割线预测、固定 0.5 V 实际求解检查点；旧输入继续使用原默认路径。
- 独立温度参数续接入口，保留所有原 PDE，以真实反求电压保存样本；需要普通定电压重求才能进入流量审核。
- 独立场/热流/电子/空穴审核，以及网格和稳定化敏感性检查。缺失、NaN、超限和未达到完整热源的点不能标为通过。
- 修复大网格边界初始化的自动临时数组问题；并行装配保持原单元散加次序。物理热源、R=0 和旧门限没有改动。

源码、项目文件和构建脚本的 **{len(snapshot)}** 个快照记录在本次汇总时逐一核对一致。
运行版本使用 `build/HV_SUPG_20261008_R17/TPSolverV2.exe`，没有覆盖原构建。

## 旧算例保护的实际证据

当前 R17 的 **16 项旧 CLI 验证全部通过，输出与任务开始前基线逐字节相同**。
最新 15 项 Python 检查通过，汇总时再次核对其脚本及测试哈希，没有代替或放宽数值门限。

最终构建另行复算同一个 3 V LQ 生产算例。新旧运行均通过原验收，
节点按冻结输入尺度归一的最大差为 **{prod['maximum_scaled_node_difference']:.6e}**。
生产文件严格字节比较结果为 **{'PASS' if prod['strict_byte_identical'] else 'FAIL'}**，保留原判定。
此前同一旧可执行程序自身重复也出现约 `5.31e-14` 的尺度差和字节差；不能把浮点末位接近称作生产文件逐字节复现。
这是一组代表性生产回归，不意味着重跑了所有历史网格。历史输入、结果和构建均保留。

## 真实电压点状态

下表是证据索引：0–5.5 V 采用 R13 L0 诊断，6 V 单列更充分的 R17 XY1_Z3 证据。
**不能将混合网格的这张表直接当成同一网格的定量曲线。**
`SOLVE_FLUX_PASS_MESH_PENDING` 仅表示完整目标求解和独立流量通过；网格/稳定化仍需另验。
`SOLVED_FLUX_GATE_FAILED` 表示已有完整代数解但流量门限失败；缺失点没有补值。

| 电压 / V | 求解与流量状态 | 独立热流误差 |
| --- | --- | --- |
{table}

另取得两个不在半伏网格上的细网格状态：**6.4593009 V** 和 **6.3729387 V**，
控制温升分别为 31.91655 K 和 37.40680 K。普通定电压重求后节点场与反求检查点
逐字节相同，独立热流误差分别为 **0.47621%** 和 **0.67897%**。
温度上升时电压下降，表明已采样到回折分支；6.4593009 V 不是已认证的精确临界值。
不同网格使用的控制节点位置不同，采样最高值也不能直接作为严格的临界电压网格误差。
这些稳态样本不证明真实动态热失稳，不排除尚未搜索的其它高温分支。

## 尚未通过的关键项

1. **载流子场敏感性。** 6 V 细网格上，SUPG 系数 0.5/1 的两个完整状态都通过流量门限。
   最高温升变化 {100*comparisons['theta_max_K']['relative_span']:.5f}%，端电流变化
   {100*comparisons['electron_right_A']['relative_span']:.5f}%，但空穴环带指标变化
   **{100*comparisons['ring_abs_p_relative_mean']['relative_span']:.5f}%**，超过预定 1%。
   系数 2 的试算未完成，已在两点明确超限后停止并记录原因；不能称为完成了三系数对照。
2. **空间误差。** 6 V 细网格与完整 L2 的温升/电流差约 0.219%/0.542%，空穴环带仍差约 9.20%。
   两网格在不同方向加密，这不是正式三层验收通过。
3. **高压可达性与模型范围。** 多种升压路径在 6.5 V 附近遇到困难，温度续接显示回折；
   尚无完整 6.5/7/7.5/8 V 检查点。当前材料仍是未实验标定的候选输入，热源仍是漂移电导热。
   屏蔽层中漂移与扩散抵消时的热源适用性是原材料审计 P05 未解决项，不能靠换一个热源公式补齐曲线。

## 确定的后续方案与通过条件

先解决 6 V 的载流子空间误差，再扩展电压。用固定面内等级的 Z2/Z3/Z4 序列分离厚度误差，
用固定厚度等级的 XY0/XY1/XY2 序列分离面内误差，重点核对表面层、孔边及接触邻域。
各状态必须重新求解；跨网格文件只用作初猜。现有 32 GB 主机不能同时运行两个大型直接分解。
当原网格门限和同网格 0.5/1/2 系数的 1% 敏感性门限都通过后，再使用该离散方案补跑半伏点。

对约 6.3–6.6 V 的分支，采用相同物理控制位置、小温度步或经过验证的弧长续接细化，
每个候选点都做普通定电压复算及独立流量检查。若仍回折，应继续区分数值分支、其它稳态分支
与动力学稳定性；不能预设一定存在一条平滑延伸到 8 V 的稳态曲线。

若目标包括真实材料在 7–8 V 的物理预测，还需核实适用温区/场强内的迁移率、热导率、
电极与散热边界，以及漂移—扩散的能量闭合。瞬态或高场材料模型需独立推导、输入和验证，
作为新模式加入，旧模式和旧结果继续保留。当前交付没有假装已完成这些物理扩展。

详细实现和逐阶段失败证据见 `docs/hv_0_8v_completion_plan_20261007.md`。
本目录的 `status_0_8V.csv`、`summary.json` 和 `stabilization_two_scale_failure.json`
提供机器可读状态；`evidence_manifest.json` 记录本次汇总引用文件的 SHA-256。
'''
    (out / '交付说明.md').write_text(report, encoding='utf-8')
    evidence = [RUN / p for p in ['regression_r17/report.json', 'python_verification_r17.json',
                'python_verification_r17.txt', 'source_build_audit_r17.json', 'source_changes_r17.patch',
                'production_r17_serial/report.json', 'production_r17_serial/numerical_comparison.json',
                'production_r17_serial/executions.json', 'fixed_voltage_branch_replay_r17.json',
                'branch_coordinate_audit_r17.json', 'six_volt_anisotropic_comparison_r17.json',
                'stabilization_fine_V060_r17.json', 'stabilization_coarse_global_V060_r17.json',
                'r17_fine_scale2_V060/stopped.json', 'r17_scale2_global_V060/stopped.json',
                'trajectory_r13_L0_serial/status.csv', 'r17_source_snapshot/snapshot_manifest.json']]
    for directory in pair_dirs + [RUN/'r17_fine_replay_branch_03', RUN/'r17_fine_replay_branch_04']:
        evidence.extend(directory / p for p in ['manifest.json', 'execution.json', 'metrics.txt', 'field_audit.json'])
    evidence.extend([ROOT/'build/HV_SUPG_20261008_R17/build_manifest.json',
                     ROOT/'build/HV_SUPG_20261008_R17/TPSolverV2.exe'])
    (out / 'evidence_manifest.json').write_text(json.dumps([
        dict(path=str(p), sha256=sha(p)) for p in evidence], indent=2))
    print(json.dumps(summary, indent=2))
    print(str(out / '交付说明.md'))


if __name__ == '__main__':
    main()
