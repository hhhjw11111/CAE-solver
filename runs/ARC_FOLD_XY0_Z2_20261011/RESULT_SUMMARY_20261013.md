# 2026-10-13 弧长与 J·E 分支状态

## 原模型（drift-conductivity heat）

- 网格：`XY0_Z2`，13,920 节点、3,136 个 HEX20 单元。
- 伪弧长接受点：11 个，最大接受残差 `9.80e-16`。
- 折叠已检测并定位：
  - 电压因子：`0.8059180665`
  - 折叠电压：`6.447344532 V`
  - 左右电压因子：`0.80591804398`、`0.80591806653`
  - 弧长括区间：`1.9531e-3`（这是状态-载荷度量区间，不是电压误差上界）
- 独立热通量误差：`1.58449%`，高于 `1%` 门限，因此该折叠点是代数/分支诊断结果，不是完整生产验收通过点。
- 记录目录：[R10 原模型弧长结果](./ARC_FOLD_XY0_Z2_20261011_R10)

## J·E 独立热源分支

- 开关：`heat_source_current_work=.true.`；原模型使用 `.false.`。
- 从同网格 6 V 接受场重新平衡，初始再平衡最终残差 `4.20e-16`。
- 伪弧长接受点：电压因子 `0.8739578450`、`1.0342950391`。
- 因此 J·E 分支已经跨过 8 V 目标因子 `1.0`；程序已进入目标 8 V 固定电压回放。
- 当前回放残差从 `2.527e-3` 降至 `1.261e-4`，尚未达到 `1e-9`，所以暂不能称为完整 8 V 验收结果。
- 记录目录：[R12 J·E 分支](./ARC_FOLD_XY0_Z2_20261011_R12_JE_FROM6)

## 编译与回归

- 构建：[ARC_JE_FINAL_20261012](../../build/ARC_JE_FINAL_20261012/TPSolverV2.exe)
- `ARC`、`ARCS`、`G02JE`、`G04JE`、`F1MJE`、原模型 `G02/G04/G07/G08` 均通过。
- 接受点保存回调已接入伪弧长分支，R6/R7/R8/R9/R10/R12 均可从保存的 `accepted_*_nodes.csv` 继续。

### 2026-10-13 J·E target replay status

The J·E run reached an accepted arclength point at voltage factor `1.0342950391` (about `8.27436 V`) with residual `6.80e-16`; this is the first accepted point above the 8 V target. The automatic fixed-voltage replay at factor `1.0` reduced residual from `2.527e-3` to `1.231e-4` before the long run was stopped. It is therefore a branch-crossing result, not yet a fixed-8-V acceptance result.
