# CAE Solver Project State

数值核查截止：2026-10-13；Git 身份及发布说明修正：2026-10-21。下文数值结论为 10 月 13 日核查快照，不能代替后续开发状态。本次发布未重编译、重跑算例或重新颁发验收。

## 1. Current baseline

- 当前程序仅为 `solver_v2/src/`；旧 Console3 不属于当前基线。
- `VERSION` 与 `src/core/build_info_m.f90` 均为 **0.16.0-trust-region**；理论版本 2.1，架构 TPS-V2-ARCH-001。旧“两处版本号不一致”已不适用。
- 最新五点 J·E 交付摘要关联构建为 [ARC_JE_FINAL_20261012](build/ARC_JE_FINAL_20261012/build_manifest.json)。本次复核其 **60 份源码及 1 个可执行文件**，61 项 SHA-256 均与清单相符。可执行 SHA-256：`405e953c8cb0b1eca2a61eccc011945892924b43fca6805e87fd8fec19006d92`。
- 清单记录 Intel ifx、优化/线程构建、MKL LP64 threaded/PARDISO；精确安装版本及本机全新构建可复现性本次未验证。版本字符串沿用至后续 ARC/J·E 开发，不能单独识别源码或结果身份。
- Git：实际仓库根为 `solver_v2/`，远程为 `hhhjw11111/CAE-solver`，已有初始提交 `a63ad1db24a5d5813afa4fbbe7eaa6e7235950bc`。此前“尚无 commit”误读了外层 Console3 仓库，现予纠正。当前本地数值开发还有未提交改动，初始提交并不覆盖上述最新构建。

## 2. Current architecture

输入/检查 (`src/io`, `src/app`) → 网格与 HEX20 (`src/mesh`) → 材料/物理 (`src/physics`) → 单元 residual/Jacobian 与全局 CSR (`src/fem`) → PARDISO (`src/linear`) → Newton/缩放/延续 (`src/nonlinear`) → 场输出、独立通量和守恒审计 (`src/post`, `src/app`)。

- 原生验证入口：`src/verification/verification_driver_m.f90`；Python 测试：`tests/`；运行与审计脚本：`scripts/`。
- Newton、线搜索、voltage/Joule continuation 已存在；信赖域、伪弧长、SUPG、log-weighted 输运为可选路径，不能据“已实现”推定全部场景通过。
- 有效命令及修改要求见 [AGENTS.md](AGENTS.md)；本页仅维护状态和证据范围。

## 3. Model contract summary

- SI，三维 HEX20，每节点 `u1,u2,u3,phi,eta_n,eta_p,theta` 七自由度，单元 140 DOF。
- 稳态力学/压电/热传导/双载流子漂移—扩散耦合；LQ/N0/NT 已存在。非线性浓度为 `n0*exp(eta_n)`、`p0*exp(eta_p)`；复合源 `R=0`。
- **默认热源**为漂移电导项；`heat_source_current_work=.true.` 使用 `E·(Jn+Jp)` 总电流功诊断热源，包含扩散电流功，改变热方程。两条路线的结果不能互相替代；J·E 路线不等于完整热电能量闭合。
- SUPG（含方向 tau）、log-weighted、信赖域、伪弧长及 J·E 开关在模型类型中的默认值为 false；实际计算路径必须读输入及 metrics。
- 未见高场迁移率、产生—复合/陷阱/施主离化、时间积分的完整实现；不得将其影响视为已包含。ZnO 材料仍标识为数值候选，非同一试样完整实验标定。
- 计划文档 `docs/model_contract_v2.md` 当前不存在；详细模型契约仍需从源码整理，尤其需明确新增实验分支的适用范围。

## 4. Verification status

### 已通过：限定于对应构建、输入和验收层级

| 证据 | 已确认状态 | 范围限制 |
|---|---|---|
| [R26c 冻结 CLI 对比](runs/HV_SUPG_20261007/regression_r26c_trust/report.json) | 16 项全部 passed：V01–V06、G02–G08、L01、F1M、A1J | 历史候选构建对比，不是最新 ARC/J·E 二进制的全量回归 |
| [JE/ARC 回归](runs/JE_ARC_REGRESSION_20261011_FINAL/summary.txt)及同目录原始日志 | G02、G02JE、G04、G04JE、G07、G08、ARC 均 exit=0 且原始输出 PASS | 保留该次运行身份；不能仅凭目录 FINAL 推定最新构建全覆盖 |
| [修复核心验证记录](runs/20260906_REPAIRED_CORE_0132/repaired_validation_report.md)、[G12](runs/20260905T173523Z_G12_REPAIRED_0132)、[G13](runs/20260905T173149Z_G13_ENDPOINTS) | G09–G13 等存在历史验证证据 | 未在本次对最新构建重跑；G13 为脚本验证，G14 为几何候选/负控制，不能列作通用物理 PASS |
| [A1](runs/A1_100NM_010NM_NE20_DT10K_20260907/A1_stdout.txt)、[A2](runs/F1_20260906/A2/A2_stdout.txt) | 历史解析验证 PASS | 仅支持各自冻结条件，不能替代高压生产验收 |
| [R25 原模型证据登记](runs/HV_SUPG_20261007/delivery_r25_stage4/evidence.json) | 0 V 精确平衡已验证；0.5、1、3、6 V 为 SOLVE_FLUX_PASS_MESH_PENDING | 非零点未获完整空间认证；方向网格、积分和稳定化研究不能合并虚构联合 PASS |

### 最新 J·E 五点：求解/审计通过，空间验收待完成

[交付摘要及输入索引](runs/JE_6_8_HALF_V_20261013/SUMMARY_20261013.md)。同一 XY0_Z2 网格，13,920 节点、3,136 HEX20；方向 SUPG，逐点固定电压重启，伪弧长关闭；五点最终 voltage/Joule 因子均为 1。

| 电压 | 原始 metrics | 独立热平衡相对误差 | 当前状态 |
|---|---|---:|---|
| 6.0 V | [V060](runs/JE_6_8_HALF_V_20261013/V060/metrics.txt) | 0.00431844% | PASS_NOT_MESH_ACCEPTANCE |
| 6.5 V | [V065](runs/JE_6_8_HALF_V_20261013/V065/metrics.txt) | 0.00583595% | PASS_NOT_MESH_ACCEPTANCE |
| 7.0 V | [V070](runs/JE_6_8_HALF_V_20261013/V070/metrics.txt) | 0.00768513% | PASS_NOT_MESH_ACCEPTANCE |
| 7.5 V | [V075](runs/JE_6_8_HALF_V_20261013/V075/metrics.txt) | 0.00989922% | PASS_NOT_MESH_ACCEPTANCE |
| 8.0 V | [V080](runs/JE_6_8_HALF_V_20261013/V080/metrics.txt) | 0.0125105% | PASS_NOT_MESH_ACCEPTANCE |

8 V 最大残差约 7.15e-16、最大更新约 2.69e-10、温升峰值约 24.3734 K。最终场文件存在。以上是 **J·E 诊断热源方程**的候选场，不证明原漂移电导热源模型通过 8 V，也不证明物理热闭合或网格独立性。目录名 HALF_V 指半伏间隔，不据此判定半域。

### FAIL / incomplete

- 原热源模型 [R10 伪弧长 metrics](runs/ARC_FOLD_XY0_Z2_20261011/ARC_FOLD_XY0_Z2_20261011_R10/metrics.txt)：折叠诊断电压约 **6.44734 V**，代数残差小，但独立热平衡误差约 **1.58449% > 1%**，最终 `case_acceptance=FAIL`。折叠括区的弧长尺度不是电压误差界，不能据此宣称物理稳态不存在。
- [同批路线摘要](runs/ARC_FOLD_XY0_Z2_20261011/RESULT_SUMMARY_20261013.md)中的 R12 J·E 弧长路线虽越过 8 V，固定 8 V 回放停滞，未形成该路线的完整验收；它与随后五点固定电压成功记录分别保留。
- R25 登记的原模型 6.5 V 为 SOLVED_FLUX_GATE_FAILED；旧 R25/R26 高压路径也不能作为原模型 0–8 V 全区间通过证据。尚未核实原模型 8 V 的完整接受态。
- [F1 fast audit](runs/F1_FAST_HOLE_L0_20260907/fast_audit.json)：代表性重复的逐字节一致门槛未通过。L0 趋势不等于最终网格验收；F1/F3 生产结论仍需逐项检查证据层级。

### 待确认 / 探索性

- 最新 J·E 五点明确标识 `REQUIRES_INDEPENDENT_MESH_AND_STABILIZATION_STUDIES`；尚无该五点的独立网格/稳定化收敛证书。
- 最新五点摘要声称 ARC、ARCS、F1MJE 及原始回归通过；本次独立定位的回归原始日志覆盖上表七项。ARCS/F1MJE 与最新二进制的完整运行身份链仍需补核，不能只引用摘要视为全量已确认。
- `tests/` 现有 19 个 `test*.py`，不再是空测试目录；上轮定向 `test_hv_runtime.py` 3 项通过。本次未运行全部 Python 测试，也未把历史测试计数套用为当前全量 PASS。
- 历史 JSON 的 RUNNING 只是当时快照；本次检查未见 TPSolverV2 进程，未据旧计划宣称后台求解仍在运行。

### Superseded / 禁作当前验收依据

- [旧 final_results 说明](final_results/README.md)涉及旧 graded-hole 包的 superseded 范围；有边界迹负 Jacobian 的旧结果不得用于当前定量物理结论，不能因其他新算例通过而恢复有效。
- 旧标准 Galerkin 高压失败仍可用于失败分析，不代表当前所有稳定化/J·E 路线失败。
- [论文包冻结通知](paper_sprint/20261008_2110/audit/FROZEN_VERSION_NOTICE.md)明确审计针对 0.15.9 / R25_CLEAN，不能当作当前 0.16.0 ARC/J·E 的源码或全量验收背书。

## 5. Current main blocker and risks

主要阻塞已从“尚未实现高 Pe 稳定化”转为：**高 Pe 空间精度及稳定化验收、热源模型适用范围、各分支证据身份闭环**。

- SUPG/log-weighted、高 Pe 验证和低 Pe 回归入口已存在；实现或 Newton 收敛不保证空间正确。J·E 五点尚缺空间认证，原模型折叠诊断还有独立热平衡失败。
- 原热源与 J·E 诊断热源控制方程不同；跨分支比较必须标注，不能用一条分支补另一条的 FAIL。
- ZnO 参数证据等级、缺失高场/复合物理限制生产解释；实验分支仍需针对性回归。
- 同一版本文字覆盖多个实现阶段、历史记录含另一工作目录绝对路径、后续源码尚未提交：必须绑定源码/构建/输入/网格/结果哈希，不假定路径相同即身份相同。

## 6. Immediate next tasks / 文档同步

1. **P0：证据与基线。** 在已有初始提交上建立后续源码与验证证据的可追踪基线；补齐最新二进制的回归证据映射。本次仅发布文档及精选小型证据，不代表本地所有源码已上传。
2. **P0：同步文档。** 本地执行计划已增加 10 月 21 日受控计算记录，超过本页数值核查截止日；本次未重新审计该批完整证据。后续需同步本页、README、模型契约与验证矩阵，保留历史文件原样。
3. **P1：空间与物理验收。** 分别制定原热源和 J·E 路线的空间/稳定化、独立通量与热源适用性验证；先明确失败指标和证据缺口，再安排新计算或算法修改。
4. **P2：生产结论。** 在对应模型验收闭环后更新 F2 高压结论，再推进 F1/F3 网格与材料证据；暂不外推全区间可信度。

本页替换旧的累计状态日志；更新前全文已备份至仓库根目录 `.codex_tmp/project_state_20261013_233205/PROJECT_STATE.md.bak`。仅本状态页更新，源码、构建、历史输入/结果及其他文档未改写。

## 7. GitHub publication scope — 2026-10-21

本次文档分支上传 AGENTS.md、PROJECT_STATE.md 及状态页直接引用的精选小型证据。源代码仍为该分支继承的初始基线；本页描述的本地实验实现不意味着已包含于该分支源码。大场文件、可执行文件、完整历史运行目录和本地备份未随本次上传；指向这些资产的相对链接须在完整本地工作区使用。此发布不是数值版本发布或完整可复现归档。
