# 已有 S1/S2 验证证据：可用于论文补充材料

本次仅整理已有原生验证日志，未运行求解器。所有图点来自 `attempt01/A1JJE_targeted.txt` 或 `attempt01/JEGRAD_targeted.txt` 的真实记录；脚本在绘图前后重新调用 `validate_acceptance` 校验冻结日志、核心源文件和构建身份。候选程序为 `JE_VERIFY_20261018_FINAL`，其 SHA-256 见 `provenance.json`。后来 A2P/driver 的修改不被冒充成这个已执行程序的源码状态。

## 数据与读图规则

- S1 完整保留 13 个已求解状态，另将 flux patch、双初值、延续和负控制列入附表。S2 完整保留 10 个已求解状态，另列局部 NT 本构检查。没有剔除任何已求解算例。
- 日志没有这两项验证的完整 FEM 节点温度场，所以不绘制所谓“FEM 温度空间曲线”。峰温是程序实际输出的数值峰温；分析解峰温另用叉形标记。标记的横向小错位只是为了区分重合点，不表示不同网格位置。
- 表内保留全部日志指标。`source_metrics_with_lines.csv` 保存原始字符串、解析数值、日志路径、原始行号和 SHA-256，可逐点追溯。
- 图上的重复网格点是确定性计算，不是统计重复；不使用误差条、置信区间或拟合收敛阶。线性电势/浓度变量和二次温升可由 HEX20 精确表示，误差主要处于浮点及迭代误差量级，不能据此宣称复杂高 Pe 场的网格精度。
- S1 零电压状态的最大场误差和通量误差恰为 0，在对数图中不画为一个任意正数；真实零值保留在表内。全部对数图遗漏清单见 `summary.json`。
- S2 日志未输出逐侧边界热流的原始值，CSV 中标为 `NOT_LOGGED`。其 `max_flux_error` 确实包含独立边界热流与热平衡检查，但不能从最大误差反推出各侧原始热流。

## S1 模型与能量表

S1 是旧 A1J 的生产装配验证接口开启 `current_work` 后的独立 ID A1JJE，不改写原 A1 或原 A1J。几何为 20 × 20 × 1 µm，电势为 `phi=-V*x/L`，载流子浓度恒定、迁移率恒定，无机电耦合，端点温升为零。独立解析式为

`sigma=q*(n0*mu_n+p0*mu_p)`，`theta=lambda*sigma*V^2*x*(L-x)/(2*k*L^2)`。

这里必须准确区分：实际 S1 源码将迁移率温度指数设为零，但 Einstein 扩散系数仍随温度变化；由于浓度梯度为零，扩散电流不参与解析解。不能将其误称为显式冻结了所有输运系数。S2 才明确开启 `frozen_transport`。

P1 为合成双载流子材料：n0=1e21 m⁻³、p0=8e20 m⁻³、mu_n=0.05、mu_p=0.02 m²/(V s)、k=10 W/(m K)。P2 为电子占优的尺度检查：n0=6e22 m⁻³、p0=(2.9e-4)²/n0、mu_n=0.0205、mu_p=0.005 m²/(V s)、k=60 W/(m K)。两者均不是对实际 ZnO 材料参数的完整验证。

以下选择协议中固定的 nx=2 两个 profile 与三个控制说明量级；全部 13 例仍在图与 CSV 中。数值按表内位数舍入，完整双精度记录见 CSV。

| 状态 | 电压 V | 热源系数 lambda | 数值峰温 K | 热源积分 µW | 左右端面总热流 µW | 原始端口功率 µW |
|---|---:|---:|---:|---:|---:|---:|
| P1/nx2 | 3 | 1 | 1.189616151 | 95.16929206 | 95.16929206 | 95.16929206 |
| P2/nx2 | 3 | 1 | 3.695019862 | 1773.609534 | 1773.609534 | 1773.609534 |
| 零电压 | 0 | 1 | 0 | 0 | 0 | 0 |
| 半电压 | 1.5 | 1 | 0.2974040377 | 23.79232301 | 23.79232301 | 23.79232301 |
| 关闭热源 | 3 | 0 | 0 | 0 | 0 | 95.16929206 |

图 S1b 用 `lambda × 原始端口功率` 与热源/热流比较，明确显示热源关闭控制的口径。CSV 同时保留未乘 lambda 的原始端口功率，不能称关闭热源时全部电功率仍转成计算中的热。

S1 温度 L2 误差按旧协议使用 `max(解析峰温,1 K)` 归一化，电势使用 `max(abs(V),1 V)`；热/电功率误差按 `sigma*max(V²,1)*面积/L` 归一化。图中特意写 recorded/scaled error，而不将其全部误称为相对误差。

## S2 模型与能量表

S2 是合成线性载流子 LQ 模型，明确开启 frozen transport；n0=p0=N=1e21 m⁻³，mu_n=0.05、mu_p=0.02 m²/(V s)，k=10 W/(m K)，T0=300 K，V=0.02 V。`eta_n=eta_p=s*x/L` 是一阶相对浓度变量，实际 LQ 浓度为 `N*(1+eta)`，不是非线性指数浓度。取 s=-0.05、0、+0.05，三个轴向网格 nx=1/2/4，加一个正梯度 2×2×2 横向检查。

令 E=V/L、VT=kB*T0/q，独立解析电流与热源为

`Jn=q*N*mu_n*(E+VT*s/L)`，`Jp=q*N*mu_p*(E-VT*s/L)`，

`QJE=q*N*((mu_n+mu_p)*E²+(mu_n-mu_p)*VT*E*s/L)`，

`theta=QJE*x*(L-x)/(2*k)`。

下表按固定的 nx=2 状态给出实际数值积分和峰温，按显示精度舍入；三个网格及横向检查全部在图和表中。

| s | JE 体积分 nW | 端口功率 nW | 带符号扩散电流功率 nW | 数值峰温 µK |
|---:|---:|---:|---:|---:|
| -0.05 | 4.361836165 | 4.361836165 | -0.1242584100 | 54.52295207 |
| 0 | 4.486094575 | 4.486094575 | 约 0，完整舍入误差见 CSV | 56.07618219 |
| +0.05 | 4.610352985 | 4.610352985 | +0.1242584100 | 57.62941231 |

这里扩散功率相对解析漂移功率为 ±2.76985712%，正负由浓度梯度决定。数值扩散功率的独立分解为“JE 积分减去按 LQ 参考浓度 N 积分的漂移项”，不是错误套用非线性实际浓度的拆分。温度误差按实际解析峰温归一化，不设 1 K 地板；功率误差按解析漂移功率归一化。

## 冻结误差门限与实际最大值

| 范围/量 | 所有相关状态的最大值 | 原冻结门限 |
|---|---:|---:|
| S1 场误差 | 2.9286453e-13 | 1e-8 |
| S1 独立通量误差 | 2.3993467e-13 | 1e-8 |
| S1 热平衡差 | 3.1665184e-14 | 1e-8 |
| S1 4/5 阶积分差 | 6.9261964e-14 | 1e-10 |
| S2 场/点电流/点热源/通量/扩散功率最大误差 | 2.8699356e-13 | 1e-8 |
| S2 4/5 阶积分差 | 2.5528353e-13 | 1e-10 |

`supplemental_controls_and_local_NT.csv` 另保留 S1 双初值/延续差约 1e-13、10% 热场扰动检出约 0.1、故意破坏电流检出约 0.4；S2 用旧漂移热源替代带符号 JE 后检出 0.0276986 的差异，温度人为乘 1.1 后检出约 0.102770 的差异。这些是验证检查能够拒绝错误状态的负控制，不是真实生产场。

局部 NT 抵消检查中，旧漂移热密度为 2.3496738e9 W/m³，而 JE 为 9.9290728e-8 W/m³（浮点抵消残量），电流抵消缩放误差 1.0902882e-16；独立解析切线误差约 1.9e-16 至 2.8e-16。这些只证明材料点公式/导数在该测试状态的实现一致，不是全局 Poisson 平衡或完整热电闭合证明。

## 英文图注草稿

**Supplementary Fig. S1 | Constant-density Joule verification of the current-work branch.** All 13 accepted A1JJE solves are retained. (a) Logged numerical and analytical peak temperature rises. (b) Independently integrated heat source, outward heat through the two end faces, and the terminal power multiplied by the prescribed heat-source factor lambda. The unscaled terminal power is retained in the source table, including the heat-source-off control. (c) Maximum recorded field error and independent flux error. (d) Difference between independent fourth- and fifth-order integration. P1/P2 denote the two synthetic material profiles; the number following the slash is the longitudinal element count. Rev, zero, half, off and 3D denote reverse voltage, zero voltage, half voltage, heat-source-off and transverse refinement. Exact-zero errors in the zero-voltage case are retained in the table but omitted from logarithmic axes. Horizontal marker offsets distinguish coincident values. These polynomial fields are exactly representable by HEX20; no spatial FEM curves, statistical intervals or asymptotic convergence rates are inferred.

**Supplementary Fig. S2 | Signed diffusion-work verification in a synthetic frozen-transport LQ model.** The three prescribed relative-concentration slopes s=-0.05, 0 and +0.05 are tested with one, two and four longitudinal elements, followed by a positive-slope 2×2×2 refinement. All ten accepted solves are shown. (a) Independently decomposed numerical diffusion-work power and its analytical value. (b) Logged numerical and analytical peak temperature rise. (c) Actual JE volume integral and terminal power. (d) Maximum recorded field/current/heat-source/flux/diffusion-work error and the fourth-/fifth-order integration difference, with the unchanged acceptance limits. The sign of the diffusion correction reverses with s. Horizontal offsets only separate coincident markers. Temperature errors use the actual analytical peak; power errors use the analytical drift power. This is a synthetic implementation benchmark, not a full nonlinear ZnO, high-Pe mesh or thermoelectric-closure validation.

## 可支持和不能支持的结论

两图支持：新 JE 分支在零扩散极限恢复已知热解；非零扩散项的正、负、零情况与独立解析式一致；完整生产 HEX20 装配/Newton 与独立通量路径能达到原冻结误差门限。结合已有 G02JE/F1MJE/G04JE 和旧 CLI 回归，可作为实现正确性的补充证据。

两图不能替代：6/8 V 的联合网格和 SUPG 验收、局部载流子收敛、材料参数校准、真实实验对照、完整热电能量模型。不能据此把原漂移模型结果与 JE 结果混用，也不能把被排除的 A2P 原 FAIL 改判通过。原 A1/A2 及所有旧数据原样保留。

## 图件质量检查

最终 `release_v2` 两图均通过渲染后面板对齐、PDF 最小文字尺寸和 PDF 碰撞检查；最小实测文字为 7.5 pt，碰撞检查均为 0 FAIL、0 WARN。已逐面板查看 PNG，首轮图例压线与标签相碰已在此版本修复。PDF/SVG 保持可编辑文字，PNG 为 600 dpi，未对图像做人工平滑或像素修补。

源代码预检的两条 WARN 已人工复核：一是导出为 PDF/SVG/PNG，未提供尚未要求的 TIFF；二是静态扫描器未识别 `values > 0` 的显式对数过滤。实际程序拒绝负误差、保留表中真实零值并列出零点遗漏，不添加伪计数或修改科学数据。图中没有统计误差条，因为没有统计重复样本。

`provenance.json` 保存脚本、协议、验收、日志、源模块与 QA 工具的 SHA-256；`generated_file_hashes.json` 是脚本运行结束时的图表清单，后加入的说明与人工 QA 则由上层最终交付清单记录。
