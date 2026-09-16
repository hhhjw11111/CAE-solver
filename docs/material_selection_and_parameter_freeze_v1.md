# 真实材料筛选与 ZnO 参数冻结方案 v1

日期：2026-09-04  
状态：**材料定标建议已完成；尚未写入生产材料模块，尚未通过真实材料算例验证**

## 1. 结论先行

第一阶段真实材料建议冻结为：

> **高质量、弱掺杂、n 型纤锌矿 ZnO 单晶或具有单晶输运特征的外延 ZnO；晶体 c 轴与全局 z 轴重合，外加电压沿 x 轴。**

选择 ZnO 的原因不是它一定产生最大的温升，而是它在当前理论所需的四组信息上最完整：

1. 具有明确的 6mm 纤锌矿弹性、压电、介电、热膨胀和热释电参数；
2. 300 K 附近有单晶 Hall 迁移率及温度依赖数据；
3. 弱掺杂 n 型样品可满足非简并统计；
4. 高场输运资料足以给当前低场模型设置清楚的适用性门禁。

这套材料**可以支持**第一阶段的“孔洞导致局部电场、电子浓度、焦耳热、温度和热变形重新分布”研究，但只能有条件地支持“强温度—迁移率反馈”：

- 当前 `600 nm × 600 nm × 20 nm、−5 V` 算例的平均场为 `83.3 kV/cm`，不能继续作为低场 ZnO 的物理主算例；
- 暗态 n-ZnO 中平衡空穴浓度几乎为零，双载流子方程可以保留作统一理论和代码验证，但不能把空穴结果解释为可观测的主要物理机制；
- 若温升达到约 10 K 以上，浅施主的不完全电离带来的 `n(T)` 变化开始与 `mu(T)` 下降竞争。当前只含 `mu(T)`、不含施主电离的模型只能作为局部有效模型。

因此必须继续保留 `TPS_G08_REFERENCE` 作为 G00--G13 数值回归材料，另建真实材料档案 `ZNO_BULK_300K_V1`，不能覆盖原参考材料。

## 2. 候选材料比较

| 候选材料 | 与当前理论的匹配 | 主要优点 | 主要问题 | 决策 |
|---|---|---|---|---|
| n-ZnO 单晶/高质量外延层 | 条件满足 | 压电性强；完整张量和输运文献丰富；300 K 弱掺杂样品可非简并 | 高热导使纳米尺度低场自热很弱；p 型和双极输运不稳定；高场时需场依赖迁移率 | **主材料** |
| n-CdS 单晶 | 部分满足 | 300 K 电子迁移率约 `390 cm^2/(V s)`；体材料热导率约 `16 W/(m K)`，较容易产生温升 | Cd 有毒；陷阱、光电导和声电流饱和明显；270 K 以上热释电/电导会受氧杂质和施主--受主过程影响 | 仅作后续敏感性材料，不作首篇主材料 |
| n-GaN 单晶 | 条件满足 | 参数体系成熟，低掺杂电子迁移率较高 | 300 K 热导率约 `220--269 W/(m K)`，同电压下温升明显弱于 ZnO；更适合高功率器件而非突出热反馈 | 不作为第一选择 |

CdS 的 297 K Hall 迁移率及散射分析来自 Zook 和 Dexter；CdS 在 270 K 以上的氧杂质和施主--受主效应见 Shaldin 和 Matyjasik。GaN 单晶 300 K 热导率的实验值见 Jeżowski 等和 Danilchenko 等。相关来源列于第 9 节。

## 3. `ZNO_BULK_300K_V1` 参数契约

### 3.1 坐标、晶向与 Voigt 顺序

- 全局 `z || [0001]`，即 ZnO 的 c 轴；
- 全局 x、y 位于基面；
- 电压沿 x 方向，因此主电流方向垂直于 c 轴；
- 程序 Voigt 顺序保持为

\[
[xx,yy,zz,yz,zx,xy].
\]

若将晶体 c 轴反向，压电张量和热释电向量的符号必须同时反向，不能只改其中之一。

### 3.2 可直接冻结的 300 K 参数

| 程序参数 | 冻结值 | SI 单位 | 证据与说明 |
|---|---:|---:|---|
| `reference_temperature` | 300 | K | 参考温度 |
| `C11=C22` | 209.7 | GPa | ZnO 单晶超声实验汇编 |
| `C12` | 121.1 | GPa | 同上 |
| `C13=C23` | 105.1 | GPa | 同上 |
| `C33` | 210.9 | GPa | 同上 |
| `C44=C55` | 42.47 | GPa | 同上 |
| `C66=(C11-C12)/2` | 44.3 | GPa | 6mm 对称性计算值 |
| `e31=e32` | -0.62 | C/m2 | 单晶共振--反共振实验值 |
| `e33` | 0.96 | C/m2 | 同上 |
| `e15=e24` | -0.37 | C/m2 | 同上 |
| `eps11=eps22` | `8.54 eps0 = 7.5615e-11` | F/m | 常用静态介电参数 |
| `eps33` | `10.2 eps0 = 9.0313e-11` | F/m | 常用静态介电参数 |
| `alpha1=alpha2` | `4.31e-6` | 1/K | 300 K XRD 热膨胀数据 |
| `alpha3` | `2.49e-6` | 1/K | 同上 |
| `k11=k22=k33` | 60 | W/(m K) | 体 ZnO 300 K 基线；先忽略很小的经典实验各向异性 |
| `pyro3` | `-6.9e-6` | C/(m2 K) | **定应变的一次热释电系数**，与当前本构匹配 |
| `mobility_n0` | `0.0205 I` | m2/(V s) | 300 K 体 ZnO Hall 迁移率 `205 cm2/(V s)`，先作有效漂移迁移率 |
| `beta_n` | 2.0 | 1 | 300 K 附近体 ZnO 高温支的局部幂律基线 |
| `n_ref` | `6.0e22` | 1/m3 | 对应 `6.0e16 cm-3` 的体 ZnO 实测载流子浓度 |

其中 `eps0=8.8541878128e-12 F/m`，`I` 为三阶单位阵。Look 等测得的高质量体 ZnO 在 300 K 时 `n≈6e16 cm-3`、`mu_H≈205 cm2/(V s)`；主施主浓度约 `1e17 cm-3`，受主约 `2e15 cm-3`。Wagner--Helbig 数据的现代汇编指出，约 `1e16 cm-3` 的体 ZnO 在高温支上迁移率近似按 `T^-2` 下降，主要涉及声学形变势和极性光学声子散射。

弹性与压电数据采用单晶实验组，而不是把当前合成参考材料的各向异性数值继续沿用。静态介电常数采用 ZnO 有限元文献常用组；正式投稿前应再对照 Kobiakov 1980 原文表格完成一次逐项核录。

### 3.3 程序中应写入的矩阵

弹性矩阵：

\[
\mathbf C=10^9
\begin{bmatrix}
209.7&121.1&105.1&0&0&0\\
121.1&209.7&105.1&0&0&0\\
105.1&105.1&210.9&0&0&0\\
0&0&0&42.47&0&0\\
0&0&0&0&42.47&0\\
0&0&0&0&0&44.3
\end{bmatrix}\ \mathrm{Pa}.
\]

压电矩阵：

\[
\mathbf e=
\begin{bmatrix}
0&0&0&0&-0.37&0\\
0&0&0&-0.37&0&0\\
-0.62&-0.62&0.96&0&0&0
\end{bmatrix}\ \mathrm{C/m^2}.
\]

介电、热导和电子迁移率张量：

\[
\boldsymbol\epsilon=\operatorname{diag}
(7.5615,7.5615,9.0313)\times10^{-11}\ \mathrm{F/m},
\]

\[
\mathbf k=60\mathbf I\ \mathrm{W/(mK)},\qquad
\boldsymbol\mu_{n0}=0.0205\mathbf I\ \mathrm{m^2/(Vs)}.
\]

### 3.4 热应力向量不能直接等同于热膨胀系数

当前代码使用

\[
\boldsymbol\sigma=\mathbf C\boldsymbol\varepsilon
-\mathbf e^T\mathbf E-\boldsymbol\beta^T\theta,
\qquad
\boldsymbol\beta^T=\mathbf C\boldsymbol\alpha.
\]

代入上面的 ZnO 参数得到程序真正需要的 `thermal_stress`：

\[
\boldsymbol\beta^T=
\begin{bmatrix}
1.687447&1.687447&1.431103&0&0&0
\end{bmatrix}^{T}
\times10^6\ \mathrm{Pa/K}.
\]

不能把 `4.31e-6` 和 `2.49e-6` 直接填入 `thermal_stress`。

### 3.5 为什么热释电必须填一次系数而不是总系数

当前代码同时包含热膨胀和压电耦合：

\[
\mathbf D=\mathbf e\boldsymbol\varepsilon
+\boldsymbol\epsilon\mathbf E+\mathbf p^\varepsilon\theta.
\]

因此 `pyroelectric` 对应的是定应变一次热释电系数 `p^epsilon`。ZnO 的文献汇编给出：

\[
p_3^\varepsilon=-6.9\times10^{-6},\quad
p_3^{\mathrm{secondary}}=-2.5\times10^{-6},\quad
p_3^\sigma=-9.4\times10^{-6}
\quad\mathrm{C/(m^2K)}.
\]

代码应填 `p3=-6.9e-6 C/(m2 K)`。如果填总系数 `-9.4e-6`，同时又让热膨胀应变经过 `e epsilon` 产生极化，就会重复计算二次热释电贡献。

以本报告冻结的 `e` 和 `alpha` 计算：

\[
2e_{31}\alpha_a+e_{33}\alpha_c
\approx-2.95\times10^{-6}\ \mathrm{C/(m^2K)},
\]

与文献二次贡献 `-2.5e-6` 同量级，差别反映不同实验数据组之间的正常离散性。

## 4. 条件性参数：载流子与空穴

### 4.1 当前代码应采用“有效离化掺杂”

Look 等报告的是：

\[
N_D\approx10^{17}\ \mathrm{cm^{-3}},\quad
N_A\approx2\times10^{15}\ \mathrm{cm^{-3}},\quad
n(300\,\mathrm K)\approx6\times10^{16}\ \mathrm{cm^{-3}}.
\]

当前 `compute_equilibrium_densities` 默认施主和受主完全离化。若直接输入上述化学掺杂浓度，程序会得到约 `9.8e16 cm-3`，与实测自由电子浓度不一致。因此在尚未加入施主电离模型之前，真实材料档案应写成：

\[
N_{D,\mathrm{eff}}=6.0\times10^{22}\ \mathrm{m^{-3}},\qquad
N_{A,\mathrm{eff}}=0,
\]

并明确它是 300 K 的**有效净离化掺杂**，不是化学掺杂浓度。

### 4.2 本征浓度和空穴参数

采用 `Eg=3.37 eV`、电子有效质量 `0.24m0` 和代表性空穴质量 `0.59m0` 的非简并估算，300 K 时：

\[
n_i\approx2.9\times10^{-4}\ \mathrm{m^{-3}},\qquad
p_0=\frac{n_i^2}{n_0}\approx1.4\times10^{-30}\ \mathrm{m^{-3}}.
\]

这不是排版错误：宽禁带 n-ZnO 的热平衡空穴确实可忽略。建议的运行参数为：

| 参数 | 建议值 | 状态 |
|---|---:|---|
| `intrinsic_concentration` | `2.9e-4 1/m3` | 由带隙和有效质量换算；数量级可靠，精确值不重要 |
| `mobility_p0` | `0.005 I m2/(V s)` | 仅数值占位，来自约 `50 cm2/(V s)` 的理论预测 |
| `beta_p` | 2.0 | 仅数值占位；没有同一样品的可靠局部拟合 |

由于 `p0` 极小，上述空穴迁移率和指数不会影响电流、焦耳热或泊松方程。任何论文结果都不应讨论该算例中的“空穴集中”。如果研究目标要求电子和空穴均显著，就必须改成光生载流子、双极注入或 p--n 结问题，并加入产生--复合、接触和可能的陷阱模型；当前 `R=0` 的暗态均匀 n-ZnO 不是合适模型。

## 5. 温度--迁移率模型的证据与限制

第一阶段可使用

\[
\mu_n(T)=0.0205\left(\frac{T}{300}\right)^{-2.0}
\ \mathrm{m^2/(Vs)}.
\]

理由：体 ZnO 的高温支在室温附近近似 `T^-2`，且文献将其归因于声学形变势和极性光学声子散射。Makino 等对非简并高质量外延 ZnO 的独立数据给出 `mu_H(300 K)=440 cm2/(V s)`，其图中 200--300 K 数据数字化后得到局部指数约 `2.78`；这说明 `beta` 是样品和温区相关的有效参数，而不是普适常数。

因此冻结：

- 主值：`beta_n=2.0`，对应体 ZnO 基线；
- 敏感性范围：`beta_n=1.8, 2.0, 2.2`；
- 不采用当前合成材料的 `beta_n=1.35` 作为 ZnO 值；
- 正式与具体试样对比时，必须用该试样 280--330 K 的 Hall 或电导数据重新拟合。

Einstein 关系在 300 K 给出：

\[
D_n=\frac{k_BT}{q}\mu_n
=5.300\times10^{-4}\ \mathrm{m^2/s}.
\]

### 5.1 施主电离造成的遗漏

60 meV 浅施主在 300 K 并非完全离化。用简单非简并施主电离模型估算，从 300 K 升到 310 K 时自由电子浓度约增加 3%，而 `beta_n=2` 给出的迁移率约降低 6.3%。于是电导净下降约 3.5%，不是固定载流子模型预测的 6.3%。

这意味着：

- 若只研究 `Delta T <= 5--10 K` 的局部响应，可把当前模型称为“300 K 附近的有效迁移率模型”，并做 `n(T)` 不确定性敏感性分析；
- 若要研究几十 K 的强反馈或温度稳定平台，必须加入施主离化/电中性约束，或直接由同一样品的 `sigma(T)` 和 `n(T)` 分别定标；
- 仅靠 `mu(T)` 下降不能证明温度最终稳定，稳态稳定还来自热传导和边界散热。

## 6. 当前几何与电压是否满足模型

### 6.1 当前 `600 nm、−5 V`：不满足低场物理解释

\[
E_{\mathrm{nom}}=\frac{5}{600\times10^{-9}}
=8.33\times10^6\ \mathrm{V/m}=83.3\ \mathrm{kV/cm}.
\]

ZnO 的脉冲高场实验已在约 `50--250 kV/cm` 观察到偏离欧姆定律的中等高场输运。当前平均场已经进入该区间，而且孔边局部场还会更高。因此：

- `−5 V / 600 nm` 可保留为数值应力测试；
- 不得用它验证低场 `J=q n mu E` 的材料真实性；
- 不能把此条件下的非线性全部归因于温度--迁移率反馈。

### 6.2 推荐的真实材料主算例尺度

在不改变 `−5 V` 的情况下，把几何整体放大到微米尺度：

- 板：`20 um × 20 um × 1 um`；
- 圆柱孔半径：`2 um、3 um、4 um`；
- 左端 `phi=0`，右端 `phi=-5 V`；
- 左、右电极同时定温 `T=300 K`；
- 其余外表面及孔壁电绝缘、绝热、力学自由；
- 左端采用消除刚体位移的最小力学约束，不把整面三个方向全部锁死作为唯一主边界；
- `c || z`，电流沿 x。

这时

\[
E_{\mathrm{nom}}=2.5\ \mathrm{kV/cm},
\]

比已报告的约 50 kV/cm 非欧姆区低一个数量级以上。仍需在计算后检查孔边 `E_95`，建议冻结门禁：

\[
E_{95}\le10\ \mathrm{kV/cm},
\]

并报告超过 `10 kV/cm` 的体积分数，而不是只报告单点最大值。

对于无孔、一维、常电导、两端定温且侧面绝热的理想条带，

\[
\sigma_0=q n_0\mu_{n0}\approx197.1\ \mathrm{S/m},
\]

\[
\Delta T_{\max}=\frac{\sigma_0V^2}{8k}\approx10.3\ \mathrm K,
\qquad
\overline{\Delta T}=\frac{\sigma_0V^2}{12k}\approx6.84\ \mathrm K.
\]

在 `Delta T=10.3 K` 处，`beta_n=2` 预测迁移率下降约 6.5%，足以作为非线性与冻结迁移率线性模型之间的可辨识差异。注意在理想一维定电压问题中，`Delta T_max` 与长度抵消；放大长度主要用于降低场强，而不会消除焦耳温升。

这些数值是前置解析预测，不是三维孔洞结果。孔洞可能同时提高孔边局部焦耳热并降低总焦耳功率，因此孔边温度的方向必须由计算和能量平衡决定，不能预先规定为升高。

## 7. 真实材料算例的物理门禁

运行 `ZNO_BULK_300K_V1` 前后至少检查：

| 门禁 | 判据 | 当前 600 nm / -5 V | 推荐 20 um / -5 V |
|---|---|---:|---:|
| 非简并 | `n/Nc << 1`；本参数约 `0.02` | PASS | PASS |
| 低场 | 名义场 `<=2.5 kV/cm`，且 `E95<=10 kV/cm` | **FAIL** | 预期 PASS，需计算确认 |
| 温区 | 大部分体积 `280--330 K` | 需换材料后重算 | 解析预测 PASS |
| 迁移率可辨识 | 非线性与冻结迁移率结果差异大于数值/网格误差 | 当前合成材料仅约 0.2% 迁移率变化 | 解析预测约 6.5% 峰值变化 |
| 热模型 | 两端定温、其余绝热与解析基准一致 | PASS | PASS |
| 薄膜效应 | 体材料参数可代表试件 | 20 nm 厚度不充分 | 1 um 仍需实验样品核验，但明显更合理 |
| 空穴可观测性 | `p mu_p/(n mu_n+p mu_p)` 非忽略 | **FAIL** | **FAIL** |
| 施主电离 | `n(T)` 变化相对目标效应可忽略或已建模 | 条件性 | 温升约 10 K 时需敏感性分析 |

## 8. 参数不确定性与建议计算矩阵

首轮不应只跑一个参数点。建议在相同网格和边界上做：

1. 主值：`mu_n0=0.0205 m2/(V s)`，`beta_n=2.0`，`k=60 W/(m K)`；
2. 迁移率指数：`beta_n=1.8, 2.0, 2.2`；
3. 热导率：`k=46, 60, 100 W/(m K)`，覆盖不同体 ZnO 表面处理/样品质量报告；
4. 线性基准：固定 `mu=mu(300 K)`、冻结载流子指数线性化，并保持相同几何、边界和焦耳热定义；
5. 非线性主模型：指数载流子 + `mu(T)` + Einstein 扩散；
6. 可选物理修正：加入 `n_eq(T)` 或施主电离后，再与第 5 项比较。

线性基准必须作为独立模型实现。现有 `M1_TWO_CARRIER_ISOTHERMAL` 仍保留指数载流子，不能直接等同于严格线性模型，也不能用“只做一次 Newton”代替线性化。

## 9. 主要文献与参数来源

### ZnO 电输运

1. A. R. Hutson, “Hall Effect Studies of Doped Zinc Oxide Single Crystals,” *Physical Review* 108, 222 (1957). 300 K Hall 迁移率约 `180 cm2/(V s)`，200 K 以上极性光学声子散射重要。<https://doi.org/10.1103/PhysRev.108.222>
2. D. C. Look et al., “Electrical properties of bulk ZnO,” *Solid State Communications* 105, 399--401 (1998). 300 K 时 `n≈6e16 cm-3`、`mu≈205 cm2/(V s)`，并给出主施主和受主浓度。<https://doi.org/10.1016/S0038-1098(97)10145-4>
3. P. Wagner and R. Helbig, “Halleffekt und Anisotropie der Beweglichkeit der Elektronen in ZnO,” *Journal of Physics and Chemistry of Solids* 35, 327--335 (1974). 4--370 K 各向异性 Hall 输运。<https://doi.org/10.1016/S0022-3697(74)80026-0>
4. T. Makino et al., “Electron transport in ZnO thin films,” *Applied Physics Letters* 87, 022101 (2005). 非简并外延 ZnO 的温度依赖输运及散射分解。<https://doi.org/10.1063/1.1991994>
5. L. Ardaravičius et al., “Electron drift velocity in wurtzite ZnO at high electric fields: Experiment and simulation,” *Journal of Applied Physics* 126, 185703 (2019). 报告约 `50--250 kV/cm` 的中等高场非欧姆区。<https://doi.org/10.1063/1.5100078>

### ZnO 弹性、压电、介电和热学

6. I. B. Kobiakov, “Elastic, piezoelectric and dielectric properties of ZnO and CdS single crystals in a wide range of temperatures,” *Solid State Communications* 35, 305--310 (1980). 室温完整张量及温度依赖实验。<https://doi.org/10.1016/0038-1098(80)90502-5>
7. U. Ozgur and H. Morkoc, “General Properties of ZnO,” in *Zinc Oxide: Fundamentals, Materials and Device Technology* (Wiley-VCH, 2009). 汇编弹性、压电、热膨胀、热导和热释电数据，并区分一次、二次和总热释电系数。<https://doi.org/10.1002/9783527623945.ch1>
8. H. Iwanaga, A. Kunishige, and S. Takeuchi, “Anisotropic thermal expansion in wurtzite-type crystals,” *Journal of Materials Science* 35, 2451--2454 (2000). 300 K ZnO 热膨胀系数来源。<https://doi.org/10.1023/A:1004709500331>
9. D. T. Morelli and G. A. Slack, “High Lattice Thermal Conductivity Solids,” 体 ZnO 300 K 实验热导率约 `0.6 W/(cm K)=60 W/(m K)`。<https://djena.engineering.cornell.edu/hws/slack_high_thermal_conductivity_solids.pdf>
10. D. I. Florescu et al., “High spatial resolution thermal conductivity of bulk ZnO (0001),” *Journal of Applied Physics* 91, 890--892 (2002). 不同表面处理的局部热导率约 `46--147 W/(m K)`。<https://doi.org/10.1063/1.1426234>
11. G. Heiland and H. Ibach, “Pyroelectricity of zinc oxide,” *Solid State Communications* 4, 353--356 (1966). 9.7--420 K 热释电测量及 300 K 定应变系数推导。<https://doi.org/10.1016/0038-1098(66)90187-6>
12. “Piezoelectric property comparison of two-dimensional ZnO nanostructures for energy harvesting devices,” *RSC Advances* 11 (2021). 提供常用 ZnO 有限元弹性、压电和静态介电参数组。<https://doi.org/10.1039/D0RA10371C>

### 候选材料排除依据

13. J. D. Zook and R. N. Dexter, “Galvanomagnetic Effects in Cadmium Sulfide,” *Physical Review* 129, 1980 (1963). CdS 的 77--297 K Hall 迁移率及散射机制。<https://doi.org/10.1103/PhysRev.129.1980>
14. Yu. V. Shaldin and S. Matyjasik, “Pyroelectric properties of the wide-gap semiconductor CdS in the low-temperature region,” *Semiconductors* 48, 562--569 (2014). 指出 270 K 以上氧杂质和施主--受主过程的影响。<https://doi.org/10.1134/S1063782614050194>
15. A. Jeżowski et al., “Thermal conductivity of GaN crystals in 4.2--300 K range,” *Solid State Communications* 128, 69--73 (2003). n-GaN 300 K 热导率约 `220 W/(m K)`。<https://doi.org/10.1016/S0038-1098(03)00629-X>
16. B. A. Danilchenko et al., “Thermal conductivity of heavily doped bulk crystals GaN:O,” *Materials Research Express* 2, 085902 (2015). 低氧掺杂样品 300 K 热导率可达 `269 W/(m K)`。<https://doi.org/10.1088/2053-1591/2/8/085902>

## 10. 检索与证据限制

- 检索日期：2026-09-04；关键词覆盖 `ZnO/CdS/GaN mobility temperature`, `Hall`, `thermal conductivity`, `piezoelectric`, `elastic`, `dielectric`, `pyroelectric`, `high field`。
- 优先采用原始实验论文和出版社书章；数值经单位换算后写入本报告。
- 本环境未提供可调用的 Scopus/ScienceDirect 学术 MCP 检索接口，因此使用出版社页面、DOI 元数据、开放预印本和可访问书章交叉核验。Kobiakov 1980 的全文表格仍需在正式投稿前通过馆藏再核录一次。
- `beta_n=2.0` 是室温附近的有效幂律，不是 ZnO 在任意温度、任意缺陷浓度下的普适常数。
- `n_i`、热应力向量、扩散系数和理想一维温升是基于已列输入的本报告换算/推导值，不是原文直接给出的实验数字。
- 本报告由 AI 辅助检索、换算和综合；所有拟用于论文的参数仍应由作者对原始论文表格、试样生长方式和测量方向进行最终人工核对。

## 11. 下一步门禁

在把 `ZNO_BULK_300K_V1` 写入源码之前，先冻结以下实现任务：

1. 新增独立材料档案，不修改 `TPS_G08_REFERENCE`；
2. 为真实 ZnO 的极小 `p0` 做材料点、缩放和线性求解专项测试；
3. 实现严格的冻结系数线性孔洞基准，而不是复用 M1 或只做一次 Newton；
4. 生成 `20 um × 20 um × 1 um` 的无孔和 `R=2,3,4 um` 匹配网格；
5. 先跑无孔一维解析温升对比，再跑三维孔洞；
6. 输出 `E95`、高场体积分数、`Delta T`、迁移率降幅、总功率及局部环带指标；
7. 若主结果依赖 `Delta T>10 K`，在论文物理结论前补充施主电离或 `n(T)` 敏感性模型。

