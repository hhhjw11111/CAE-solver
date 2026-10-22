# 当前源码模型契约（只读审计，2026-10-08）

> 冻结时刻限定：本文“当前”指本轮启动审计时与R25_CLEAN匹配的0.15.9-log-inverse-checked源码及二进制。21:53交付核查发现工作树另有8个源码文件发生并行变化，版本改为0.16.0-trust-region，具体前后哈希见`../concurrent_source_changes.json`；该变化非本审计写入。下文59项匹配和小型核验结论均为冻结时刻记录，不适用于这些后续变化，也不认证0.16.0。保留当时哈希，不把已完成验证追溯归属给新版本。

本文件依据 `solver_v2/src/` 当前文件编写。源码 SHA-256、大小、构建逐文件对应和本轮运行命令见 `source_audit_runtime_retry.json`；文末源索引给出真实行号。当前 58 个 Fortran 文件及 R25_CLEAN 可执行文件共 59 项与其构建清单逐项哈希一致。本轮未修改生产源码、输入、材料或历史运行；未启动三维计算。

## 变量、单位和材料

每节点依次为 u1,u2,u3,phi,eta_n,eta_p,theta，HEX20 每单元 140 自由度；SI 单位。T=300 K+theta，E=-grad(phi)，小应变 Voigt 顺序 xx,yy,zz,yz,zx,xy，剪切分量是工程剪应变。ZnO 候选材料直接以 c||+z 编码，本条路径没有晶向旋转输入；不可将一般旋转张量功能写入本算例。

| 量 | 当前值及单位 |
|---|---|
| C11=C22；C33 | 209.7；210.9 GPa |
| C12；C13=C23 | 121.1；105.1 GPa |
| C44=C55；C66 | 42.47；44.3 GPa |
| e15=e24；e31=e32；e33 | −0.37；−0.62；0.96 C/m² |
| kappa11=kappa22；kappa33 | 8.54 epsilon0；10.2 epsilon0，epsilon0=8.8541878128×10^-12 F/m |
| k | 60 I W/(m K) |
| alpha | (4.31,4.31,2.49,0,0,0)×10^-6 K^-1 |
| pi | (0,0,−6.9×10^-6) C/(m² K) |
| mu_n0；mu_p0 | 0.0205 I；0.005 I m²/(V s) |
| beta_n=beta_p | 2 |
| n0；ni | 6×10^22；2.9×10^-4 m^-3 |
| p0；ND | ni²/n0=1.4016666667×10^-30 m^-3；n0−p0 |

弹性矩阵的对称项由源码显式赋值，其余项由默认材料类型置零。b=C alpha，不单独引入一个热应力系数。以上是混合来源的有效候选输入，不是同一样品的实验标定。极小 p0 使“程序含双载流子自由度”与“此 ZnO 算例具有显著双极导电”成为不同命题。

## 本构、热源和模型区分

sigma=C epsilon−e^T E−b theta；D=e epsilon+kappa E+pi theta。rho=q(p−n+ND−NA)。N0/NT 用 n=n0 exp(eta_n)、p=p0 exp(eta_p)：Jn=q n(mu_n E+D_n grad eta_n)，Jp=q p(mu_p E−D_p grad eta_p)。热流 q_h=−k grad T。

Qd=heat_factor q[n E·mu_n E+p E·mu_p E]。这是漂移电导型热源，既不是 J_total·E，也不是完整的漂移—扩散能量闭合。当前复合 R及其导数严格为零。没有陷阱、不完全电离、高场迁移率、显式热扩散载流子项、时间积分或挠曲电。

N0 同时冻结 mu=mu0 和 D=kB T0 mu0/q；四个温度导数全零。NT 用 mu=mu0(T/T0)^(-beta)，D=kB T mu/q，一致导数为 mu_T=−beta mu/T、D_T=(1−beta)D/T。因而 NT−N0 是恢复整个温度相关输运的差值，不能称为纯迁移率贡献。

LQ 中 eta 是一阶相对浓度而非对数：n=n0(1+eta_n)、p=p0(1+eta_p)，非正密度会被拒绝。电流使用参考密度，Jn=q n0(mu_n0 E+D_n0 grad eta_n)，Jp=q p0(mu_p0 E−D_p0 grad eta_p)，去掉 eta E 和 eta grad eta 乘积；Qd 也使用 n0/p0。三阶段依次为：（1）theta 固定且热源为零，求其余场；（2）固定前六场，仅解温度；（3）固定 phi/eta_n/eta_p/theta，仅解位移。最终状态没有返回电—载流子反馈，不能冒称与 NT 同义的全场同时平衡。N0−LQ同时混合非线性密度和反馈路径差异。

## 强式、弱式与一致切线

无体力及 R=0 时，div sigma=0，div D=rho，div Jn=div Jp=0，div q_h=Qd。普通测试空间的单元残差按源码符号为

R_u=∫B^T sigma dV；R_phi=∫(grad N·D+N rho)dV；R_n=∫(grad N·Jn+qNR)dV；R_p=∫(grad N·Jp−qNR)dV；R_T=∫(−grad N·q_h−NQd)dV。

状态与测试算子分开构造。Jacobian 由材料一致导数与状态算子相乘装配；热源对E、eta及T导数均在材料点实现。以上残差符号与给定强式相容，不用任意调符号解释求解器输出。

当前另有默认关闭的 log-weighted 输运选项。其测试函数为 w_i=exp(−eta)N_i，载流子弱式项变成 exp(−eta)[(grad N_i−N_i grad eta)·J+s q N_i R]，电子 s=+1、空穴 s=−1。额外切线包含测试函数对 eta 的导数。原物理电流、力学、静电、热传导和 Qd 定义不变；但有限维测试空间改变，离散解及代数守恒含义会变，必须独立重做物理通量和空间检查。该选项与 SUPG、LQ 互斥，可选 3³ 或 5³ 积分；默认普通 Galerkin 为 3³。不能把加权自由行和称为原物理电流守恒。

## 全域与半域边界分别登记

全域包围盒为 [−10,10]×[−10,10]×[−0.5,0.5] μm。左 x=−10 μm 的 phi=0，右 x=10 μm 的 phi=输入 voltage（本研究负号）；左右 eta_n=eta_p=theta=0。A=(−10,−10,−0.5) μm：u1=u2=u3=0；B=(−10,+10,−0.5) μm：u1=u3=0；C=(−10,−10,+0.5) μm：u1=0。合计六个机械约束消除刚体自由度。其他表面与孔壁采用零自然通量/牵引，未求解孔内空气场、表面电荷或接触热阻。

y半域是单独的对称分支模型：y≥0，y=0 整面 uy=0；A=(−10,0,−0.5) μm 固定 ux、uz，B=(+10,0,−0.5) μm 固定 uz。不能沿用全域A/B/C再加对称约束。半域积分是真实半域量，只有明示的重构总量才能乘2；其峰值不能乘2。新半域状态不能默认为与旧全域节点位移参考一致。

## 既有验收及可写范围

`run-f1` 条件为最大阶段balance<1e−8、残差和增量各<1e−9、独立热流及电子电流相对缺口<1%、热源4/5阶积分差<0.1%。启用 SUPG/log-weighted 还检查空穴电流<1%及完整电压/热源因子=1；半域再查物理外边界和对称切面。通过字符串明确为 PASS_NOT_MESH_ACCEPTANCE。收敛和通量通过都不是空间、器件或材料适用性证明。

本轮7次成功调用为版本、V01、V02、G02、F1M、G03与既有R22生产HEX20条带重放，原始输出和退出码均保留于 validation/runtime_retry。首次直接调用因未设置 Intel DLL 搜索路径返回3221225781（0xC0000135），原失败日志保留；追加编译器bin到本任务子进程PATH后全部退出0。这是运行环境补全，不是数值阈值修改。没有重编生产程序。

生产条带采用真实材料点与 HEX20 单元残差/Jacobian，自定义一维缩约装配和稠密解法；不是完整CSR/Newton生产入口验证。已有正Pe电子/空穴镜像不能写成同种电子正负E测试。误差是每单元11个点上的最大浓度/电流偏差，不是相对L2积分。零场、负场、独立积分加密稳定性、完整剖面/端面通量文件均未由此得到。附件第7节规格仍为 PARTIAL_REUSED_EVIDENCE，未冒称完成。


## 源码精确引用

以下文件全部含SHA-256，完整表也见 source_index.md。

| 路径与行 | 内容 | SHA-256 |
|---|---|---|
| `src/core/field_layout_m.f90:5–15` | 七场排序 | `0b02e79e04674d322966a1076e64dbb944ffb1634d6d007a7e24e7acde335c96` |
| `src/core/constants_m.f90:6–7` | q与kB | `dc6872db24efbbb10627ef53a0b46fb182a9b61940a66c2497456d99c31c0d15` |
| `src/physics/zno_candidate_m.f90:8–45` | 全部ZnO候选常数与模型开关 | `407519d9ec931d9d29ed1a32d2204f34a74238f1ed0ccb2cbeaf3b2c63f81026` |
| `src/physics/transport_m.f90:40–58` | N0冻结及NT温导数 | `7215021b45f41f4e0cb7654d1bbe1318275581af2eaade621b9d0cd5dee48bd1` |
| `src/physics/thermo_piezo_m.f90:11–24` | 热压电本构 | `0acecab76b6e6468b90107c75d32b69b74a08fc0599e61584a3cf5f71ea0e65d` |
| `src/physics/material_point_m.f90:107–149` | 密度、载流子电流、热源 | `649601bea7dd70704f9eecfe3d4b28de7cc80cb7b33c4ab6a8ddb8568e2ee0aa` |
| `src/physics/material_point_m.f90:151–204` | 一致材料导数及LQ替换 | `649601bea7dd70704f9eecfe3d4b28de7cc80cb7b33c4ab6a8ddb8568e2ee0aa` |
| `src/physics/recombination_m.f90:10–18` | R=0 | `06133b7b899ca05634e56b40afc9fe3f7509d3b15033a340b30d9a88200cf800` |
| `src/fem/element_kernel_m.f90:63–140` | 积分、生产残差/切线和格式入口 | `ce69b2a924e58a91f2edde78bd0568fc3689371ad9a3e65dd4c30f70db0fbda8` |
| `src/fem/element_kernel_m.f90:143–176` | log-weighted测试及其状态导数 | `ce69b2a924e58a91f2edde78bd0568fc3689371ad9a3e65dd4c30f70db0fbda8` |
| `src/fem/element_kernel_m.f90:178–236` | 普通弱式状态和测试算子 | `ce69b2a924e58a91f2edde78bd0568fc3689371ad9a3e65dd4c30f70db0fbda8` |
| `src/app/f1_case_runner_m.f90:157–178` | LQ三阶段 | `47fe98f1b659509ffcf0cfd1aec8141c7d5ea44aa61617a56b3bee0be010a6c6` |
| `src/app/f1_case_runner_m.f90:315–329` | 原验收门槛 | `47fe98f1b659509ffcf0cfd1aec8141c7d5ea44aa61617a56b3bee0be010a6c6` |
| `src/app/f1_case_runner_m.f90:543–595` | 全域电极与A/B/C | `47fe98f1b659509ffcf0cfd1aec8141c7d5ea44aa61617a56b3bee0be010a6c6` |
| `src/app/f1_case_runner_m.f90:598–657` | 半域对称与刚体约束 | `47fe98f1b659509ffcf0cfd1aec8141c7d5ea44aa61617a56b3bee0be010a6c6` |
| `src/app/f1_case_runner_m.f90:659–680` | LQ固定字段 | `47fe98f1b659509ffcf0cfd1aec8141c7d5ea44aa61617a56b3bee0be010a6c6` |
| `src/io/native_input_m.f90:59–62` | 新格式默认关闭 | `ed1b5979bf08fb8a9dd4261a8868e9dd40bbcbed87c3bd68b27ad3e2e72b8451` |
| `src/verification/verification_driver_m.f90:41–76` | 真实verify入口 | `ab204284f55a2924ca8e48f9520cf3e75bae2a7d1ef45c025048070b7b1c8de9` |
| `src/app/main.f90:29–90` | version/verify/run-f1命令分发 | `acec6d3f0674d523205d22b5d0dfeb2cd95b6416d9c8e37945b9d78a52185fc8` |
| `tests/hv_log_weighted_production_bar.f90:89–113` | 生产核条带采样误差定义 | `fbb1472fd7ce1b416bb208bf24ae34f660c8d564acc1861cb16c9c9234a3d189` |
| `tests/hv_log_weighted_production_bar.f90:123–153` | 自定义缩约装配调用生产单元 | `fbb1472fd7ce1b416bb208bf24ae34f660c8d564acc1861cb16c9c9234a3d189` |
| `tests/hv_log_weighted_production_bar.f90:187–213` | 实际Pe/网格及容差 | `fbb1472fd7ce1b416bb208bf24ae34f660c8d564acc1861cb16c9c9234a3d189` |
