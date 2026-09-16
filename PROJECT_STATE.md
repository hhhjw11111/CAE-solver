# CAE Solver Project State

## 1. Current baseline

- 基线日期：2026-09-15。
- 当前程序：`solver_v2`；理论版本 `2.1`，架构标识 `TPS-V2-ARCH-001`。
- `src/core/build_info_m.f90`：`0.14.9-cauchy-restart`。
- `VERSION`：`0.14.5-f3-radius-threaded`。两者不一致，当前不擅自修改。
- 当前构建环境：Intel oneAPI ifx、Intel MKL PARDISO；源码构建脚本支持 sequential/threaded MKL 及 IC/OOC PARDISO。
- Git：当前 `main` 无提交，`solver_v2/` 尚未形成正式可追踪基线。

## 2. Current architecture

输入 (`io`) → 网格/HEX20 (`mesh`) → 材料与物理 (`physics`) → 单元 residual/Jacobian (`fem/element_kernel_m`) → 全局 CSR (`fem/global_assembly_m`, `sparse_pattern_m`) → PARDISO (`linear`) → Newton/缩放/Armijo (`nonlinear/newton_solver_m`, `scaling_m`) → voltage/Joule continuation (`continuation_solver_m`) → 独立通量、守恒、场指标和结果写出 (`post`, `app`)。

## 3. Model contract summary

- 每节点 7 DOF：`u1,u2,u3,phi,eta_n,eta_p,theta`；HEX20 局部 140 DOF。
- 当前耦合包括力学、压电、热传导、双载流子漂移—扩散和 Joule 热。
- `LQ`、`N0`、`NT` 均已实现；具体定义以 `src/physics` 和输入契约为准。
- 非线性浓度使用指数映射；当前 Joule 热为漂移电导项；复合为 `R=0`。
- 尚未包含：输运稳定化、高场迁移率、产生—复合/陷阱/施主离化、时间积分等。

## 4. Verification status

### Frozen / PASS

V01–V06、G02–G08、G09–G13、L01、A1、A2、A1J 均有对应 PASS 验证证据；G14 作为几何负控制使用。解析验证、材料点/单元切线、全局装配、Newton、延续和几何/网格审计均有专门入口。

### Exploratory

F1 ZnO H0/H1 × LQ/N0/NT 的 L0 结果支持趋势和探索；ZnO 参数是数值候选输入，不是同一试样完整实验标定。F2 的部分 NT 低电压点有求解/守恒记录，但高 Pe 限制其定量解释。

### FAIL / incomplete

- F2：3.5 V 及以上已出现标准 Galerkin 高 Pe 风险；4–6 V 有热流或完整门禁失败，7–8 V 没有完整 NT 状态。
- F1 的代表性重复运行未满足逐字节一致门槛；原始 L2 最终验收未完整闭环。
- 高 Pe 稳定化、一维高 Pe 基准和低 Pe regression 尚未实现。

### Superseded / forbidden evidence

`final_results/` 中旧 graded-hole 结果因边界迹负 Jacobian 已 superseded。不得将其用于当前物理结论或论文定量结果。

## 5. Current main blocker

当前最主要数值阻塞是高 Péclet 漂移—扩散空间离散稳定性。标准 Galerkin 在 F2 网格上处于漂移占优区；Newton 残差小不能证明载流子空间离散正确。下一阶段应研究有理论依据的 Scharfetter–Gummel、exponential fitting 或其他稳定化，并同步推导 residual/Jacobian，完成一维高 Pe 验证和低 Pe 回归。解决前，高压 F2 定量结论可信度受限。

## 6. Known risks

- 高 Pe 输运离散风险；
- Joule 热仅覆盖漂移电导项；
- ZnO 参数证据等级有限；
- `build_info_m.f90` 与 `VERSION` 标签不一致；
- 旧孔洞结果可能被误读；
- Git 尚无正式版本基线；
- 常规单元测试目录不足；
- F2 高压路径尚未闭环。

## 7. Immediate next tasks

### P0

- 建立可追踪 Git 初始基线；
- 固定模型、验证、构建和 provenance 文档；
- 制定版本号同步规则。

### P1

- 建立一维高 Pe 基准；
- 研究并实现有一致 Jacobian 的输运稳定化；
- 建立低 Pe regression。

### P2

- 在稳定化和回归通过后重新评估 F2 高压状态；
- 继续完善 F1/F3 生产算例的网格与证据闭环。
