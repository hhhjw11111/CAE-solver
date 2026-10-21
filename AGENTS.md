# TPSolverV2

## 场图风格复用

- 2026-10-14 用户指定的参考图为 `output/F3_焦耳热_下排再次放大_20260909/comparisons/compare_Q_midplane.png`（相对于项目根目录）。二维云图与放大三维曲面对照复用 `scripts/styles/joule_map_height_larger_surfaces.json` 的完整布局、字体、视角、配色和导出约定，不能只读取配色后重新写死版式。
- 新工况按实际数据采用共用色标范围，保留热源定义、变量及单位差异；不要套用旧数据的 0–17 数值范围，也不要为复制旧图光滑度而抹去原始单元间导数场跳变。
- 绘图修订写入新目录并保存源数据与风格哈希，历史图件保留。当前可复用入口为 `scripts/plot_je_map_height_20261014.py`。

## 架构与目录

- 三维 HEX20 多物理场 FEM：力学、压电、热传导、双载流子漂移—扩散及 Joule 热；七场联合 Newton，CSR/PARDISO 稀疏直接求解。
- `src/physics` 为本构与输运，`src/fem` 为单元 residual/Jacobian 和全局装配，`src/nonlinear` 为 Newton/延续，`src/linear` 为线性求解，`src/post` 为场量及独立通量审计。
- `cases/` 存输入和网格；`src/verification/`、`tests/` 存原生验证与 Python 测试；`scripts/` 存构建、运行和审计入口。`runs/`、`final_results/` 是结果与证据，`build/`、`x64/` 是构建产物；结果有效范围按对应审计记录判断。

## 模型边界

- SI 单位；节点顺序 `u1,u2,u3,phi,eta_n,eta_p,theta`，HEX20 为 140 DOF。非线性浓度为 `n0*exp(eta_n)`、`p0*exp(eta_p)`；LQ 使用一阶相对浓度，N0 冻结迁移率和扩散系数，NT 包含温度相关输运。
- 默认热源为漂移电导项；可选 `heat_source_current_work` 是完整漂移—扩散电流做功诊断，不等于完整热电能量闭合。当前复合为 `R=0`；解释结果时核对实际模型开关，不宣称包含未实现物理。

## 命令

以下从 `solver_v2/` 执行；`$exe` 指向本次核对过的构建，尖括号参数需替换。依赖 Windows Intel oneAPI ifx、MKL 和 Visual Studio x64 工具链。

- 构建：`./scripts/build_debug.ps1 -BuildSubdirectory 'build/<new-build>'`；优化/线程构建按需加 `-Optimize -Threaded`。
- 身份/输入检查：`& $exe version`；`& $exe inspect <case.nml>`。
- 参考材料运行：`& $exe run <case.nml> <new-result.dat>`；ZnO：`& $exe run-f1 <case.nml> <new-output-directory>`（先创建空输出目录）。
- 原生定向验证：`& $exe verify <ID>`；如 V01、G02、G03、G04、G06、G07、G08、L01、F1M、A1、A2、A1J。G09–G12 需输入路径，参数见 README；实验验证 ID 见 `src/verification/verification_driver_m.f90`。
- Python 定向测试：`python -m unittest discover -s tests -p 'test_hv_runtime.py'`，按改动替换匹配文件；整组 HV 测试用 `test_hv_*.py`。
- 冻结 CLI 对比：`python scripts/run_hv_regression.py --baseline <baseline.exe> --candidate <candidate.exe> --output <new-directory>`。
- `scripts/run_verification.ps1 -ExecutablePath <exe>` 是包含双运行、生产算例和 Visual Studio 重建的广泛验证入口，按影响范围选用。

运行线程构建前确保 oneAPI 编译器运行库目录在 PATH；OOC 参数与隔离 scratch 由 `scripts/hv_runtime.py` 管理。构建来源和版本规则见 `docs/build_and_provenance.md`。

## 核心修改验证

- 修改前明确受影响方程、本构或数值模块，以及 residual/Jacobian 是否变化；修改后编译并运行相关验证，报告残差、守恒及结果变化。
- 按改动选验证：本构/切线对应 G02/F1M、G03/G06；装配/边界对应 G04 和相关解析算例；线性求解对应 L01；Newton/延续对应 G07/G08 及受影响实验分支；调度脚本对应 Python 测试。解析或切线验证不能仅用生产趋势替代。
- 输运稳定化需检查 residual、Jacobian、高 Pe 基准及低 Pe 回归。Newton 收敛不能替代空间离散、独立通量或网格验收；高 Pe 问题不能只靠调整 Newton/线搜索解决。
- 不得隐式放宽容差制造 PASS；区分数值收敛、守恒、空间精度和物理可信度，保留失败记录。
