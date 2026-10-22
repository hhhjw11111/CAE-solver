# 纯数学参照，不是生产求解器验证

运行脚本：`audit/source_validation/formula_reference_check.py`。本轮实际执行 `python -B solver_v2/paper_sprint/20261008_2110/audit/source_validation/formula_reference_check.py`，进程返回码0。脚本首先保存frozen_contract.json，再计算数值；同名输出目录已存在会拒绝覆盖。

`drift_diffusion_reference.csv`包含A=0、1、10、40、−40的闭式浓度和漂移/扩散/总电流。坐标是x/L，浓度是n/nL，电流尺度是qDnL/L，nR/nL=2。A=muEL/D，不是三维网格单元Pe。正参数使用expm1；负参数使用指数平移避免exp(|A|)溢出；A=0显式使用线性极限。这些是解析参考样本，不是FEM解。

端点误差均为0，恒流关系最大归一化误差3.56×10^-16以内。4096面板Simpson独立积分检查Pd+Pdiff=Pterm=∫J·E，最大归一化缺口1.24×10^-12以内；正负A均按进入材料电功率为正。另以A=±10^-8检查零场极限，电流绝对偏差1.51×10^-8以内，低于事前2×10^-8容差。这些误差全部只属于公式实现检查。

`thermal_reference_3V.csv`是各向同性常欧姆参考在20 μm一维坐标的解析温度，sigma0=197.067725982 S/m、k=60 W/(m K)，3 V峰温3.695019862 K。它不是LQ结果，也不是圆孔三维求解。

`local_equilibrium_heat_counterexample.csv`以无量纲A=1、n=exp(−x)显示漂移与扩散电流严格抵消，总电流做功为零而漂移热仍大于零；这是约化热源适用性的数学反例，不是对热源的代码改动。

本文件夹不补齐附件第7节的生产路径缺项。真实生产核重放位于`../runtime_retry/`，其归档版本、采样误差定义和未完成规格另列。
