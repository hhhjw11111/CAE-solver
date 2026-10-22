# 3 V 零新增求解后处理报告

本轮重构10个已完成状态：6个F1全域L0（H0/H1×LQ/N0/NT）、3个R23半域NT网格层级及1个R25 A64半域NT。没有调用求解器，没有增加电压、几何或物理。F1的同模型几何对照与同几何模型对照均为历史粗网格探索，未完成空间验收。R23/R25采用试验LOG_WEIGHTED格式，不能当作与F1相同离散格式的复算。

## 证据与重构

input、mesh、nodes、field_samples、metrics及exe逐文件SHA-256保存在audit/postprocess/各case.json；源码快照逐文件哈希与对应归档manifest复核，详见provenance_gate_recheck.json。原始结果所在runs是工作区链接，审计保存resolve后的实际路径；原文件仅只读。10例原exit code为0，残差/增量/阶段平衡/独立通量/4–5阶热源积分按原阈值重新数值核查。本次归档身份检查不等价于从源码重新编译。

F1早期H0三例缺少当时逐case provenance.json；metadata_correction.json中的环境是后补回溯，未补造实时记录。因此这些case的完整运行来源状态为PENDING，归档输入/网格/二进制及54个源码文件哈希匹配为已完成子项。F1 H1有当时provenance；R23/R25各58个归档源码文件匹配。没有缺少两份源码的情况。

全体第一5阶采样点电场/热源/温度/浓度/位移/迁移率/应力重构最大尺度误差4.72e-09；独立5×5×5积分漂移热相对差最大2.81e-13。形函数先插值η再指数；LQ使用1+η和参考密度电流/热源，N0冻结μ与D，NT使用当地T。

几何通过5阶体点及每单元六面5×5点正Jacobian抽样；这是一项必要的局部映射检查，不是全域正性数学证明。废弃旧graded-hole不在输入清单。厚度线与孔周线对物理取样点求逆等参映射，不作域外外推；中面图直接在原始HEX20面片自然坐标上求值。

## 功率与参考偏离（F1 L0，探索性，空间待验收）

| Case | Pd / mW | Pdiff / mW | ∫abs(Jdiff·E) / mW | Pterm / mW | Qout / mW | 采样峰温升 / K | rT RMS / K |
|---|---:|---:|---:|---:|---:|---:|---:|
|FAST_H0_LQ_L0|1.77360953|-0.00000000|0.00000000|1.77360953|1.77360953|3.69318528|0.00000000|
|FAST_H0_N0_L0|1.88626833|-0.11265909|0.11268516|1.77359718|1.87222467|4.02514966|0.22331740|
|FAST_H0_NT_L0|1.85021113|-0.10719967|0.10722486|1.74299983|1.83695267|3.96088472|0.18221670|
|FAST_H1_LQ_L0|1.53951235|-0.00000000|0.00000000|1.53950941|1.53951346|3.69463035|0.00663098|
|FAST_H1_N0_L0|1.62465325|-0.08514112|0.08515643|1.53949940|1.61405416|3.96721008|0.16688930|
|FAST_H1_NT_L0|1.59447271|-0.08118642|0.08120124|1.51327119|1.58444010|3.90687795|0.13321748|

H1 NT相对N0的采样峰温升变化-1.520770%，漂移功率变化-1.857660%；同为NT的H1相对H0采样峰温升变化-1.363503%。这些是新核对原始场后的历史粗网格趋势，不进入正式摘要的定量结论。H1 NT的绝对扩散做功/Pd=5.092671%，提示不能忽略其物理解释；功率比不等于温度误差，需独立能量模型论证，不能用约化漂移热源直接证明真实器件的精细热反馈机制。

常系数参考θref=σ0(φ−φL)(φR−φ)/(2k)，其3V理论峰值为3.695019862 K。rT以5阶体积权重计算；该偏离包含模型与离散差异，不是纯半导体贡献。原metrics中的theta_max是保存的积分采样最大值（不等价连续场严格最大值）。

半域功率保留真实半域值，表中字段分别列全部边界和物理外边界，不自动乘2。P_JdotE与Pterm的差单独保存；未将它强行设为零。漂移/扩散重构没有修改热方程。

## Pe与空间证据

Pe使用20节点电势范围，非二次场内部精确极差。N0使用固定T0对应D/μ；NT使用节点平均T；LQ的Pe仅为电势背景指标。区域的最终定义为：接触区优先选择接触x=±10μm的第一层单元；其余按中心|z|≥0.45μm为表层，再按r≤3.5μm为孔周，其余为体内粗区。上下表层合并；每区计数与体积比例分别输出。初始|x中心|≥9μm出现空区，已在region_definition_correction.json记录诊断修正；全域Pe和物理指标没有变化。

| Case | Pe max | 元素计数Pe>1 | 体积权重Pe>1 | 空间状态 |
|---|---:|---:|---:|---|
|FAST_H0_LQ_L0|13.78037|28.571%|95.000%|EXPLORATORY_L0|
|FAST_H0_N0_L0|13.88921|28.571%|95.000%|EXPLORATORY_L0|
|FAST_H0_NT_L0|13.81822|28.571%|95.000%|EXPLORATORY_L0|
|FAST_H1_LQ_L0|17.74523|82.143%|95.363%|EXPLORATORY_L0|
|FAST_H1_N0_L0|17.85884|82.143%|95.363%|EXPLORATORY_L0|
|FAST_H1_NT_L0|17.66894|82.143%|95.363%|EXPLORATORY_L0|
|r23_log_half_L0_V030|17.66900|82.143%|95.363%|MESH_PENDING_LOCAL_FAIL_3_PERCENT|
|r23_log_half_L1_V030|13.88552|75.000%|94.906%|MESH_PENDING_LOCAL_FAIL_3_PERCENT|
|r23_log_half_L2_V030|11.85285|71.528%|94.875%|MESH_PENDING_LOCAL_FAIL_3_PERCENT|
|r25_A64_V030|8.54514|63.889%|94.261%|SOLVE_FLUX_PASS_MESH_PENDING|

R23三层全局/环带指标原通过不能覆盖局部场失败。L1→L2映射电子与空穴扰动最大差分别3.0065437%与3.0097661%，超过原2%门槛，保留FAIL。Pe>1也是风险诊断，不是任意HEX20的误差界；LOG_WEIGHTED结果不能只凭Pe获得PASS。

## 厚度、力学与图源

厚度线为(x,y)=(0,3.25μm)，距r=3μm孔壁0.25μm；H0使用相同坐标。每条线含401个均匀z点并加入真实网格层面，保存ηn/ηp、n/p相对扰动、φ、Ez、θ及参考偏离。场梯度为所选单元内迹，接缝保留。孔周应力在r=3.25μm、z=0取样并变换成σrr和σθθ，不作节点平滑。厚度变化uz(top)−uz(bottom)保留真实均匀应变；相对位移仅减厚度线中点平移。半域与全域刚体约束不同，原始位移绝对峰值不直接作版本差异。

H1 NT厚度变化10.604187 pm，线上相对位移最大5.304714 pm；孔外中面σθθ绝对采样最大0.533870 MPa。它们均未做独立空间验收，不据此声称ZnO屈服或失效。

图件清单（均有PNG/SVG/PDF）：

- fig01_geometry_boundary：模型/边界/孔外厚度线；示意与结果分开。
- fig02_analytic_errors：从保存A1/A2日志读取真实误差，无新求解；合成解析例不替代生产验证。
- fig03_peclet_regions：同3V不同已存网格的风险分布。
- fig03b_spatial_failure：R23 3V局部载流子3.01% FAIL。
- fig04_midplane_fields：F1 H0/H1 NT中面E/J/θ共同色标，三角片不跨原单元和孔洞。
- fig06_thickness_profiles：F1六例真实材料线。
- fig07_mechanics：未平滑应力与相对位移。
- fig08_power_and_reference：漂移/扩散/端口功率与电热参考偏离。

CSV图源位于data/metrics_summary.csv、peclet_regions.csv、*_element_peclet.csv、*_midplane.csv及独立triangles.csv、profiles/*.csv、A1_error_source.csv与A2_error_source.csv。fig03b原始来源为runs/HV_SUPG_20261007/field_shape_r23_L1_L2_V030.json，已在图件清单记录哈希。

## 阻塞与允许表述

1. F1六例没有完整三层空间证据；R23局部载流子门槛明确FAIL。正式局部场与小比例反馈结论BLOCKED_SPATIAL；图用于定位与方法展示。
2. H0旧环带积分欠采样约28.48%，本轮删除H0/H1环带增强比，不拿色图替代区域积分验证。
3. NT/N0在新R25网格缺少同条件N0状态，不能从旧L0模型差值推断新网格模型差已收敛。
4. 材料参数未实验标定，热源未包含完整漂移–扩散能量闭合；真实器件机制BLOCKED_MODEL_SCOPE。
5. 本轮没有求解独立常导电生产例，解析参考只是解释/诊断工具；LQ不得冒充严格常导电验证。
6. 资源与现有活动任务限制下未新增三维求解，未中断既有高压；未扩新物理/几何/电压。

## 实际执行与改动

运行postprocess_3v.py（退出0；153.77秒）、supplemental_3v.py（退出0）、refine_pe_regions.py（退出0）、本脚本finalize_postprocess.py。全部为新目录脚本；没有修改src、历史输入/结果、exe、分支或标签。首次草稿脚本N0/LQ热电压数组广播报错已修复并完整重跑；该后处理实现错误不改写原算例状态。可重复命令如下：

```text
python -B solver_v2/paper_sprint/20261008_2110/scripts/postprocess_3v.py
python -B solver_v2/paper_sprint/20261008_2110/scripts/supplemental_3v.py
python -B solver_v2/paper_sprint/20261008_2110/scripts/refine_pe_regions.py
python -B solver_v2/paper_sprint/20261008_2110/scripts/finalize_postprocess.py
python -B solver_v2/paper_sprint/20261008_2110/scripts/qa_postprocess.py
```
