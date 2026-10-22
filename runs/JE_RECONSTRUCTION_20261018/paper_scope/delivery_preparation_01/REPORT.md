## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: run and evidence reuse
- Origin Date: 2026-10-19
- Verification Status: PARTIAL — completed preparation is verified; mesh study is running
- Version Label: paper_preparation_01

## 本快照已完成的工作

原A1、A2、原漂移热源默认路径、40个生产核心源码及历史输入/结果保留。
归档已在此前完成；本轮不再复制或替换这些历史结果。

| 工作 | 实际结果 | 证据 |
|---|---|---|
| 必要验证复用 | 16项旧CLI逐字节一致；10项专项PASS；源码/构建/日志身份重新核验 | ../verification_reuse/acceptance.json |
| 厚度方向复用 | 6/8 V的Z14/16/18三层审计及4组局部场比较重新执行，与原结果一致 | ../directional_reuse/REPORT.md |
| S1/S2补充图表 | 已有13+10个真实状态；完整峰温、能量、误差和来源表；无新增求解 | ../supplemental_figures/release_v2/README.md |
| 三孔径网格生成 | r=2/3/4 µm几何检查通过；r3三组网格逐字节复现旧生成器 | ../production_preparation/geometry_validation_01/report.json |
| 九工况输入准备 | 28项软件检查及2项实际转移身份检查通过，支持半伏点和三个孔径 | ../production_preparation/input_validation_01/report.json |

旧11项专项中的A2P仍为FAIL，旧整体验证仍为FAIL；本轮用户明确删除A2扩展后建立的是单独的10项论文范围清单，未改写原结果。
S1只固定迁移率，浓度梯度为零，D(T)不参与其解析电流；S2明确冻结输运。两项都不是完整三维高Pe解析证明。

## 首个新试算

FULL A40/R14/Z14、SUPG=1、6 V的34,800节点试算正在
`../../stage4_paper_mesh/A40_R14_Z14_S1_V060/`运行。
旧A32同电压JE场映射到新网格，外推节点数为0；映射只是初值。
此快照时已进入Newton迭代，尚无完成/通量/功率/局部场证书。
具体最新运行状态以该例的stdout.txt、execution.json、summary.json和completion.json为准；缺少完成文件不能视为PASS。

资源为单原生进程、8线程、OOC 6000 MB、单次7200秒。
磁盘和系统可用内存按预先冻结限制监控，仅在本次进程退出后清理该新例可再生因子；旧数据不清理。

## 尚未通过的门槛

已通过的厚度方向不替代周向、径向、联合细网格或SUPG检查。
现有周向6/8 V、径向6 V、SUPG敏感性失败记录均保留。
新试算完成后仍须检查实际载荷、原生残差和更新、独立物理通量、功率差额、局部载流子和环带变化。
局部载流子2%、SUPG全局/环带1%及其他原门限不变；不选较小功率差额的系数来代替通过。

九个生产工况仍为待验收：r3的6/6.5/7/7.5/8 V，以及r2/r4各6/8 V。
几何通过、输入准备完成、已有粗网格解均不能直接升级为生产或空间精度证书。
每个孔径需要适用的网格/稳定化证据。正式高电压生产图尚未生成。

## 图表使用范围

S1/S2最终图为release_v2的PDF/SVG/PNG；已检查面板对齐、最小文字尺寸、碰撞及实际面板内容。
S2d的原PDF包含对数刻度和10^-8、10^-10两条门限；另有独立局部渲染核查，图件未修改。
初版排版失败图保留为superseded，不能取代release_v2。
日志缺少全场FEM温度节点时，只展示真实峰温点，未用解析曲线伪造FEM曲线。

新生产图按已确认模板和实际数据共用色标。总等效应力与本构热膨胀项分开标注；
JE为简化总电流做功热分支，不能宣称完整热电闭合或已经解决高Pe问题。
