# 结论与证据对应

本表的PASS只覆盖同一行所述命题，不代表整个算例全层通过。源码和本轮小CLI验证冻结于启动时0.15.9/R25_CLEAN，不覆盖之后外部并行更新的0.16.0版本；源码精确副本见frozen_r25_source/。

| ID | 具体命题 | 状态 | 证据及可写边界 |
|---|---|---|---|
| C01 | 启动时58份Fortran及R25_CLEAN二进制与清单对应 | PASS | source_validation/source_audit_runtime_retry.json；不宣称工作树干净或外部新版本已验证 |
| C02 | LQ/N0/NT、本构、材料和全域约束已按冻结源码核录 | PASS | source_validation/source_contract.md；有效模型定义，不是实测材料证明 |
| C03 | 六个历史F1 3 V L0完整运行来源可追溯 | PENDING（H0三例） | 现存54源/exe/input/mesh档案身份匹配，H0缺当时逐case进程/环境记录；H1有记录；保留空间PENDING和重复字节FAIL |
| C04 | 冻结R25小CLI V01/V02/G02/F1M/G03通过 | PASS | validation/runtime_retry原始日志；不是全部生产或空间测试 |
| C05 | R22归档生产HEX20条带可重放 | PASS（限定） | 61项身份匹配；自定义缩约装配，采样最大误差，非当前R25重编/完整CSR |
| C06 | 附件完整一维生产规格已经完成 | PENDING | 缺同电子负场、零场、L2独立积分稳定性和完整自由行记录；不能以公式检查替代 |
| C07 | 五种A的解析参考实现符合端值/恒流/零场极限 | PASS（纯数学） | validation/formula_reference；production_solver_validation=false |
| C08 | 3 V R23对数加权半域局部载流子场通过2%空间门槛 | FAIL | data/spatial_failure_source.csv；电子/空穴约3.0065%/3.0098% |
| C09 | NT−N0精确温度修正已被网格分辨 | PENDING | 缺同电压同格式多网格模型差值；摘要和结论不写精确百分比 |
| C10 | 所有局部热源代表真实器件热量 | BLOCKED | 漂移电导约化、扩散抵消、材料和边界未标定；禁止器件真实性强结论 |
| C11 | H0/H1环带比值可直接做定量几何结论 | FAIL（既有积分定义） | 旧H0区域体积偏差约28.48%；已删除其增强比值 |
| C12 | 孔外同坐标厚度、应力与功率可由原场重构 | PASS（后处理一致性） | data/profiles、metrics_summary及逐case重构检查；数值场空间认证仍未完成 |
| C13 | 常系数绝缘绝热含孔域温度参考关系成立 | PASS（在假设内） | 报告第4节自洽推导；参考峰值3.695019862 K非新FEM结果 |
| C14 | 已完成实体实验及完整器件标定 | NOT_APPLICABLE | 本轮实体实验数0，稿件没有实验验证结论 |

详尽逐case状态使用provenance_status、geometry_status、solve_status、flux_status、spatial_status、model_scope_status与claim_status七字段，见case_registry.csv。
