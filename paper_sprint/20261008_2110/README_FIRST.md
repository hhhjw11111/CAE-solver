# 2026年10月8日研究报告与论文工作初稿

**状态：WORKING_DRAFT。** 已完成正文、现有数据核查、10个3 V状态后处理、8组图及必要的小型验证；尚未完成全部空间/模型验收，因此不是可投稿定稿。

**版本范围：** 本轮新验证冻结于0.15.9 / R25_CLEAN。终检发现外部并行工作将源码更新到0.16.0-trust-region；本任务保留这些更改，不宣称已验证新版本。见[冻结说明与变化清单](audit/FROZEN_VERSION_NOTICE.md)。

全部58份启动源码已从既有冻结快照按SHA-256精确恢复到`audit/frozen_r25_source/src/`，用于核读本次报告中的源码行号，不替换工作区源码。

本地根目录：`C:\Users\Administrator\Desktop\codex-workspaces\Console3\solver_v2\paper_sprint\20261008_2110`

## 先阅读两份正文

- [详细研究报告 PDF](reports/study_report_final.pdf) · [可编辑 Markdown](reports/study_report_final.md)
- [论文工作初稿 PDF](manuscript/manuscript_working_draft.pdf) · [可编辑 Markdown](manuscript/manuscript_working_draft.md)
- [首版报告](reports/study_report_v0.1.md)与[首版初稿](manuscript/manuscript_v0.1.md)保留，证明先成文再持续核验的工作过程。

## 证据与结果

- [证据登记 JSON](audit/evidence_register.json) · [CSV索引](audit/evidence_registry.csv) · [逐结论可用范围](audit/claim_matrix.md)
- [10个算例的七层状态](data/case_registry.csv) · [完整指标 CSV](data/metrics_summary.csv) · [功率分项](data/power_audit.json) · [空间失败记录](data/spatial_comparison.json)
- [后处理详细报告](analysis/postprocess_report.md) · [图件与取样契约](analysis/figure_contract.md)
- [当前源码模型契约](audit/source_validation/source_contract.md) · [基线与当前差异](audit/source_validation/MODEL_CONTRACT_DIFF.md) · [源码行号和哈希](audit/source_validation/source_index.md)
- [本轮小型验证](validation/runtime_retry/verification_summary.json) · [纯数学参考检查](validation/formula_reference/formula_reference_summary.json) · [文献核验](audit/literature/literature_verification.md)

`figures/`保存8组PNG/SVG/PDF，共24个图文件；图1是几何示意，图2是历史A1/A2实际误差，其余6组为原始有限元场的后处理。`data/`包含40份CSV图源和剖面，另有JSON指标；`scripts/`保留后处理、解析参考、来源核对和排版脚本。完整大体积节点场未复制，在逐case审计中提供原位置、大小和SHA-256。

## 使用边界和未完成事项

历史六组F1结果均为L0探索性数据，空间验收未完成；H0三例现存档案身份匹配，但缺运行当时逐case进程/环境记录，来源状态PENDING。3 V R23对数加权半域局部载流子约3.01%差异超过2%门槛，保留FAIL。NT/N0模型差值尚未有匹配多网格证据；当前漂移电导热源与材料输入不能证明真实器件的精细温度反馈。

本轮新三维求解为0。21:15资源核查时可用内存约2.36 GiB，另有两项用户求解；本任务不干预这些进程。五项当前CLI验证与一个R22归档生产核条带重放成功；后者不满足附件完整一维规格，纯数学检查也不能替代生产验证。

- [OPEN_ITEMS及最小下一步](OPEN_ITEMS.md)
- [本轮改动说明](CHANGELOG.md) · [执行决策记录](audit/decision_log.md)
- [源码/未提交修改/基线保护检查](audit/protection_check.json)
- [PDF视觉检查](audit/render_qa.md) · [论文逐页检查](audit/manuscript_pdf_qa.md)
- [交付清单](audit/run_manifest.json) · [文件SHA-256清单](audit/delivery_files_sha256.csv)

Markdown为可编辑正文，PDF为阅读稿。当前未找到可用Word/LibreOffice渲染器，未交付未经视觉核验的DOCX。新分支为`codex/paper-sprint-20261008-2110`；无提交、push、合并、标签移动或生产源码修改。

## 复现与本地打包

后处理具体入口及输入见`analysis/postprocess_report.md`和`audit/postprocess/`，不要将纯数学参考脚本当成求解器。PDF生成脚本是`scripts/render_documents.ps1`，使用已有Pandoc/XeLaTeX并禁用安装和shell escape。成功验证命令、退出码及初始DLL环境失败都在`validation/`保留。

完整本次成果包为[CAE_solver_paper_sprint_20261008.zip](CAE_solver_paper_sprint_20261008.zip)。它包括文稿、图件、源数据、脚本与审计；不包括原大型节点场、逐页QA位图和启动时用户差异补丁。差异补丁仍保留在本地audit目录。压缩包未向任何外部服务上传。
