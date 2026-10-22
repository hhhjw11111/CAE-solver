# S1/S2 补充验证图与数据

**交付版本为 `release_v2/`。** 本目录根部的首轮图和碰撞审查是保留的排版失败记录，不用于论文。没有修改或覆盖任何历史求解结果，没有新增原生求解。

- [S1：常浓度焦耳热验证](release_v2/S1_A1JJE_summary.pdf)
- [S2：非零扩散电流做功验证](release_v2/S2_JEGRAD_summary.pdf)
- [图注、能量/误差摘要及解释边界](release_v2/README.md)
- [S1 全部 13 例数据](release_v2/S1_all_cases.csv)
- [S2 全部 10 例数据](release_v2/S2_all_cases.csv)
- [所有原始指标及日志行号](release_v2/source_metrics_with_lines.csv)

图源为 `stage3_verification/attempt01` 已通过的 A1JJE、JEGRAD 日志；生成前后均调用冻结复用审计，重新核验真实文件及构建身份。原 A1/A2、A2P FAIL、原整体 FAIL、复用 acceptance 和 Fortran 全部保留。

绘图脚本为 `solver_v2/scripts/plot_je_paper_supplements_20261018.py`。重新生成时必须指定一个新的输出目录；脚本拒绝覆盖已有图件和科学证据。

```powershell
python solver_v2/scripts/plot_je_paper_supplements_20261018.py --output solver_v2/runs/JE_RECONSTRUCTION_20261018/paper_scope/supplemental_figures/new_reproduction
```

需要 Python、NumPy、Matplotlib、PyMuPDF，以及本机 `nature-figure` 技能的审查脚本。PDF/SVG 为可编辑矢量图，PNG 为 600 dpi 预览；每幅图 183 × 154 mm。未针对某个尚未指定的期刊声称符合其全部投稿规范。
