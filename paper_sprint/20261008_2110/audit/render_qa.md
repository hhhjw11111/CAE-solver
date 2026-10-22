# 最终PDF视觉核查

状态：FINAL_VISUAL_QA_COMPLETE，仅指本地阅读稿的呈现检查，不代表科学结果、空间验收或期刊模板已通过。

| 文件 | 页数 | 最终SHA-256 |
|---|---:|---|
| reports/study_report_final.pdf | 12 | 6e42657e7cf9b10e8046c88c88780d9c8139f15170902ede780fcbbb81094015 |
| manuscript/manuscript_working_draft.pdf | 15 | cae0444cf43bbdd95d19cb34a55f436b28ad88a3110643739daf2f815dc47714 |

两稿以已有Pandoc和XeLaTeX生成，使用Cambria、Microsoft YaHei与Cambria Math。禁用包安装和shell escape。Poppler以105 dpi渲染所有页面，退出码0；页面数与PDF页数一致。对应命令、页图和联系图保留于audit/render_qa/各文件名/上述哈希前12位目录。该大体积视觉检查目录不收入精简ZIP。

研究报告最终版检查了3张联系图覆盖全部12页，并全页复核第7、8、12页。标题、正文、公式、表格和8幅图均可见，页码连续，没有观察到裁切、图文重叠、乱码或缺图。第6—7页表格跨页时重复表头，属于正常续表；第12页为自然文末留白。早一版dd099628b083在第7页底部出现Pe图裁切，已以独立图注及不可拆分图块修复；中间版ac09541fa63c剩余一行段落跨页问题，最终以club/widow控制消除。旧位图保留，可核对修复轨迹。

论文最终版由独立审查逐一检查全部15张全页图，图表、中文、公式和文献均完整。修复图注重复、图注与图分离、验证表题孤立及近空白页后的最终记录见manuscript_pdf_qa.md。图5是稿内明确标注的待补证据项，没有生成虚假占位曲线；8幅已交付图件与正文引用保持对应。

渲染日志保留MiKTeX未由管理员检查更新的提示；命令退出码为0且实际PDF已目视检查。没有为消除该提示安装或更新软件。最初沙箱缓存故障和参数故障日志保留，不计作最终成功运行。DOCX因未找到可用Word/LibreOffice渲染环境未交付；可编辑正文为Markdown。

WORKING_DRAFT、H0来源PENDING、R23空间FAIL及其他科研限制均未因排版完成而升级。
