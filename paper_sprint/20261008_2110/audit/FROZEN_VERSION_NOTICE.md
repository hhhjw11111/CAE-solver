# 证据冻结范围与外部并行修改

本轮源码审计、行号引用及当前小型CLI验证实际针对0.15.9-log-inverse-checked / R25_CLEAN，可执行SHA-256为a47da18f47ceca0209ef9b6fcc17bb13c1175f4e3a3181ffc1eb4d7b7f2f6814。核查时58份Fortran及该二进制59项匹配。R22条带使用其另行归档核，F1、R23和R25原场分别按各自身份登记，不互相替代。

2026-10-08 21:53:25终检发现外部并行写入将工作区VERSION改为0.16.0-trust-region，8个Fortran文件与启动哈希不同，详见concurrent_source_changes.json。这些不在本任务写入范围；本任务没有回滚、清理或验证这些更改。src/physics相对基线仍无差异，不能由此推断新增求解算法已通过本轮验证。

因此source_validation/source_contract.md与source_index.md中的“当前”及行号，是审计时的冻结源码身份，可能不再适用于终检后的工作文件行号。基线标签和HEAD检查独立进行。启动/结束差异补丁保留；后续程序验证应重新建立自己的源码—构建—结果链。

为保持可核读性，已从既有R25/R23冻结快照找回全部58/58份启动源码，逐文件哈希与启动记录完全一致，缺项0。它们保存于`frozen_r25_source/src/`，恢复来源及复制后核验见`frozen_r25_source/frozen_source_manifest.json`。没有使用变化后的工作源码冒充旧版，没有重新编译或修改Git引用。
