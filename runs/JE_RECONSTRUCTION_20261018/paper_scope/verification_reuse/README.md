# 已有验证证据复用审计

本目录只读重验已有日志、退出状态、逐字节回归及源码/构建哈希，没有重新运行求解器。

用户批准范围 JE_PAPER_MINIMAL_20261018：10 项专项全部通过，16 项旧 CLI 基准/候选逐字节一致。40 个物理/核心源文件与原协议及两份构建清单一致。

原 11 项专项中的 A2P 仍为 FAIL，原整体 result.json 的 all_passed=false 保留。A2P 因用户明确缩小论文范围而不再是必做项，并未获得数值豁免或改判。

执行身份：baseline ARC_JE_FINAL_20261012 与 candidate JE_VERIFY_20261018_FINAL 的可执行文件均重新核对 SHA-256。候选完整源码清单与当次 provenance 一致；当前 A2P 模块及验证 driver存在后来修改，逐一记录于 current_source_differences，不能声称当前全部源码等于已运行候选源码。所有其余源文件（包括 A1/A2、A1JJE、JEGRAD 及物理依赖）重新核对未变。

两份 build manifest 的核验范围是执行文件、冻结核心，以及候选执行时完整源码清单；其中历史文档条目不被误称为当前执行计划仍未变化。

acceptance.json 只接受本范围；网格、局部载流子、SUPG 和9个生产工况仍需后续验收。下游应调用 validate_acceptance(path, expected_exe=...) 重新检查实际文件，而非只读 PASS 布尔值。
