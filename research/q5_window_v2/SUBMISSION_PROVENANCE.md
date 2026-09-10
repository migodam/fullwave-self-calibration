# 本轮提交来源与并行分支保护

本轮科学基线为 418ef53a7d7d2f79d6228d023baf8b101e3b3e2e。最初协议在 research/q5-certified-windows-v2 分支登记，PROTOCOL.md 的内容 blob 为 fa9af6f7524fc44f2e293e20d132a865c5507301。

提交审计期间，原分支出现并行提交；因此本交付单独使用 research/q5-window-v2-audited-delivery 分支，直接以基线提交为父提交，只增加 research/q5_window_v2，不覆盖、强推或删除原分支，也不修改 A5 / Q5 V1。

PROTOCOL.md 和 IMPLEMENTATION_FREEZE.md 描述已执行的主数值流；ATTRIBUTION_DIAGNOSTIC_PROTOCOL.md 是检查主开发结果后另行冻结的事后归因诊断，不冒充最终测试。其他历史登记文本原样保留，仅作来源记录。旧 REPORT_ZH 已保存在 review_inputs/REPORT_ZH_pre_audit.md；以本目录当前 REPORT_ZH、数学证明、原始数据和审计快照为准。

仓库保存可运行源代码、数学段落、证明证书、48 项逐场景/方法结果及原始数据哈希。完整原始复数组、每个优化初值、中断检查点和运行日志通过本次对话的 Q5_V2_RESEARCH_AND_EVIDENCE.zip 交付。GitHub 中的轻量摘要不替代原始证据，未声称已上传该归档至 GitHub。

提交设为 draft：受限数学证明可审阅；C 类有限材料分离、实际连续模型误差和读出硬件验收未闭合，不宣称 TAP 终稿完成。没有启动 GitHub CI 或自动合并。
