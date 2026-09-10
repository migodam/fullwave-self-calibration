# A4 收尾 → GPT Pro 下一轮交接

主 Prompt 为 `Theory/Questions/Q5.md`，中文汇报为
`communication/A4_COMPLETE_RESEARCH_REPORT_ZH.md`。先阅读汇报，再把 Q5
全文作为任务发送；不需要把旧 Q4 再作为执行指令发送。

最新配套上传包：`communication/A4_GPT_PRO_HANDOFF_20260910_v2.zip`。
其中 Q5 已按材料归因因果链、两独立材料区域、四机制消融和 L/F/R 三层
Gate 修订；请使用标有 v2 的 Prompt。原 `A4_GPT_PRO_HANDOFF_20260910.zip`
保留为历史版本，不含这次修订。包内原 A4 收尾审计记录的是当时 v1
文档哈希，不是 v2 的文档审计，不能拿它验证新版 Q5。
里面有 Q5、本报告、A4_2、原英文稿、实验与文献审计、A4 原始代码包、
最新英文增补、机制验证代码／原始结果／审计、依赖的 A3 Maxwell 适配器。
不包含凭据、环境目录或第三方论文 PDF。

包内 `Theory/Questions/A4_fullwave_research_v1.zip` 是保留不动的原始包，
解压后得到 `research/a4_reliability_v1/`。新机制脚本仍沿用本地工作区
`public_release/research/a4_reliability_v1/` 导入布局；若在网页端重跑，
需要按脚本导入路径布置原包或明确修改适配路径，并记录修改。
这不是已验证的全新环境一键安装包。代码执行还需 NumPy、SciPy、Treams；
原 A4 的测试环境说明随原包保留。

若 GPT Pro 不支持读取 ZIP，直接上传 Q5、中文汇报、A4_2 和英文增补，
再按它实际需要上传原 A4 代码／理论。Q5 已包含核心事实与反例，
不应因缺失全文而虚构文献核验或实验执行。

最新文件本轮未推送 GitHub。此前公开仓库不自动包含这些本地成果。
阶段交接完整不代表 TAP 投稿就绪；原生模型科研 pipeline 本轮未重启。
