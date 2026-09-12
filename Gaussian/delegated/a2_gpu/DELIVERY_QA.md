# GPU 交付机械 QA

审计时间：本地读取；未连接或操作远端。检查范围是 `code/a2/gpu/summarize_campaign.py`、三份交付入口与已同步的本地结果。根任务提供的运行状态为：4060 主运行曾在 case 18 因 SSH 断开中断，已以相同源码/配置恢复 cases 18--29（PID 25432）。这不是本地远端状态核验。

## 必须修复

1. `RESEARCH_REPORT_ZH.md` 第 45--55 行的表明示为 **3000 RHS 以内、每组十例中位数**；`figures/a2/extensions/manifold_budget.png` 由 `summarize_extensions.py` 第 40--44、65--71 行绘制的是每个预算点的 **均值**。两者数值可实测不同：分离组 spectral 在 3000 RHS 的材料误差中位数为 0.48%，均值为 1.91%。报告中的图注必须明确“图为均值，表为中位数”，或将图改为中位数后再把它作为表的可视化；目前“同预算上限内的实际质量”会让读者误以为统计量相同。

2. GPU 主结果仍不得生成或引用。当前本地 `main128/` 只有 18 个 `case_*_summary.json`、0 个 method reconstruction arrays、1 个 data archive；尚非 30×4。`summarize_campaign.py` 第 24--50 行会正确拒绝输出。任何状态更新都必须写“SSH 中断后以相同 hash/config 恢复运行，尚未完整同步/验收”，不能使用“连续完成”或“120 次大规模反演已完成”。

## 通过的口径核对

- 报告表的 3000 RHS 边界与原始 `manifold_v3` 的 `history.rhs <= 3000` 最后记录一致；表中的中位数和相邻组 held-clean 中位数均可复算。预算图对 1000/2000/3000/3400 也只取 `rhs <= budget` 最后点，且明确无插值/额外求解。
- `summarize_campaign.py` 对未来 GPU 主摘要检查 N128、120 steps、四种参数数 96/384/864/16384、config hash、source hash、RHS 计数、完整轨迹、arrays 与 data archive；完成门槛在代码语义上是严格的。其未来候选段落仅在检查通过后写入。
- `RESEARCH_REPORT_ZH.md` 当前 GPU 段落明确写“正在执行”“不计为 120 次已完成”，没有把 GPU 当作已完成证据。应只按上面的中断/恢复事实补充状态，不改变这一边界。
- 研究报告列出的四题“流形原生 SOM、可复用 n-port、可信结构增长、联合概率与 NN 推断”与 `GPT_PRO_QUESTIONS_ZH.md` 的问题一至四逐项一致；`START_HERE.md` 也链接同一四题入口。未发现旧三题包仍被作为当前问题包的冲突。

## 非修复性备注

本 QA 不作科学验收。概率与流形结果在本地已完整，但 GPU 恢复运行的进度、PID 和远端文件完整性没有在此处远程验证。
