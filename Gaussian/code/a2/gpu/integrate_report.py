"""Integrate validated, completed GPU evidence into the existing A2 documents."""
from pathlib import Path
import json,re
ROOT=Path(__file__).resolve().parents[3]
RUN=ROOT/'runs/a2/gpu/representation_campaign/main128'
def section(path,tag,content):
    s=path.read_text();pattern=f'<!-- {tag}_START -->.*?<!-- {tag}_END -->'
    assert re.search(pattern,s,re.S);path.write_text(re.sub(pattern,lambda _:f'<!-- {tag}_START -->\n{content}\n<!-- {tag}_END -->',s,flags=re.S))
def main():
    s=json.loads((RUN/'SUMMARY.json').read_text());receipt=json.loads((ROOT/'runs/a2/gpu/execution_receipt.json').read_text());assert receipt['complete_case_methods']==120 and receipt['resumed_ssh_exit_code']==0
    names={'off_grid_gaussian_mixture':'离网格多Gaussian','sharp_multi_inclusion':'多夹杂尖边','curves_with_holes':'细曲线与局部挖空'}
    pretty={'gaussian_K16':'Gaussian K16','gaussian_K64':'Gaussian K64','gaussian_K144':'Gaussian K144','voxel':'体素'}
    lines=['## 4060大规模分支：参数规模扩大后，哪些优势保留下来','',
    '**已完整执行30个对象×4种表示＝120次反演。** 反演网格为128×128（16384个复状态单元/频率/发射），数据来自192×192单元积分模型。三频1.25、2.0、2.75 GHz，12个圆周线源、64个半圆接收，发射半径34cm、接收半径30cm，成像视场40cm×40cm。偶数接收拟合、奇数接收留出。', '',
    '三类各10个对象：16–32个离网格平滑Gaussian、多圆盘/矩形尖边夹杂、相交细曲线与局部挖空。K16/K64/K144分别有96/384/864个实材料参数；体素有16384个。每例120次Adam更新，使用同一个全波后端和同形式的离散平滑项。全波状态仍需求解，Gaussian基函数解析并没有使全波场变成解析。', '',
    '**下表为各组10例均值。** 留出误差相对无噪声生成场，和上文概率实验的带噪留出指标不同。', '',
    '|目标|表示|材料相对L2|留出场误差|','|---|---|---:|---:|']
    for fam,v in s['families'].items():
        for m,x in v.items():lines.append(f"|{names[fam]}|{pretty[m]}|{100*x['material_l2_mean']:.2f}%|{100*x['held_clean_mean']:.2f}%|")
    lines+=['',f"![三类目标的实际重建，同一色标与视场]({ROOT}/figures/a2/gpu/representation_main128_representative.png)",'',
    '**配对结果：** '+'；'.join(f"{names[f]}中K144的材料误差低于体素为{x['k144_material_win_count']}/10例" for f,x in s['paired_k144_vs_voxel'].items())+'。这些是本次对象和固定设置下的计数，不是跨分布成功率。','',
    '视觉上，更多Gaussian能保留更多小结构，并让细曲线更连贯；K16容易将它们合并成宽带。尖边仍会变圆，较低留出场误差并不总对应较低材料误差。K144还会出现很淡的长条/斜条背景伪影，不能解释成已经确认的材料细节。具体观察和改进假设见[视觉审查](../delegated/a2_gpu/VISUAL_REVIEW_ZH.md)。', '',
    f"![每组K144相对体素的最好和最差案例]({ROOT}/figures/a2/gpu/representation_main128_k144_vs_voxel_best_worst.png)",'',
    '上图每组按K144材料L2减体素材料L2的差值取最小/最大，展示相同对象的真值和两种重建；差图蓝色表示该像素K144误差较小。选择只用于展示，不进入优化或超参数选择。', '',
    '### 参数少了，成本怎样变化','',
    '|表示|实材料参数|优化与逐步存档平均秒/例（29例）|相对体素的配对耗时比中位数|','|---|---:|---:|---:|']
    for m,x in s['paired_cost_excluding_interrupted_case18'].items():lines.append(f"|{pretty[m]}|{s['methods'][m]['parameter_count']}|{x['seconds_including_checkpoint_io_mean']:.2f}|{x['paired_seconds_ratio_to_voxel_median']:.2f}×|")
    c=s['cost_accounting']
    lines+=['',
    '第19个对象（case18）的两种方法在系统重启前完成、另两种在重启后完成，因此只在配对计时中排除，质量统计仍保留全部30例。方法顺序固定、未做独立重复计时；该表是执行策略下的时间，包含每步checkpoint I/O，不能当作纯CUDA速度。', '',
    f"保存的120次完整优化时间合计{c['completed_optimization_seconds_including_checkpoint_io_sum']:.1f}秒；记录到的数据生成/首次构建合计{c['recorded_data_generation_seconds_sum']:.1f}秒，反演算子构建合计{c['recorded_inverse_setup_seconds_sum']:.1f}秒。这不是无中断端到端总耗时：原case18被中断的部分工作、部分重复初始化、最终评价/绘图/同步与停机时间没有被完整分项计量。跨两次进程记录到的PyTorch分配显存峰值为{c['peak_cuda_allocated_bytes_across_invocations']/2**20:.1f} MiB；这是共同过程峰值，不能用于不同表示间的峰值比较。", '',
    '每个完整反演包含4356个正演RHS和4320个伴随RHS（含最后36个正演RHS评价），各表示相同。隐式梯度不需要逐材料参数构造完整Jacobian；Gaussian渲染仍有O(NK)开销，体素映射为O(N)。因此864对16384、约19倍的参数差并没有转化为速度优势。完整复杂度见[复杂度说明](../Theory/a2/extensions/COMPLEXITY_ZH.md)。', '',
    f"![所有完成对象的成本与质量，不隐藏离群值]({ROOT}/figures/a2/gpu/representation_main128_cost_accuracy.png)",'',
    '另有256×256＝65536状态单元的单频正演/伴随规模检查。当前GPU后端已重新对照SciPy场参考（相对差约5.30×10⁻⁷）、隐式梯度有限差分（约1.51×10⁻⁴）及批量多频一致性；这些验证求解实现，不证明材料恢复或加速。', '',
    f"![单频已知材料的GPU规模检查，不是反演速度比较]({ROOT}/figures/a2/gpu/representation_gpu_scaling.png)",'',
    '### 这些结果意味着什么','',
    '**Gaussian的紧凑表示值得继续研究，当前证据尚未支持SOM、概率或NN的新主贡献。** 平滑多结构和细曲线给出正面线索；尖边与低幅伪影说明需要同时研究表示与数据支持。下一步应检验普通自适应Gaussian、TV/界面模型与完整GN，再决定是否需要NN。torch.nn.Module在这里仅封装待反演参数，不是跨对象训练的神经网络。', '',
    '本组没有额外测量噪声，仍含网格与积分失配；两种正演属于同一VIE软件体系，不能称完全独立软件验证。不同K从不同的初始材料图开始，体素从K64初始图开始；相同更新次数/步长/罚项不等于有效先验或调参难度相同。因此这是有意义的规模开发实验，尚非满足重复噪声、强基线调参、实测和3D Maxwell要求的论文验收。', '',
    '原运行在系统重启后中断；保留已完成结果，隔离损坏的case18 K144检查点后重算该方法并继续。重建数据与旧数据逐数组一致，最终源码/配置哈希、120步轨迹、求解计数和图像误差全部核对。详见[执行记录](../runs/a2/gpu/execution_receipt.json)。']
    section(ROOT/'deliverables/RESEARCH_REPORT_ZH.md','GPU_EXTENSION','\n'.join(lines))
    p=ROOT/'deliverables/RESEARCH_REPORT_ZH.md';t=p.read_text().replace('并开始真正扩大材料参数与状态规模','并完成了材料参数与状态规模的扩大比较');p.write_text(t)
    sentence='GPU规模分支完成30对象×4表示。K144以864个材料参数对照16384参数体素，在当前平滑多结构与曲线任务有质量收益，尖边仍有不同指标的冲突，且未取得配对速度优势。这支持继续研究紧凑表示，不自动支持原生SOM、概率或NN的新机制。'
    section(ROOT/'reviews/a2/FINAL_VERDICT_ZH.md','GPU_VERDICT',sentence+' 运行跨一次系统重启，质量统计保留全部对象，配对成本排除混合阶段case18；完整成本记录仍有声明缺口。')
    p=ROOT/'deliverables/GAUSSIAN_VS_CALIBRATION_ZH.md';t=p.read_text();t+='\n## 完成规模分支后的投入建议\n\n'+sentence+' 相较继续谱截断调参，更值得优先解决自适应材料表示与错误细节判别；Calibration仍按实际采集失配需要引入。本轮不能给出Gaussian与Calibration的公平算法排名。\n';p.write_text(t)
    p=ROOT/'deliverables/PAPER_WORKING_DRAFT_ZH.md';t=p.read_text().replace('所有数字、停止条件、图和样本边界以','8. GPU规模：N128反演/N192积分数据、30对象×4表示、120步隐式伴随优化，含系统重启与计时边界。\n\n所有数字、停止条件、图和样本边界以');t=t.replace('GPU规模实验、可选PyTorch实现和远端硬件记录加入同包，具体完成范围以主报告与持久化run记录为准。本稿不把运行中的任务写成已完成数值结果，也不把GPU加速当作Gaussian专有收益。',sentence+' GPU已完成，相关代码、数据数组、图和执行边界加入同包。');p.write_text(t)
    p=ROOT/'deliverables/GPT_PRO_QUESTIONS_ZH.md';t=p.read_text();pos=t.index('## 问题一：');summary='## 新增大规模证据及必须解释的现象\n\n'+sentence+' 数字、完整图、等RHS定义和跨重启计时范围以统一报告的4060章节及 `runs/a2/gpu/representation_campaign/main128/SUMMARY.json` 为准。这里用普通参数表示和Adam，不能称已经验证的新SOM；固定无噪声开发配置不能代替强基线与确认性测试。问题一需要击败现有隐式梯度实现，问题三需同时处理曲线连贯与低幅细长伪影。\n\n';p.write_text(t[:pos]+summary+t[pos:])
    print('Integrated completed GPU evidence into existing A2 documents.')
if __name__=='__main__':main()
