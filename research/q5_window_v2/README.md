# Q5：有限校准窗口与独立恢复审计 V2

**科学状态：受限类 M/M2 的有限误差条件与区间计算已完成；完整类 C 尚未满足 TAP 投稿门槛。** 本提交没有新的自适应采集算法优越性结论。A5、Q5 V1 及失败选择器保持原样。

## 阅读顺序与版本冲突

以本 README、`MANUSCRIPT_ZH.md`、`PRIOR_WORK_ZH.md`、`RESULTS_SUMMARY.json` 和绑定的源码/原始数组为本次执行依据。分支已有 `REPORT_ZH.md` 另版文字报告，其中的2/12→7/12、已知真实增益+位移 oracle、模态凹性等说法，不对应本次已核验的源码及数组；本次不采用、不合并，也不把它们宣称为独立复现。该报告及较早的多个 registration 文件保留，防止覆盖失败历史。

本次可核验结果是：无参考3/12、带噪参考GLS5/12、固定/随机附加EM各3/12；次级丢弃参考相位6/12、只给定真实位移仍5/12。相同seed并不足以保证不同实现的随机数组相同。这些结果始终只作开发诊断，不作未触碰最终集或确认性统计声明。正式后续试验必须另注册独立版本及数据流。

已有协议主文件为 PROTOCOL.md、EXECUTION_SPEC.md、M2_SPEC.md；本次代码与参数明确绑定 SHA-256。`run_v2.py`：`d9a9930cdfd1016f4ef1255796ee49d3c4c64387021f9c84e4cbcd91d4b368bb`；`multipole.py`：`88d7c62451d79ef03b02b1b6dd1edd208ecf6c888d085ff7b94318fb789750c7`。

## 交付与证据

`modal_interval.py`：5500盒名义材料覆盖、精确尺寸参数与外向有理端点；`finite_window.py`：5100盒有限差分覆盖、550盒尺寸/小损耗复条带、有限参考精度条件及集合中点估计；`radial_window.py`：理想电偶极径向窗口的精确有理二分；`multipole.py`：独立的经典向量 Maxwell 多球求解器；`run_v2.py`：12场景、4主基线、2次级消融；`checks.py` 与 `audit_tests.py`：物理一致性和证据完整性检查。

完整证据包包含原始复观测、所有初值结果、全部区间trace、未完成前缀及恢复记录。Git 中保存可读源码、报告与紧凑结果；**约15MB的完整trace与详细数组随本次交付ZIP提供，不假称已经全部写入Git**。也可在新目录运行以下命令生成全套。`MANIFEST.json` 是交付文件清单；JSON里的显示小数不是认证端点，真正边界使用完整证书中的有理数字符串。

## 运行

从仓库根目录运行，建议 Python 3.13；本次环境 NumPy2.3.5、SciPy1.17.0、mpmath1.3.0。不需要Treams。无网络环境请预装依赖；本次pip安装Treams/python-flint因DNS失败，按协议使用SciPy备选实现。不要用 `python -O` 禁用断言。

```bash
python -m pip install numpy==2.3.5 scipy==1.17.0 mpmath==1.3.0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python research/q5_window_v2/checks.py
python research/q5_window_v2/modal_interval.py
python research/q5_window_v2/finite_window.py
python research/q5_window_v2/radial_window.py
python research/q5_window_v2/run_v2.py
python research/q5_window_v2/audit_tests.py
```

modal_interval、finite_window、run_v2 默认拒绝覆盖已有实验输出。完整ZIP已有输出，可直接运行audit_tests；重做时对这三个脚本分别使用 `--out /新的对应目录`。audit_tests 默认审计交付位置，不会自动指向替代目录。finite_window 被外部时间限制终止时可用同一个目录加 `--resume` 恢复；恢复前应保存旧前缀哈希和源版本，不从未覆盖前缀宣布成功。

检查结果：物理一致性5组通过；证据/边界审计5组通过；M2有界误差100个数学样例通过。100样例不是100个独立电磁恢复实验。边界残差、阶数差和网格差都是数值诊断，不是 Maxwell 连续误差上界。

## 计算与失败记录

V2独立恢复运行约21.89秒、单BLAS线程、峰值约640.10MiB，DDA网格312/526/1016单元。区间名义完整覆盖约22.79秒。M2初次执行因200秒工具时限中止，后续使用保存前缀恢复；其JSON中seconds约10.69仅是恢复段，不能报告为全部证明用时。较早名义区间执行的2761行未完成前缀也保留为aborted，不是覆盖证书。具体时间依赖硬件，不构成方法速度优势。

## 不能跨越的结论边界

16次独立重复方案给出了非空的条件性模态误差预算；0.2%模态泄漏、0.1%相对通道增益、0.03%尺寸、参考系统偏差0.0005仍是待实测验收要求。逐样本逆函数用数值brentq，未把估计端点外向舍入。类C的可靠模型误差上界、覆盖有限分离、真实模态读出/漂移验证与最近邻全文原创性审计仍未闭合。本轮没有把失败证书变成普遍不可能性，也没有把PR完成当成论文完成。
