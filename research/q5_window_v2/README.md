# Q5 V2：有限误差信息窗口、参考精度与独立恢复审计

## 阅读顺序与科学状态

1. [FINAL_AUDIT_ZH.md](FINAL_AUDIT_ZH.md)：最终重放、版本对账、已完成与未闭合门槛。
2. [MANUSCRIPT_ZH.md](MANUSCRIPT_ZH.md)：中文论文段落、模型、定理、证明、预算、Information Window Splitting及天线设计含义。
3. [PRIOR_WORK_ZH.md](PRIOR_WORK_ZH.md)：最近邻全文公式比较与未取得全文的项目。
4. [RESULTS_SUMMARY.json](RESULTS_SUMMARY.json)：绑定源码与观测的结果摘要。

本轮M/M2的受限有限误差条件和区间计算已完成，类C没有达到TAP强正结论门槛。没有新的自适应采集优势声明，不恢复已撤回选择器，不修改A5/V1。

`REPORT_ZH.md`是保留的另一版历史报告；其中2/12→7/12、oracle增益加位移和模态凹性等结论没有与当前源码/数组对齐，不并入本轮证据。当前主计数为无参考3/12、带噪参考GLS5/12、固定/随机EM各3/12；次要消融为丢弃参考相位6/12、已知真位移5/12。所有计数都是开发诊断，不是最终测试或总体成功率。

## 协议、实现与证据

实施范围见PROTOCOL.md、EXECUTION_SPEC.md、M2_SPEC.md。其他注册文件保留历史，不据其拼接不同数组。

- `modal_interval.py`：5500盒名义材料覆盖，精确有理尺寸，外向有理端点和级数尾界。
- `finite_window.py`：5100盒有限增量、550盒尺寸/小损耗复条带、参考条件和集合中点估计。
- `radial_window.py`：理想电偶极材料任务距离窗口的精确有理二分。
- `multipole.py`：经典向量Maxwell多球展开，不依赖Treams，不作为算法创新。
- `run_v2.py`：12例独立DDA数据，4主基线和2次要消融；原始复观测、所有初值和失败均保存。
- `checks.py`、`audit_tests.py`：物理一致性与证据检查。

实际环境：Python3.13.5、NumPy2.3.5、SciPy1.17.0、mpmath1.3.0，单BLAS线程。最终重放5项物理检查和5项证据检查全部通过；100个M2例子是数学误差盒测试，不是电磁恢复场景。

## 重放包

对话交付的`Q5_V2_REPLAY_EVIDENCE.zip`为重新执行后生成的证据包，29个文件，压缩后2,667,910字节。源码、结果和日志约10.79MB，不包含付费书籍或字体。它不是此前运行目录的磁盘备份；原始观测和两份完整覆盖trace与已有记录逐字节匹配。完整清单见包内MANIFEST.json。

ZIP SHA256：`0fd1672cb6220cec517b35f8f9b3b09133d36585b134fea88fca522b3b0c759d`。

`run_v2.py` SHA256：`d9a9930cdfd1016f4ef1255796ee49d3c4c64387021f9c84e4cbcd91d4b368bb`。

`multipole.py` SHA256：`88d7c62451d79ef03b02b1b6dd1edd208ecf6c888d085ff7b94318fb789750c7`。

Git中提交可读源码、论文段落和紧凑摘要；完整trace与原始数组在上述交付包内，也可用脚本重新生成。不要将JSON显示用浮点数当作证明端点，证明端点采用trace中的精确有理字符串。

## 运行

从仓库根目录、使用一个尚无输出的新结果目录运行：

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p research/q5_window_v2/results
python research/q5_window_v2/checks.py
python research/q5_window_v2/modal_interval.py
python research/q5_window_v2/finite_window.py
python research/q5_window_v2/radial_window.py
python research/q5_window_v2/run_v2.py
python research/q5_window_v2/audit_tests.py
```

交付包已含结果，重新计算前应先把results目录改名保留；直接审计包内证据只需运行audit_tests.py。modal_interval.py、finite_window.py和run_v2.py拒绝覆盖既有主输出，也可用--out指定新目录；audit_tests.py默认审计上述标准目录。

有限覆盖可能需要数分钟。本次两次触及200秒工具时限后，以相同源码/精度使用`finite_window.py --resume`完成，最后续跑22.62秒不是总耗时。中断摘要见重放包REPLAY_RECORD.json，不能把续跑时间解释为完整计算成本。

## 不能据此宣布的结论

16次独立重复时，0.2%相对模态泄漏、0.1%相对通道误差、0.03%尺寸误差等可形成非空条件预算；这些是待验收规格，不是已经达到的硬件指标。逐样本Brent反演端点未外向舍入。类C连续前向误差、全域有限分离、真实读出/漂移验证及最终原创性审计仍未闭合。保持Draft PR，不标为TAP-ready。
