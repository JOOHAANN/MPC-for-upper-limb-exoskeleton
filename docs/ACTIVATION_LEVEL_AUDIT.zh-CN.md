# 240528Pilot、240617Pilot 的 EMG 与 activation 幅值核查

这两组旧 Pilot 的可核对弯腰记录中，**没有出现新 Exp 那种 activation 明显高于 EMG 的转换**。六个对应背部通道都满足形状参数为 **+0.1** 的指数映射，曲线接近恒等关系，activation 比其模型输入 EMG 略低。

| 数据组 | 找到的弯腰导出文件 | 能够完整核对 | 符合 +0.1 映射 | activation / EMG 的 P95 比值中位数 |
|---|---:|---:|---:|---:|
| 240528Pilot | 30 | 30 | 30 | 0.975 |
| 240617Pilot | 34 | 30 | 30 | 0.969 |

这里按文件计数，包含补录；不是 30 或 34 个独立工况。P95 比值先在同一文件里，对六通道均值分别取第 95 百分位数并相除，再在文件之间取中位数，不能解释为每个采样点都是固定的百分比。

![旧 Pilot 与新 Exp 的激活转换差异](activation_level_audit/pilot_EMG_activation_comparison.png)

## 核对方法与原因

读取各文件的原生时间戳，将 `emg.sto` 与 `Activations.sto` 精确配对，通道按 IL_L1、LTpT_T7、LTpL_L1 左右六个代表肌肉对应。逐通道检查

\[
a=\frac{\exp(Au)-1}{\exp(A)-1},\qquad A=0.1.
\]

60 份可核对文件全部吻合，最大 RMSE 约 **7.24×10⁻¹⁶**。单通道 activation 与其 EMG 输入的最大绝对差约 **0.01250**，即 0–1 数值范围的 1.25 个百分点。

同样输入经过两种转换的数值如下，说明两张图为什么有明显差别：

| 单通道输入 u | 旧两组 Pilot：A=+0.1 | 240730Exp：A=−2.7789 |
|---:|---:|---:|
| 0.20 | 0.192 | 0.455 |
| 0.40 | 0.388 | 0.715 |
| 0.60 | 0.588 | 0.865 |

`240617Pilot/subjectUncalib.xml` 中确实保存 +0.1。240528Pilot 本目录中未找到对应的已标定 XML，因此其 +0.1 关系是直接从现有 EMG/activation 数值验证的，没有宣称已经找到了当时实际加载的完整配置文件。

对完整的 `emgFilt.sto` 也作了同时间段、同六个传感器通道对照。传感器时间戳插值到模型输入时刻后，两组按文件统计的 EMG 对 EMG 差值 RMSE 中位数约为 **0.00354、0.00301**；不是新 Exp 那种系统性大幅抬高。插值只用于对照，不能证明硬件完全同步，也不修复来源断档。图中两阶段统一使用 5 ms + 300 ms 平滑，不按各自峰值重新归一化。

## 四份过短的补录

240617Pilot 的以下文件只有 1–2 个 EMG/activation 导出样本，不能判断动作曲线，也没有混入上面的 30 份统计：

- `CalibrationData/dynamics_stoop_5kg_53_1`：1 个样本。
- `CalibrationData/dynamics_stoop_5kg_53_2`：1 个样本。
- `CalibrationData/dynamics_stoop_5kg_53_3`：2 个样本。
- `CalibrationData/dynamics_stoop_5kg_53_4`：1 个样本。

这次结论只针对 EMG 与 activation 幅值处理关系，不是重新宣告全部试次都满足 BSE/MAP 训练质量要求。原文件、旧 BSE 参数和 MPC 参数均未修改。

## 后续模型输入

新 Exp 与旧 Pilot 使用不同的激活转换约定，不能直接凭图中 activation 更高就判断新受试者发力更大。统一输入定义后，需要针对各 subject 校准转换与 BSE 参数；已有 activation 不能再次套转换。参数重新训练属于后续标定，不是本轮为让曲线靠近而改数值。

复查运行：

```powershell
python experiments/multi_subject/compare_pilot_signal_levels.py
```

逐文件结果见 [CSV](activation_level_audit/pilot_level_audit.csv)，原文件 SHA-256、短文件样本数及详细计算见 [JSON](activation_level_audit/pilot_level_audit.json)。新 Exp 的黄色与蓝色说明见 [Exp 数据复核](EXP_ONE_CYCLE_REVIEW.zh-CN.md#黄色为什么低于蓝色)。
