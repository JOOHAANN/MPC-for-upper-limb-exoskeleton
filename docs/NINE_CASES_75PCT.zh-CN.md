# 九工况 BSE 与 75% 峰值阈值 NMPC

在 MATLAB 中进入项目根目录，运行：

```matlab
demo_nine_case_75pct          % 全部 subject
demo_nine_case_75pct('725')   % 只查看 725 与总览
demo_nine_case_75pct('730')
```

这个入口直接显示已计算的结果，不重新训练。完整图集位于
`output/pdf/nine_case_bse_nmpc_75pct.pdf`，逐工况数值、模型哈希和问题清单位于
`docs/nine_case_75pct/`。以该目录的 `summary.json` 和 `nmpc_cases.csv` 为数值依据，
其中未达标的行也保留。这里报告的是同九工况的联合标定与回放，不是独立验证集表现。

本轮按用户要求停止继续训练和调参，固定目前选定的版本。528、725、730 的
54 条工况/来源结果都达到辅助后背部峰值位于阈值下方 0–5 N·m 的目标。
725 的两种来源在未辅助时九工况均有两个弯腰波谷；辅助后 15 kg / 45 各有一个
小附加波谷。730 仍有少量附加波谷，实验来源 10 kg / 65 的第二次弯腰偏浅。
这些不足保留在图集和问题清单中，不以继续试验的参数覆盖当前选择。
`results/nine_case_75pct/release_manifest.json` 记录冻结版本的文件哈希；同级
`candidates/` 内的其他训练属于
未选用的试验，尤其 `divisive_intent`、`soft_tonic` 不属于交付模型。

## 使用哪些工况

| Subject | 5 kg 的速度代码 | 10 kg 的速度代码 | 15 kg 的速度代码 |
| --- | --- | --- | --- |
| 240528Pilot | 35 / 45 / 55 | 35 / 45 / 55 | 35 / 45 / 55 |
| 240725Exp_A01 | 35 / 45 / 65 | 25 / 45 / 55 | 45 / 55 / 65 |
| 240730Exp_A01 | 35 / 45 / 65 | 25 / 35 / 45 / 65 | 45 / 55 |

730 是 3+4+2 个不同工况，合计九个。编号重复试次只作为该工况的数据候选，
不会重复计数。选择依据是时间戳、激活平滑度、完整动作和已有 CEINMS 力矩的覆盖情况，
选择发生在 BSE 拟合前。完全无 CEINMS 力矩覆盖的动作不进入这次训练；部分缺失保留
缺失掩码，画图时留空，不补造力矩。

## 参数固定到什么范围

每个 subject 有两个固定的 BSE 向量：实验 activation 一套，MAP activation 一套。
每套由 32 个物理参数和 10 个神经控制/激活解码参数组成，在该 subject 的九种工况中
完全固定。M、L、肌肉强度、肌腱/肌纤维、激活动力学、阻尼等都进入优化。
身体/负载几何和软限位的 8 个物理参数在两种 activation 来源之间进一步共享。
`fixed_parameter_audit.csv` 可以逐行核对，九工况中的参数哈希不会改变。

每个 subject 单独选择 NMPC 的权重、预测误差裕量及 activation 预测器。
同一 subject 在九工况中复用这套策略，没有按某一个工况重新选择控制权重。
两种输入来源的预测误差统计不同，所以误差补偿分别由其标定残差计算。
这是固定策略下的统计量，不是逐工况手工加大辅助。

## 角度、力矩与等待平衡

角度仍使用原项目 `Scale` 生成的参考，选取动作后不再次缩放。
实验来源使用 A=+0.1 后的实验激活，MAP 来源使用已训练的 MAP activation 预测器。
730 的输入在原处理之后另外统一使用 310 ms Gaussian 窗口；两种来源都采用相同
附加窗口，不能用不同工况的角度反过来修改 activation。

两条肌肉链都包含 Hill 收缩元件、被动元件、肌纤维阻尼和弹性肌腱。它们与关节方程
在每个时间步耦合求解：

\[
\tau_{net}=\tau_{back}-\tau_{abdomen},\qquad
J\ddot\theta=\tau_{net}+\tau_g+\tau_{boundary}+u-b\dot\theta.
\]

这是同一个肌肉-肌腱-关节前向模型的状态更新；没有另开角度跟踪方程来替换输出，
也没有把算出的 theta 截断到参考范围。等待阶段的姿态 PI/PD 请求经过肌肉激活和
同一个 Hill 模型产生真实模型内肌肉力矩，不直接叠加一个关节 PD 扭矩。
NMPC 的 u 是唯一的外骨骼辅助力矩；软限位单独记录，绝不算进肌肉力矩。

当前 Exp 的拿起/放下时间由实验 EMG 的两次主峰估计，负载接触使用 80 ms 的平滑过渡。
它是明确写入的实验协议假设，不是测得的手部外力。连续图在等待间隙传递全部状态，
不会把三个单动作曲线拼接冒充连续积分结果。

## 阈值与“刚刚好”的定义

每种工况、每个 activation 来源先固定不辅助模型的背部肌肉峰值，再设置

\[
T_{limit}=0.75\max_t\tau_{back,unassisted}(t).
\]

图上的黑色水平虚线就是这个阈值，绿色区域是其下方 0-5 N·m。
考察辅助后峰值是否位于这个区域，同时检查整条曲线有无超阈值以及辅助是否产生
额外角度波谷。自然低需求段不被抬高到阈值，不能要求整条动作曲线都贴着横线。
辅助上限为 80 N·m，变化率上限为 200 N·m/s。

NMPC 使用当前模型状态、历史 activation、训练好的 MAP 预测及已声明的动作/负载协议。
实验来源额外使用固定的多步 ARX 预测器预测“实测 activation 减 MAP activation”的未来
残差。系数在全部九工况上固定，正则化通过按工况留一的 activation 预测误差选择。
在线预测不读取未来实测 EMG、未来 IK 或 CEINMS 力矩。ARX 的留一预测检查不等同于
整套 BSE/NMPC 已做独立工况验证。

控制器同时惩罚超阈值、过度辅助、运动改变和辅助能耗，并为未知输入误差保留预测裕量。
它使用两阶段受限优化：先在执行器上限和变化率允许的范围内减少预测超限，再保护运动
趋势。保留求解失败标志与回退结果；不把实际 torque 曲线裁剪到阈值。

这些减少肌肉扭矩的结果依赖模型内的人体对辅助的激活响应假设。目前 Exp 没有独立的
力矩或手部力参考，因此只能说明模型内的行为和数值一致性，不能据此证明人体真实的
腰腹分担比例或硬件辅助效果。

## 重算与核对

Python 需要 numpy、scipy、matplotlib、reportlab，以及项目现有的 C 编译器环境。
在项目根目录运行以下命令可核查当前选择，不重新标定旧的 617/627：

```powershell
python experiments/multi_subject/validate_nine_cases.py
python -m unittest discover -s tests -p 'test_nine_case_control.py' -v
python experiments/multi_subject/report_nine_cases.py
```

使用已保存的各 subject 参数重新运行完整 NMPC（结果另存，不覆盖所选版本）：

```powershell
python experiments/multi_subject/replay_selected_nine_cases.py
python experiments/multi_subject/replay_selected_nine_cases.py --subjects 240730Exp_A01
```

实际选择的 NMPC 版本由 `results/nine_case_75pct/selected_nmpc.json` 指定。
对应目录中的 `evaluation.json` 保存所有权重、预测器和逐工况结果；BSE 的训练配置
在各 subject 的 `fit_policy.json`，数据选择与处理记录在 `quality_gate.json`。
原有 528 模型、617/627 模型和之前的 Exp 结果保留；`preserved_models.json` 与
`validation.json` 记录受保护文件的哈希检查。
