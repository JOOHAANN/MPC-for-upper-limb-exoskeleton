# MAP 使用 excitation；BSE/NMPC 只做 5、10、15 kg

**725 BSE 后续版本：** 现已扩展到有 torque 支持的 10 个负载/速度工况，排除 5 个零支持工况，见 [跨速度方法及结果](EXP_SPEED_GRID.zh-CN.md)。本文保留此前三动作标定和 NMPC，控制结果仍对应此前模型。

按用户进一步澄清，MAP 与 BSE 的数据准入分开：**MAP 用 muscle excitation，BSE/NMPC 只需要 5、10、15 kg 三个负载的代表动作。** 没有要求两者共享全部工况，也不再用一个固定的“至少五段合格动作”条件阻止三负载 BSE 标定。

617、627 和已完成的 730 Exp 模型保留并跳过。本轮补做 725 Exp 的三负载 BSE/NMPC。730 之前额外回放过的 12.5 kg 文件保留在历史目录，本报告只显示用户要求的三种负载。

## MAP 训练独立于肌骨导出

两组现有 MAP 已训练，参数本轮不重训；例如 **240725Exp 的 5kg/25 已经参与 MAP 训练**，没有因其配套力矩提前结束而剔除。

新训练入口 `excitation_only_map.py` 对每份试次只读取 `emgFilt.sto`：检查 excitation 的时间戳、断档、完整 M 峰和平滑度，按既定 `A=0.1` 映射产生 activation 标签，学习四个非负协同、三个 Gaussian 基元及负载/节律到权重的映射。它不需要 IK、LMT、MuscleForces 或 Torque 导出。训练得到的参数再用于 activation 预测；BSE 的缺文件情况不决定 MAP 是否可以训练。

已在 725 的实际 excitation 文件上单独复核，25 个工况都有可选输入；另有只含 excitation 文件的自动回归检查。已有预测器的参数和输入窗口指纹不变，新入口没有为了更换准入方式而重新拟合已完成模型。旧缓存的背肌融合权重保留；全新 excitation-only 训练使用六通道等权融合，不需要从肌骨几何计算融合权重。

## 725 的三负载 BSE 标定

| 负载 | 标定动作 | 角度主要谷点 | 实际 CEINMS torque 覆盖 |
|---|---|---:|---:|
| 5 kg | 5kg/45 | 2 | 100% |
| 10 kg | 10kg/45 | 2 | 100% |
| 15 kg | 15kg/55 | 2 | 约 71% |

这三段都保留完整 excitation 和实测 IK，所以能够计算并评价整段角度趋势。15 kg 缺失的力矩时段明确标成缺失，从 torque 拟合和 torque RMSE 中排除，没有用 BSE 输出、端点重复值或其他试次曲线补成 reference。

角度保留原项目的 Scale。已有原生 LMT 支持区间的 reference 与旧处理完全相同；IK 比肌骨导出更长的部分沿用同一缩放系数及偏置，不针对所选动作重新做 min/max，也不裁剪角度。因此 15 kg 全动作 reference 约为 0.386–1.523 rad，这是原记录延伸后的同一比例，不是手工调出来的目标范围。

原来的 angle/CEINMS-torque 相关性筛选改作警示：CEINMS torque 本身是尚无独立测力验证的模型估计，不能仅因它与 IK 相关性弱就否定完整 excitation 或实测角度。原生导出仍检查时间戳和 `Torque = Σ(force × moment arm)`；硬件延迟没有另行实测。

![三段标定输入及缺失力矩](exp_three_loads/three_anchor_inputs.png)

BSE 沿用同一耦合 Hill、肌腱和神经调节架构，所有 32 个物理参数及 10 个神经参数可参与优化；M、L、负载几何和 soft boundary 在实验/MAP 两分支间共享。每分支只有一组固定参数，覆盖三种负载。实际 θ 和 torque 来自同一个模型，没有用 reference 直接驱动输出角度。

两轮标定后，725 实验 activation 分支三个标定动作的角度趋势相关性约为 **0.913、0.865、0.821**；MAP 分支约为 **0.844、0.878、0.724**。六条标定回放都保留两个主要谷点。CEINMS 有效区间的 torque RMSE 约 57–89 N·m，说明力矩一致性仍有明显误差，不能称作绝对力矩已校准。

![三种负载的标定动作对比](exp_three_loads/240725Exp_A01_calibration.png)

图中各列为 5、10、15 kg，上半部分是实验 activation，下半部分是训练 MAP 生成的 activation；分别比较 θ 与净肌肉 torque。红虚线 torque 只画实际存在的 CEINMS 估计，不补缺失尾段。这张图展示标定动作，不当作独立验证。

## 额外回放与结果含义

另留出 5kg/35、10kg/25、15kg/65 三段做 BSE 额外回放；它们没有参与本轮 BSE 参数优化。部分记录此前已经用于冻结 MAP 的训练，所以这里是 **BSE 参数留出检查**，不是整条 MAP→BSE→NMPC 流程的独立验证。并非必须有这些额外动作才能启动三负载标定。

额外试次的趋势拟合弱于标定动作，尤其 10kg/25、15kg/65；有的输出出现额外谷点。这部分会展示，不因三段标定表现较好就宣称泛化已解决。只有三个标定动作，也不足以唯一辨识全部生理参数，仍依赖参数范围、先验和初始化。

NMPC 在冻结的 BSE/MAP 上仅用这三个标定动作选参数，保留 0–80 N·m、200 N·m/s 的执行器限制。背部设计目标仍为 100/140/180 N·m。结果图同时给出实验/MAP activation、θ、净肌肉 torque、背部/腹部 torque、阈值虚线、助力及 soft-boundary torque。具体峰值、超限及数值检查见 [最新状态和指标](exp_three_loads/completion_status.json)。

725 的 6 段动作 × 两种 activation 共 12 组回放均未超出背部目标。三个标定动作的 MPC 峰值约为实验 **91.7 / 133.3 / 176.1 N·m**、MAP **91.2 / 131.8 / 172.7 N·m**。实验分支 15 kg 标定段的 soft-boundary 瞬时峰值约 **1.71 N·m**，其余这些回放为零。步长减半的最大角度差约 0.0033 rad、净力矩差约 0.53 N·m，四个站立扰动恢复检查通过；改变 reference 后控制器输出差为零。

**无助力 BSE 仍有一项明显限制：** 实验 activation 的 15 kg 标定动作 soft-boundary 峰值约 **74.17 N·m**，额外 15 kg 动作约 **69.37 N·m**；MAP 分支这些动作为零。这个恢复力矩不是肌肉 torque，不能把它计作肌肉生理有效性的证明。连续图最下面同时画无助力边界力矩（点线）与 MPC 后边界力矩（实线），避免只看到控制后很小的数值而误判原 BSE 已没有边界依赖。这轮已完成三负载训练和控制回放，绝对力矩、生理参数与跨节律泛化仍未解决。

## 使用

```matlab
demo_exp_subjects                 % 最新三负载图与五负载 MAP 图
demo_exp_subjects('status')       % 四组 subject 的分阶段状态
demo_exp_subjects('run')          % 跳过完成的模型，只补齐缺少的阶段
p = predict_exp_activation('240725Exp_A01',10,9);
plot(p.timeSeconds,p.backActivation);
```

最后一例是输入 10 kg、完整动作频率 9 BPM 来生成 activation。文件名的速度代码与这个完整动作 BPM 不等同。

训练和 NMPC 仍由 Python/C 执行，MATLAB 提供调用及查看器，未运行旧 Simulink 模型。原始数据目录及已完成模型不覆盖；这次缺失力矩的 mask 只属于新的准备数据。原模型结果保留在旧目录，新增中间结果在 `results/exp_three_loads*`，展示结果在 `docs/exp_three_loads/`。
