# 725 Exp：有力矩参考的跨速度 BSE 标定

本轮只处理 `240725Exp_A01` 的 BSE，负载为 **5、10、15 kg**，先检查速度代码 **25、35、45、55、65**，再按用户要求排除所选完整动作中完全没有 CEINMS torque 参考的工况。**保留 10 个工况，优先保证其中 35–55 的 6 个工况。**

| 负载 | 保留的速度代码 | 完全排除的速度代码 |
|---|---|---|
| 5 kg | 25、35、45、65 | 55 |
| 10 kg | 25、45、55 | 35、65 |
| 15 kg | 45、55、65 | 25、35 |

排除项既不进入角度训练，也不进入力矩训练、指标均值或成功工况计数。图中保留灰色位置，仅说明数据为何没有纳入。部分有 torque 的动作仍保留完整的 activation/IK；缺少 torque 的样本保留 NaN，并用显式支持掩码从 torque loss 和 RMSE 排除。没有端点补齐，也没有用模型输出构造 reference。

这里“没有参考”指所选完整动作时间范围与 CEINMS 导出没有重叠，并非说整个试次没有 `Torque.sto` 文件。完整动作的选择沿用事先审查的结果；例如 15kg/45 的早期候选存在超过 250 ms 的输入断档，因此没有为了增加 torque 覆盖而换成那段数据。

## 模型与优化

沿用既有的神经调节与弹性肌腱 Hill 联立模型。每个 activation 来源仅有一组固定的 **32 个物理参数 + 10 个神经参数**，覆盖这 10 个工况，不按速度或负载切换参数。所有参数均可参与有界优化，包括 M、身体和负载力臂；M、身体/负载几何及 soft boundary 在实验与 MAP 来源之间共享。

实验分支使用既有 A=+0.1 activation。MAP 分支使用已训练预测器生成的 activation，不直接拿实验背肌 activation 充当预测值；辅助腹肌输入生成器由保留工况的 activation 标签拟合，不使用角度/力矩作为其输入标签。已有 MAP 预测器保持冻结，故这次 BSE 工况排除不改变之前 MAP 的训练数据或参数。

神经层通过肌肉 excitation 施加调节，关节角度与肌肉 torque 由同一个 Hill 系统产生。运行时不将 reference theta、reference torque 或其未来样本输入神经控制器。这里包含神经运动意图与本体反馈假设，不是从原始 EMG 单独推断绝对角度的无反馈 Hill 模型。

35、45、55 的各工况采用 2 倍损失权重，25、65 为 1 倍。每个动作按样本数归一化，避免慢速动作只因为采样点多就主导目标。训练联合考虑角度形状、宽松谷点深度、可用区间 torque、腹背共同收缩、等待段恢复和 soft-boundary torque。角度不裁剪，soft boundary 是弹性恢复项，其扭矩单独输出，不并入肌肉 torque。

后续阶段对较大的标准化角度形状误差增加平滑四次尾项，使最差工况得到更大权重；另提高宽松的角度值/谷点深度惩罚，避免只得到高相关性、弯腰却偏浅的曲线。没有逐工况修改目标或参数。实现仍是 Python 有界最小二乘训练、C 联立求解和 MATLAB 查看器，本轮没有运行 Simulink。

相对之前三动作标定，将 soft-boundary 残差系数从 0.25 提高到 2.0，以减少用边界弹力代替肌肉控制的倾向。这是优化权重，不是把边界硬锁死。原有角度 Scale 完整保留；已有肌骨导出支持区间逐样本与旧 reference 一致，外延 IK 沿用同一比例和偏置，不重新归一化选定动作。

后续标定允许运动段神经反馈占比在 0.05–1 范围优化，替代原先最低 0.7 的搜索范围；等待段仍保留原来的肌肉激励平衡机制。因部分神经参数到达旧搜索范围上限，也检查了更宽的因果 EMG 导数增益和屈曲请求饱和范围。导数增益只使用过去输入，不是读取未来 EMG。实际采用的范围和优化阶段记录在 [参数策略](exp_speed_grid/fit_policy.json)。这些上界是参数搜索范围，不等于图中的实际腹肌 torque。

## 结果与解释

定量前后对比见 [结果汇总](exp_speed_grid/summary.zh-CN.md)，完整逐工况指标见 [评价数据](exp_speed_grid/evaluation.json)，数值与 reference 独立性检查见 [验证数据](exp_speed_grid/validation.json)。图上的 `r` 仅在真实动作段计算；0.7 只是便于浏览的工程趋势指标，不能单独证明 W 形正确或生理参数真实。还需同时检查谷点数、弯腰幅度、腹背 torque 和边界出力。曲线包括动作前后站立等待，竖点线标记动作起止；等待期没有实测 reference 的部分不补作实测值。边界峰值同时统计完整序列，避免遗漏动作后尖峰。

本轮选定参数的全部保留工况平均趋势相关性：实验 **0.529→0.752**，MAP **0.529→0.760**；35–55 两分支均为 **5/6** 个工况达到 0.7。**15kg/45 仍是主要未解决项**，实验/ MAP 分别约 0.277/0.615，不能宣称全部 35–55 已满足；10kg/25 还会出现额外谷点，15kg/65 的时序也有偏差。这些问题完整保留，不因其难拟合而删掉工况。未证明其数据错误，目前能确认的是这套共享方程及固定参数对这些动作仍有局限。

角度趋势改善伴随 torque 一致性代价：有效区间的逐工况 torque RMSE 均值从实验/MAP **63.5/66.1** 增至 **70.3/72.6 N·m**。因此不能称力矩也全面改善。全序列软限位峰值为实验 **7.05 N·m**、MAP **0**；腹/背平均 torque 比在各工况中最大约 **10.7%/11.8%**，背部仍为主要出力来源。这些比例是模型结果，尚无独立力矩测量验证。

步长减半最大差为角度约 **0.0076 rad**、净 torque **1.26 N·m**；四个站立扰动恢复检查通过，改变评价 reference 后仿真输出差为零。模型与已完成 MAP/其他 subject 的保护校验、各次标定策略和选择记录分别见 [保护清单](exp_speed_grid/preserved_models.json)、[本轮参数选择](exp_speed_grid/model_selection.json)。

![实验 activation 的角度](exp_speed_grid/experimental_theta.png)

![训练 MAP activation 的角度](exp_speed_grid/map_theta.png)

![实验 activation 的净肌肉 torque](exp_speed_grid/experimental_torque.png)

![MAP activation 的净肌肉 torque](exp_speed_grid/map_torque.png)

![实验 activation 的背部和腹部 torque](exp_speed_grid/experimental_muscles.png)

![MAP activation 的背部和腹部 torque](exp_speed_grid/map_muscles.png)

![独立列出的 soft-boundary torque](exp_speed_grid/soft_boundary.png)

CEINMS torque 仍是没有独立测力验证的模型估计，不称作实测 ground truth。全部保留工况都参与了这轮标定，所以这些结果不是独立泛化测试，也不能由它们唯一辨识全部生理参数。背部/腹部比例限制属于工程先验，尚无独立肌肉力矩来验证绝对数值。

## 查看和复现

在 MATLAB 中打开本项目根目录：

```matlab
demo_725_bse_speed_grid          % 查看保存结果，不重新训练
demo_725_bse_speed_grid('run')   % 模型缺失时训练；已完成时跳过训练
```

Python 入口为 `experiments/multi_subject/run_exp_speed_grid.py --fit`。中间结果在 `results/exp_speed_grid_torque/`，可携带模型在 `data/processed/exp_speed_grid/`，展示图和 PDF 在 `docs/exp_speed_grid/`。`--refine` 只在明确需要继续训练时使用，并且校验准备数据指纹后才恢复 checkpoint。

本轮没有重新调 NMPC。`demo_exp_subjects` 中 725 的三负载 NMPC 图仍对应上一版 BSE；不能将它解释为本轮跨速度模型的控制结果。617、627、730 及已有 MAP/MPC 缓存保持原样。
