# 减少过度助力，使背部力矩接近阈值下方

新结果位于 `docs/mpc_threshold_tracking/`。先看 `240627Pilot_torque_threshold.png`，或打开三页的 `threshold_tracking_results.pdf`：橙色为之前的 MPC，蓝色为本轮 MPC，黑色横虚线为阈值，绿色为阈值下方 0–5 N·m。

本轮保持同一套已训练 Hill/BSE/MAP 参数、activation、reference 和连续动作输入。改变的是控制器目标和预测余量。每个 subject 的同一套控制器覆盖两种 activation 来源及全部合格动作，角度和力矩仍由同一个 Hill 系统计算。

## 为什么之前离阈值较远

旧版主要优先避免超限，同时使用较大的预测误差储备。240627 的 400 ms 预测末端储备约为 59.4 N·m，容易提前给出偏大的助力。旧目标还缺少对“把背部肌肉力矩降低得过多”的明确惩罚，因此可能以较多助力换取较大的阈值余量。

这属于过度助力；仅仅把一个增益调小，不能同时处理预测误差、助力代价和肌肉响应滞后。本轮在完整训练段上比较参数，再冻结配置复算留出段。

## 本轮控制目标

记同一个未助力模型预测的背部力矩为 \(T_0\)，实际助力模型预测为 \(T_b\)，负载对应阈值为 \(T_{\lim}\)。需要助力时，目标取阈值下方约 3 N·m：

\[
T_{\mathrm{goal},j}=\min(T_{0,j},T_{\lim}-3).
\]

通过需求门控，只在未助力模型预测超过阈值时惩罚过度降低：

\[
d_j=\operatorname{clip}\left(\frac{T_{0,j}-T_{\lim}}{10},0,1\right),
\qquad
J_{\mathrm{low}}=5\sum_jd_j
\left[\frac{\max(T_{\mathrm{goal},j}-T_{b,j},0)}{10}\right]^2.
\]

同时提高对助力本身的惩罚：

\[
J_u=4\sum_\ell\left(\frac{u_\ell}{80}\right)^2.
\]

该项原来的权重为 0.01。调大这项权重使控制器更倾向于减少不必要的助力，而不是增加肌肉驱动力。

原训练数据得到的 99 分位误差储备记为 \(g_{99,j}\)。本轮完整训练比较了 25% 和 50% 的储备强度，两组 subject 都选择 25%，保留至少 3 N·m：

\[
g_j=3+0.25\max(g_{99,j}-3,0).
\]

预测的力矩约束仍有显式可达性松弛量：\(T_{b,j}\le T_{\lim}-g_j+s\)。先减少必要松弛量，再优化运动跟踪、过度助力和助力代价。保留 80 N·m 助力幅值上限、200 N·m/s 变化率上限及原有回正退出规则。240617 / 240627 的力矩约束范围仍分别为 100 / 400 ms。

这里的 25% 是用这批训练动作选出的经验配置，不再代表原始 99 分位储备的覆盖保证。未知动作仍可能产生新的预测误差。

原本就不需要较大背部力矩的时段会自然低于阈值；控制器不会为了使整条曲线贴线而把这些时段的肌肉力矩抬高。开始、结束及力矩上升/下降过程也保留肌肉与执行器的动态响应。

## 全部合格动作的结果

240627 的 14 段动作 × 两种 activation 来源，共 28 个组合：

| 负载 | 阈值 | 新版动作峰值范围 |
|---|---:|---:|
| 5 kg | 100 N·m | 95.6–97.5 N·m |
| 10 kg | 140 N·m | 137.3–138.0 N·m |
| 15 kg | 180 N·m | 178.0–178.5 N·m |

在未助力峰值超过阈值、确实需要助力的动作中，240627 的平均峰值余量由 **16.87 降到 2.47 N·m**，助力矩时间积分减少 **34.8%**。该积分单位为 N·m·s，不能直接当作电机能耗。

240617 的 10 段动作 × 两种来源，共 20 个组合：5 kg 峰值为 97.1–98.0 N·m，10 kg 为 137.3–137.5 N·m。部分 15 kg 动作未助力峰值本来就低于 180 N·m，因此没有为了贴线增加肌肉力矩。需要助力动作的平均峰值余量由 **5.83 降到 2.66 N·m**，助力矩时间积分减少 **21.9%**。

合计 **48 个组合无背部力矩超阈值、无新增角度谷点、无 QP 拒绝回退**；助力幅值和变化率检查通过。固定参数、零助力一致性、reference 扰动独立性、数值半步及站立恢复检查通过，七项关键 Python 回归检查通过。

图中的力矩是耦合 Hill 模型的直接输出，没有对输出力矩进行裁剪或绘图后处理。BSE/MAP 模型文件与本轮之前逐字节相同；reference、activation 和负载输入也逐项核对一致。逐动作数据在 `all_cases_before_after.csv`，训练选择在各 subject 的 `threshold_candidates.json`，复算检查在各自的 `validation.json`。

## 物理模型与验证范围

本轮减少助力后，240617 的软限位瞬时峰值约为 **65.15 N·m**，图集单独显示该力矩。240627 的助力复算未触发软限位。减少过度助力改善了阈值利用率，但没有修复或证明 BSE 模型本身的生理真实性。

这仍是使用原先假设的人体募集响应进行的离线仿真；没有新的真实助力实验 EMG 用来验证人体是否按模型减少肌肉力矩。居中 activation 滤波、数据质量限制及实验分支对离线肌骨处理的依赖见 [质量优先说明](QUALITY_FIRST.zh-CN.md)。留出动作属于开发验证，不是新增独立实验，更不是未知输入不超限的保证。

## 查看与复算

MATLAB 在项目根目录运行 `demo_hill_subjects`，默认显示最新结果；每组新增“Controller before / after”页。完整图集包含 torque、theta、腹部肌肉力矩、助力及软限位。

```matlab
demo_hill_subjects                    % View the current figures.
demo_hill_subjects('threshold-replay') % Replay the saved near-threshold MPC.
demo_hill_subjects('tune')             % Retune MPC only; BSE/MAP stay frozen.
```

也可在项目目录运行：

```powershell
python experiments/multi_subject/run_threshold.py
python experiments/multi_subject/run_threshold.py --tune
```

复算使用 `results/multi_subject_quality_v2/` 中的质量合格缓存。若尚未生成缓存，需要先运行 `demo_hill_subjects('train')`，再运行 `demo_hill_subjects('tune')`。原始数据不被覆盖，前一版结果仍保留。
