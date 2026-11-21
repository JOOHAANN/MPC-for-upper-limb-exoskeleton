# CEINMS reference 力矩核查与 A=0.1 回算

**可以重新计算并核查 reference，但目前不能把回算结果称为经过实验校准的真实力矩。** 用户确认，这两组 Exp 没有独立的逆动力学力矩、手部外力或测力台记录。本轮因此做了内部一致性核查、原求力规则辨识和 A=0.1 的敏感性回算，没有用 BSE 输出反过来调整 reference，也没有替换原始 `Torque.sto`。

## 已发现的具体差异

共检查 **54 份记录**：240725Exp 的 25 份，240730Exp 的 25 份基本记录及 4 份补录。另一个补录仍因 IK 无效跳过。这是原生肌骨导出时段的文件计数，不是 54 段完整、合格动作。

54 份导出全部满足：

1. `Torque.sto = sum(moment_arm * MuscleForces)`，最大 RMSE 约 **2.68×10⁻¹³ N·m**。因此它是肌肉模型产生的力矩估计，不是另一个独立力矩测量。
2. `MuscleForces = max(MuscleActiveForces, 0)`，逐值误差为 **0**。
3. `MusclePassiveForces` 非零，但并未进入上式导出的总力和 Torque。
4. 用记录中的肌纤维状态及 subject 参数检查，主动项等价于 **fV=1** 的力长模型，被动输出等价于**纯被动弹性项、无速度阻尼**。主动项回算最大 RMSE **7.16×10⁻⁸ N**，被动项最大 RMSE **4.88×10⁻¹³ N**。

这些是对导出数值的判定。没有当时实际编译的 DLL/源版本，不能据此断言某段历史源码一定写错；可能存在有意关闭某些项的简化版本。当前公开的 CEINMS-RT 完整实现会合并主动、被动和阻尼项，并计算力速关系，两者有明确差别。[官方求力实现，固定版本](https://github.com/CEINMS-RT/ceinmsrt-core-cpp/blob/5e53f3293ff4a5662ce575e83d2610c146de35ff/lib/NMSmodel/MTU.cpp#L296)

## 力矩不是 activation 乘一个固定比例

完整形式可写为

\[
F_i^{MT}=\left[s_iF_i^{max}\left(a_i f_L\left(\frac{l_i^m}{l_{0i}(a_i)}\right)f_V(\tilde v_i)
+f_P\left(\frac{l_i^m}{l_{0i}}\right)+d_i\tilde v_i\right)\cos\alpha_i\right]_+,
\qquad \tau_{L5/S1}=\sum_i r_iF_i^{MT},
\]

其中 \([x]_+=\max(x,0)\) 约束肌肉肌腱单元只提供拉力，**不对最后的净关节力矩截零**。负的净力矩表示屈曲方向，不单凭符号判定结果错误。

原导出表现等价于关闭被动项和速度效应：

\[
F_{i,old}^{MT}=\left[s_iF_i^{max}a_i f_L\left(\frac{l_i^m}{l_{0i}(a_i)}\right)\cos\alpha_i\right]_+.
\]

由于 \(l_{0i}(a_i)=l_{0i}[1+0.15(1-a_i)]\)，修改 A 还会改变有效最佳长度和刚性肌腱下的纤维几何。不同肌肉力臂有正负、通道的 activation 也不同，不能把旧 torque 简单乘以背部均值的下降比例。

## 实际回算的候选与图

![原 reference、A=0.1 与完整公式候选](ceinms_reference_audit/reference_torque_diagnostics.png)

上排固定使用前次已选中的两个完整动作窗口：

| 示例 | 原力矩峰值 | A=0.1、原等价规则峰值 | A=0.1、完整公式候选峰值 |
|---|---:|---:|---:|
| 240725Exp：5 kg / 45 | 321.8 N·m | **206.2 N·m** | 253.8 N·m |
| 240730Exp：10 kg / 45 | 330.4 N·m | **256.9 N·m** | 313.6 N·m |

- **灰线**：原导出的 `Torque.sto`。
- **蓝线**：保留原导出所体现的求力规则，仅把 A 改为 +0.1，同时重算激活相关长度和肌肉力。它是与原规则一致的新估计。
- **紫线**：加入力速关系、被动弹性及阻尼的完整公式诊断候选。它还出现较大变化、尖峰和负向力矩；不能据此宣称更接近真实值。
- **中排**：原文件已经记录、却没有进入原总力矩的被动项，按背部、腹部和其他组拆成**力矩，单位 N·m**。两个示例所选窗口的被动净力矩最低约 −98.5、−56.0 N·m。
- **下排**：完整公式候选的背部、腹部和其他组力矩，伸展为正、屈曲为负。

所有导出时段中的完整候选范围更宽：725 约 **−128.2～575.8 N·m**，730 约 **−160.5～514.9 N·m**。这些极值包括未被前次筛选为完整有效动作的原生时段，不能作为全工况验证成绩。它们提示旧参数与完整求力结构、纤维速度数据需要一起检查；直接补上公式项并不自动得到可信 reference。

## 回算方法与限制

EMG/activation 和肌骨导出是异步的；本轮不拟合时间平移。为保持原骨骼时刻上的相同状态，背部激活从已导出的刚性肌腱纤维几何恢复，零羽状角腹部通道通过已记录的主动肌肉力方程恢复。再用旧 A 的逆映射恢复同状态输入，改用 A=0.1 计算。这属于**基于已导出状态的敏感性回算**：因为使用了原肌纤维状态和部分主动肌肉力，前述小误差只是内部一致性，不是独立前向预测精度。

完整候选的纤维速度采用“原记录速度 + 新旧纤维长度差的 15 点因果样条导数”，与公开实时曲线缓存规则对照。没有用平滑把候选尖峰隐藏起来；新输入中的小负值仍需按输入质量规则处理。当前脚本也不是当时 CEINMS.exe 的逐指令重放。

还有一个需要保留的单位问题：记录的纤维长度和羽状角符合参数以 **degree** 解释的运行约定；原未标定 XML 的注释又写 radians。回算为了复现旧输出采用已核对的 degree 约定，**不能据此证明解剖参数的单位本来就正确**。正式标定前需要统一 XML 来源、角度单位及肌肉几何。

几何只覆盖原 `lmt.sto` 和力臂实际导出的区间。完整 EMG 中后续动作仍缺这些数据，A 的修改和本次回算不会自动补出那些参考力矩。

## 现在能做到的“校准”边界

当前可以修正处理不一致、检查单位和生理范围，并在明确假设下重估力矩。若继续做弱约束标定，可把已知负载、重力平衡、合理的肌肉参数范围和背腹共同收缩约束引入，报告假设和不确定性；这仍不是独立验证后的真实力矩。

真正 subject-specific 标定应以独立动力学目标约束肌力系数、最佳纤维长度、肌腱松弛长度等。当前没有这个目标，无法唯一识别绝对肌肉力量与力矩尺度；也不能通过让 reference 更接近 BSE 来证明 BSE 更准确。[官方标定流程](https://ceinms-docs.readthedocs.io/en/latest/Calibrate.html)

因此本轮**没有把蓝线或紫线自动替换为新的真实 reference**。需要先确定模型结构及标定依据，再统一后续 MAP/BSE/MPC 的输入和评价。

## 查看和重算

MATLAB 项目根目录运行：

```matlab
demo_ceinms_reference_audit;
```

重算核查：

```powershell
python experiments/multi_subject/audit_ceinms_reference.py
```

逐文件指标见 [records.csv](ceinms_reference_audit/records.csv)，完整参数来源、误差和用户确认的数据限制见 [audit.json](ceinms_reference_audit/audit.json)。中间候选保存在 `data/processed/ceinms_reference_audit/`，三个 STO 文件分别是原规则、加被动弹性、完整公式候选，文件头明确标记为未校准的模型估计。派生数据已从 Git 排除，原始文件未更改。

## 与原 Data 并列保存的新目录

按用户后续要求，现已将新数据打包至原 Data 目录内的两个并列目录：

- `D:\matlab\Journal Paper of Upper-Limb Exoskeleton\YOUHAN LI code and project description\Data\240725Exp_A01`：25 份记录。
- `D:\matlab\Journal Paper of Upper-Limb Exoskeleton\YOUHAN LI code and project description\Data\240730Exp_A01`：25 份基本记录和 4 份补录。

新目录 `CalibrationData/` 外层的 `Activations.sto`、肌肉力/状态及 `Torque.sto` 已按 A=0.1 重算；**默认 Torque 采用原等价主动项规则，仍是未独立标定的模型估计**。原结果保存在每试次的 `Original_CEINMS/`，原 XML 在 `Original_configuration/`。完整公式与加被动弹性候选另用文件名标记。

新旧目录彼此独立，原文件未覆盖；派生标记避免它们被当成新增受试者重复统计。各新目录的 `README.zh-CN.md` 和 `derived_dataset.json` 记录内容、时钟、来源和逐文件哈希。复建打包目录可运行 `python experiments/multi_subject/package_exp_datasets.py`。
