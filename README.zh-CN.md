# 项目整理及使用说明

**当前推荐入口：`demo_hill_compliant`。** 两条作用相反的弹性肌腱 Hill 肌肉共同产生 theta 和 torque，对比实验来源与 MAP 来源 activation；32 个参数参与优化，包括 M、L、阻尼和软限位，显示速度留出、历史 BSE 和各项力矩。reference 保持原项目 Scale。[运行、参数、完整方程和图表](docs/HILL_COMPLIANT.zh-CN.md)。

`demo_hill_forward` 是早期单肌肉、刚性肌腱基线；其拟合仍不理想，保留用于对照。[旧版本说明](docs/HILL_FORWARD.zh-CN.md)。

力矩诊断：运行 `demo_bse_torque_breakdown`，查看旧 BSE 边界力矩，以及新模型的惯性、重力和内部虚拟阻抗。[设置、图表和解释](docs/BSE_TORQUE_BREAKDOWN.zh-CN.md)。

历史角度平衡点 surrogate：`demo_bse_angle_priority`，其有界角度和逆动力学 torque 不符合当前“同一 Hill 系统”要求，仅保留为历史记录。[历史方法和限制](docs/BSE_ANGLE_PRIORITY.zh-CN.md)。

双输出对比：运行 `demo_bse_torque_theta`，查看原 BSE 与新降阶候选模型的
torque、theta 对 reference 的拟合。[方法、图表与限制](docs/BSE_TORQUE_THETA_COMPARISON.zh-CN.md)。

新增：运行 `demo_bse_generalization` 可查看 BSE 标定目标诊断及跨负载/速度的
力矩估计对比。[优化建议与实际结果](docs/BSE_OPTIMIZATION_REPORT.zh-CN.md)。

本项目已按“共享代码、模型变体、实验配置、数据、结果、文档”整理。
410 个有效 MATLAB 文件合并为 174 份独立实现，保留了 38 个实验配置。
同名但内容不同的函数分别保留，避免误用其他实验的模型。

原文件夹名称涉及上肢外骨骼，但现有说明和主要方程针对 BSE 背部支撑外骨骼、
L5/S1 运动，因此没有擅自更改模型的解剖学含义。

## 运行一个示例

在 MATLAB 中打开本项目根目录，然后执行：

```matlab
setup_project;
workDir = prepare_workflow('map_predictor');
cd(workDir);
[weights, XP, MEP] = MAP_model(45, 10, 'SM', 'baseline_data.mat');
plot(MEP);
```

运行原来的 MAP 对比图脚本：在同一工作目录执行 `test_MAP`。
其他主要实验入口见英文 [README](README.md)，全部配置见
[WORKFLOWS.md](docs/WORKFLOWS.md)。`list_workflows` 也可列出全部配置。

每次准备实验都会在 `results/` 创建新工作目录，复制匹配的代码和输入数据。
不要对整个项目使用 `addpath(genpath(...))`，否则同名模型函数可能冲突。
修改代码时请修改 `src/` 中的正式文件，再重新准备实验；`results/` 中是运行副本。
各标定阶段不会自动传递新参数，需要将前一阶段输出明确提供给后一阶段。

## 数据与 GitHub

- `data/processed/`：运行所用的 MAT 数据和已保存参数，合并了完全重复的文件。
- `data/raw/`：本地原始采集数据，默认不提交 Git。
- `docs/reference/`：原论文、图表、说明文档等，本地保留，默认不提交 Git。
- `results/`：运行副本与结果，默认不提交 Git。

原始项目另有完整备份。`docs/file_inventory.csv` 可查询每个原文件的去向。
GitHub 上传包不包含原始采集数据；别人可以运行已准备数据的示例，重建原始
信号时则需要按 `data/raw/README.md` 配置采集数据。

## 整理范围与验证

统一了代码缩进、空格、换行和 UTF-8 编码，补充了英文功能说明，移除了部分
超长废弃数值注释；整理了 MAP 公共函数，增加预分配及明确的 MAT 变量加载。
原来写死的个人电脑路径已改为统一的数据目录配置。

MATLAB R2023b 下通过了代码分析、全部 38 个配置的目录准备检查和 MAP 示例。
MAP 输出、低通滤波及 500 点几何模型计算与原版结果完全一致。
没有重跑所有 GA 标定或完整 NMPC 仿真；历史脚本的已知问题见
[KNOWN_LIMITATIONS.md](docs/KNOWN_LIMITATIONS.md)。
