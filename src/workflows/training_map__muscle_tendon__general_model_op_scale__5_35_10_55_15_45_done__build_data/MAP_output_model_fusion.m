% MAP_OUTPUT_MODEL_FUSION Combine predicted muscle activation with the reference movement data.
% Use prepare_workflow to select matching code and input data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main_Script.m
%%% 这里只放参数，没有乱七八糟的逻辑
%%% 35_55_45
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;

% 1. 加载基础数据
load('baseline_data.mat');
load('Combination_raw.mat');
% 假设 Reference 数据已经加载或在这里加载...

% 2. 【核心修改】定义每一段的“身份证” (长度 P, 起始点 S, 实验条件)
% ------------------------------------------------------------------
% Segment 1
tasks(1).P_length   = 6449;           % P1: 插值后的长度
tasks(1).start_idx  = 5298;           % S1: 从第几个点开始 (对应原本的 zero(1:5297))
tasks(1).speed      = 35;             % 速度
tasks(1).weight     = 5;              % 重量
tasks(1).norm_max   = 0.45;           % 全局归一化最大值
tasks(1).ref_min    = 0.08;           % 参考最小值 (可以传具体数值，也可以传数组让函数算)
tasks(1).local_norm = true;           % 是否需要那个特殊的 50:150 二次归一化
tasks(1).local_idx  = 7244:10141;     % 局部归一化的参考索引

% Segment 2
tasks(2).P_length   = 4035;           % P2
tasks(2).start_idx  = 22014;          % S2: 需要你自己算一下：上一段结束位置 + 中间间隔
tasks(2).speed      = 55;
tasks(2).weight     = 10;
tasks(2).norm_max   = 0.60;
tasks(2).ref_min    = 0.16;
tasks(2).local_norm = true;
tasks(2).local_idx  = 8029:9620;

% Segment 3
tasks(3).P_length   = 4392;           % P3
tasks(3).start_idx  = 36725;          % S3
tasks(3).speed      = 45;
tasks(3).weight     = 15;
tasks(3).norm_max   = 0.75;
tasks(3).ref_min    = 0.24;
tasks(3).local_norm = true;
tasks(3).local_idx  = 5512:7740;

% ------------------------------------------------------------------

% 3. 调用通用缝合函数
% 告诉它总长度是多少 (即 Combination 的行数)
Total_Length = size(Combination, 1);

% 【一键生成】
[stitched_signal] = run_model_stitcher(tasks, Combination(:, 2), 'baseline_data.mat');

% 4. 保存结果
Combination(:, 2) = stitched_signal;
Combination(:, 4) = Combination(:, 1);
Combination(:, 1) = 310 * stitched_signal;
save('Combination_Optimized.mat', 'Combination');

figure;
plot(Combination(:, 2));
title('Stitched Result');
