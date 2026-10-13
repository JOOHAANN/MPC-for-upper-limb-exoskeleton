clc;
clear;
close all;

%% 1. 生成一个模拟的 M 形信号 (Original Signal)
x = linspace(-10, 10, 200);
% 这是一个由两个高斯函数叠加再减去一点偏移构成的 M 形
y_raw = 10 * exp(-(x + 4).^2 / 8) + 10 * exp(-(x - 4).^2 / 8);
% 让它更像 M，中间低一点
y_original = y_raw - 0.5;

%% 2. 设置你的目标参数
target_max = 25;   % 你希望峰顶变成 5
target_min = -2;  % 你希望中间的低谷变成 -2

%% 3. 寻找原信号的关键特征点 (Peak & Valley)
% 3.1 找到原信号的全局最大值 (峰顶)
P_old = max(y_original);

% 3.2 找到中间的低谷 (难点)
% 注意：min(y) 通常会找到 M 的两个脚底，而不是中间的谷。
% 方法：使用 islocalmin 寻找局部极小值
TF_min = islocalmin(y_original);
local_mins = y_original(TF_min);

if isempty(local_mins)
    error('未检测到局部低谷，信号可能不是 M 形或过于平滑。');
end

% 假设最深的那个局部极小值就是 M 的中间谷底
% (如果有多个局部低谷，这里取最小的那个作为 Anchor)
V_old = min(local_mins);

% 输出原特征以供检查
fprintf('原信号峰值: %.4f\n', P_old);
fprintf('原信号谷值: %.4f\n', V_old);

%% 4. 计算线性变换参数 (核心数学部分)
% 我们有两个方程:
% 1) a * P_old + b = target_max
% 2) a * V_old + b = target_min
% 求解 a (缩放) 和 b (偏移):

scale_factor = (target_max - target_min) / (P_old - V_old);
offset = target_max - (scale_factor * P_old);

%% 5. 应用变换
y_rescaled = y_original * scale_factor + offset;

%% 6. 验证结果
% 检查新信号的峰值和谷值是否符合要求
[new_peak, ~] = max(y_rescaled);

% 同样用局部极小值来验证谷底
TF_min_new = islocalmin(y_rescaled);
new_valley = min(y_rescaled(TF_min_new));

fprintf('--------------------------\n');
fprintf('目标最大值: %.2f -> 实际结果: %.4f\n', target_max, new_peak);
fprintf('目标最小值: %.2f -> 实际结果: %.4f\n', target_min, new_valley);

%% 7. 绘图对比
figure('Color', 'w');
subplot(2, 1, 1);
plot(x, y_original, 'LineWidth', 2, 'Color', [0 0.4 0.7]);
yline(P_old, '--', 'Old Peak');
yline(V_old, '--', 'Old Valley');
title('原始 M 形信号');
grid on;

subplot(2, 1, 2);
plot(x, y_rescaled, 'LineWidth', 2, 'Color', [0.8 0.2 0.2]);
yline(target_max, '--', ['Target Max: ' num2str(target_max)]);
yline(target_min, '--', ['Target Min: ' num2str(target_min)]);
title('加工(Rescale)后的信号');
grid on;
