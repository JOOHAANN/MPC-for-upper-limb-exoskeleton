clc;
clear;
close all;

% 数据点
x = 1:25;  % x 是一个 1x25 的向量
y = rand(25, 6, 4);  % y 是一个 25x6x4 的随机矩阵，用实际数据替换

% 获取 y 的维度
[dim1, dim2, dim3] = size(y);

% 初始化存储拟合系数的矩阵
coefficients = zeros(3, dim2, dim3);

% 对每个 6x4 矩阵的每个元素进行二次多项式拟合
for i = 1:dim2
    for j = 1:dim3
        % 获取 y 中每个位置的时间序列
        y_series = squeeze(y(:, i, j));

        % 进行二次多项式拟合
        p = polyfit(x, y_series, 2);

        % 存储拟合系数
        coefficients(:, i, j) = p;
    end
end

% 打印结果
disp('拟合得到的系数 (每个位置对应一个响应变量):');
disp(coefficients);

% 可视化拟合结果（例如对某个位置 (1,1) 的拟合结果进行可视化）
figure;
y_series = squeeze(y(:, 1, 1));
p = coefficients(:, 1, 1);
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(p, x_fit);
plot(x, y_series, 'bo', 'MarkerFaceColor', 'b'); % 原始数据点
hold on;
plot(x_fit, y_fit, 'r-'); % 拟合曲线
xlabel('x');
ylabel('y');
title('位置 (1,1) 的二次多项式拟合');
legend('数据点', '拟合曲线');
grid on;
hold off;
