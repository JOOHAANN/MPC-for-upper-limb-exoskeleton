clc;
clear all;

% 生成一些示例数据
x = 1:0.1:30;
y = 3 * exp(-((x - 5) / 1.5).^2) + 2 * exp(-((x - 15) / 2).^2) + randn(size(x)) * 0.2; % 三个高斯峰和一些噪声

% 绘制原始数据
figure;
scatter(x, y);
hold on;

% 拟合数据
% 使用内置的 gauss3 模型
ft = fittype('gauss3');
startPoints = [3 5 1.5 2 15 2 1.5 25 1]; % 初始猜测值
fitResult = fit(x', y', ft, 'Start', startPoints);

coeffs = coeffvalues(fitResult);

% 显示提取的参数
a1 = coeffs(1);
b1 = coeffs(2);
c1 = coeffs(3);
a2 = coeffs(4);
b2 = coeffs(5);
c2 = coeffs(6);
a3 = coeffs(7);
b3 = coeffs(8);
c3 = coeffs(9);

% 绘制拟合结果
plot(fitResult, x, y);
legend('Data', 'Gaussian Fit');
hold off;

% 显示拟合参数
disp('Fitted parameters:');
disp(fitResult);
