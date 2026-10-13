clc;
clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');
% DEMO_THETA_LS_TWOPEAKS
% 目标：
%  1) 构造一个“倒着的麦当劳标志(W形)”的真实角度 theta_true，范围[0.5,1.6]。
%  2) 用系统模型生成扭矩 tau_model，并归一化为 [20,120] N·m 的“两座山峰”扭矩 tau（输入 model_torque）。
%  3) 只给 tau，采用最小二乘（lsqnonlin）反演整段 theta（未知 dtheta, ddtheta 用差分近似）。
%  4) 不使用卡尔曼滤波，仅最小二乘 + 二阶差分平滑 + 角度上下界。
%
% 依赖：Optimization Toolbox（lsqnonlin）
load('Combination.mat');
%% ------------------ 0) 基本设置 ------------------
rng(0);                                % 固定随机种子，便于复现
N  = 400;                              % 采样点数
dt = 0.01;                             % 采样周期(s)
t  = (0:N - 1)' * dt;                    % 时间向量

% ---- 系统参数（按你提供的数值） ----
c   = 10;                              % 粘性阻尼系数
m   = 5.0;                             % 质量(项1)
M   = 22.9819;                         % 质量(项2)
g   = 9.81;                            % 重力加速度
L1  = 0.345;                           % 几何长度1
L2  = 0.756835;                        % 几何长度2
phi = 0.950692;                        % 夹角（弧度）
I   = m * L2^2 + (4 * M * L1^2) / 3;           % 转动惯量（给定公式）

% ---- 角度限位 ----
theta_lb = 0.4;                        % 角度下界
theta_ub = 1.6;                        % 角度上界

%% ------------------ 1) 构造“真实”W形角度 θ_true ------------------
% 思路：用两个高斯“下陷”叠加到一个高位基线，得到 W 形（两谷、三高）。
%       再截断到 [0.5, 1.6] 范围。
t1 = t(round(N * 0.30));                 % 第一个谷中心（≈30%处）
t2 = t(round(N * 0.70));                 % 第二个谷中心（≈70%处）
w  = 0.18;                             % 谷的宽度(标准差)，越小下陷越尖锐

g1 = exp(-0.5 * ((t - t1) / w).^2);        % 高斯1（数值范围 0..1）
g2 = exp(-0.5 * ((t - t2) / w).^2);        % 高斯2（数值范围 0..1）
r  = g1 + g2;                          % 两个下陷的叠加
r  = r - min(r);                       % 平移到 >=0
r  = r / (max(r) + eps);               % 归一化到 [0,1]

theta_peak   = 1.55;                   % 希望的“高位”角度（接近上界但预留余量）
theta_valley = 0.55;                   % 希望的谷底角度（接近下界但预留余量）
theta_true   = theta_peak - (theta_peak - theta_valley) * r; % r=1时到谷底，r=0时到高位
theta_true   = min(max(theta_true, theta_lb), theta_ub);     % 保险：截断到范围内

% 计算真值的速度/加速度（仅用于生成扭矩；估计时不会用）
dtheta_true  = dtheta_from_theta(theta_true, dt); % 与 theta 等长的一阶差分
D2           = D2_operator(N, dt);                % N×N 二阶差分矩阵
ddtheta_true = D2 * theta_true;                   % 二阶差分近似加速度

%% ------------------ 2) 先造 “两峰” 目标扭矩 tau_des ，再解出 Tac 以保持一致 ------------------
% 目标形状：两个高斯峰，峰值 120 N·m，谷底 20 N·m（平滑）

% 关键：令 (theta_true, tau_des) 对模型是可行的，反推出 Tac

% tau = I*ddtheta_true + c*dtheta_true ...
%     + m*g*L2*cos(theta_true - phi) + M*g*L1*cos(theta_true);

tau = Combination(7966:13261, 1);
tau_x = (1:length(tau))';
tau_x_new = (linspace(1, length(tau), 400))';
tau = interp1(tau_x, tau, tau_x_new, 'linear');
% 之后优化阶段，用 ‘tau = tau_des’、‘Tac = 上面这个Tac’
% tau = tau_des;

%% ------------------ 3) 最小二乘反演 θ（只用 tau，未知 dθ, ddθ 用差分） ------------------
% 残差 r(th) = [ w_tau*(model(th) - tau);
%                w_smooth*(D2*th) ]
% 其中 model(th) = I*(D2*th) + c*dtheta(th) - Tac + m*g*L2*cos(th - phi) + M*g*L1*cos(th)

% 速度近似函数（中心差分，首尾用单边差分）
dth_fun = @(th) dtheta_from_theta(th, dt);

% 建立模型函数：用 th 估计扭矩
% model_fun = @(th) I*(D2*th) + c*(dth_fun(th)) - Tac ...
%                  + m*g*L2*cos(th - phi) + M*g*L1*cos(th);
model_fun = @(th) I * (D2 * th) + c * (dth_fun(th)) ...
                 + m * g * L2 * cos(th - phi) + M * g * L1 * cos(th);

% 权重设置（可根据数据量级微调）
w_tau    = 1.0;                        % 扭矩一致性权重
w_smooth = 1e-2;                       % 二阶平滑权重（抑制抖动）

% 拼接残差向量：长度 N（扭矩项） + N（二阶平滑项）
residual = @(th) [w_tau   * (model_fun(th) - tau)
                  w_smooth * (D2 * th)];

% 初值：取在界内的常值（中点），对非线性收敛友好
theta0 = ((theta_lb + theta_ub) / 2) * ones(N, 1);

% 上下界：角度限位
lb = theta_lb * ones(N, 1);
ub = theta_ub * ones(N, 1);

% lsqnonlin 选项
opts = optimoptions('lsqnonlin', ...
                    'Display', 'iter', ...                   % 显示迭代信息
                    'MaxFunctionEvaluations', 5e5, ...
                    'MaxIterations', 1000, ...
                    'SpecifyObjectiveGradient', false);      % 不提供显式雅可比

% 开始优化
fprintf('开始最小二乘反演 theta（角度限位 [%.3f, %.3f]）。\n', theta_lb, theta_ub);
theta_est = lsqnonlin(residual, theta0, lb, ub, opts);

%% ------------------ 4) 结果评估与可视化 ------------------
% 根据估计的 theta 计算速度/加速度（用于展示）
dtheta_est  = dth_fun(theta_est);
ddtheta_est = D2 * theta_est;

% 用估计的 theta 回代模型，得到拟合扭矩
tau_est = model_fun(theta_est);

% 角度/扭矩 RMSE（演示：theta_true 仅用于参考）
rmse_theta = sqrt(mean((theta_est - theta_true).^2));
rmse_tau   = sqrt(mean((tau_est   - tau).^2));
fprintf('角度 RMSE = %.4f rad, 扭矩 RMSE = %.4f N·m\n', rmse_theta, rmse_tau);

% 画图窗口
figure('WindowStyle', 'normal', 'Color', 'w', 'Position', [100 100 1100 780]);

% (1) 角度：真值 vs 估计（含限位线）
subplot(3, 1, 1);
plot(t, theta_true, 'k-', 'LineWidth', 1.5);
hold on;
plot(t, theta_est,  'r--', 'LineWidth', 1.5);
yline(theta_lb, 'b:', 'lb');
yline(theta_ub, 'b:', 'ub');
xlabel('Time (s)');
ylabel('\theta (rad)');
legend('\theta true (W-shape)', '\theta estimated', 'Location', 'best');
title('Angle: True (W) vs Estimated');

% (2) 扭矩：两峰观测 vs 回代模型
subplot(3, 1, 2);
plot(t, tau,     '-', 'LineWidth', 1.0);
hold on;
plot(t, tau_est, '--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Torque (N·m)');
legend('measured \tau (two peaks)', 'model(\theta_{est})', 'Location', 'best');
title('Torque fit (Two-peak input)');

% (3) 估计的速度（供参考）
subplot(3, 1, 3);
plot(t, dtheta_est,  'r-', 'LineWidth', 1.0);
hold on;
plot(t, dtheta_true, 'k--', 'LineWidth', 1.0);
xlabel('Time (s)');
ylabel('d\theta (rad/s)');
legend('d\theta estimated', 'd\theta true (ref)', 'Location', 'best');
title('Angular velocity (reference only)');

sgtitle('Least-Squares Inversion of \theta with Two-Peak Torque Input');

% ====== 主函数结束 ======

% ===================== 辅助函数 =====================

function dthN = dtheta_from_theta(th, dt)
    % dtheta_from_theta
    % 说明：把 theta (N×1) 用差分转成同长度的 dtheta (N×1)。
    N = numel(th);                    % 序列长度
    dthN = zeros(N, 1);                % 预分配
    dthN(1)     = (th(2)   - th(1))   / dt;         % 首点：前向差分
    dthN(2:N - 1) = (th(3:N) - th(1:N - 2)) / (2 * dt);   % 中间：中心差分
    dthN(N)     = (th(N)   - th(N - 1)) / dt;         % 末点：后向差分
end

function D2 = D2_operator(N, dt)
    % D2_operator
    % 说明：构造 N×N 的二阶差分稀疏矩阵，近似 ddtheta。
    e  = ones(N, 1);                                   % 单位向量
    D2 = spdiags([e -2 * e e], -1:1, N, N) / (dt^2);    % 三对角稀疏矩阵
    % （边界采用相同模板的简单处理，若需更精细可自行改进）
end
